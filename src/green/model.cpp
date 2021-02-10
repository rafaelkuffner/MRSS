
#include "model.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <chrono>
#include <mutex>
#include <thread>
#include <algorithm>
#include <memory>
#include <charconv>
#include <string_view>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/euler_angles.hpp>

#include <cgu/shader.hpp>

#include <imgui.h>

#include "imguiex.hpp"
#include "main.hpp"
#include "curvature.hpp"

#include "model.glsl.hpp"

namespace {

	template <size_t NBins>
	float histogram_entropy(const std::array<float, NBins> &hist) {
		float atot = 0;
		for (auto &a : hist) atot += a;
		const float iatot = 1.f / atot;
		float s = 0;
		const float nilog2 = -1.f / log(2.f);
		for (auto &a : hist) {
			// limit of p*log(p) as p tends to 0 is 0
			// so we can safely discard anything very small
			const float p = a * iatot;
			if (!(p > 1e-6f)) continue;
			s += p * log(p) * nilog2;
		}
		return s;
	}

	template <size_t NBins>
	void print_histogram(const std::array<float, NBins> &hist) {
		using namespace std;
		float atot = 0;
		for (auto &a : hist) atot += a;
		const float iatot = 1.f / atot;
		for (int i = 0; i < NBins; i++) {
			cout << setw(3) << i << " | ";
			const string_view maxbar = "================================================================================";
			int c = maxbar.size() * hist[i] * iatot;
			cout << maxbar.substr(0, c) << endl;
		}
	}

	// by area, normalized
	template <size_t NBins = 256, typename URNG>
	inline std::array<float, NBins> maxdon_histogram(
		const green::PolyMesh &mesh,
		OpenMesh::VPropHandleT<float> vertexAreasProperty,
		float normalPower,
		float maxval,
		URNG &&rand,
		const std::vector<int> &shuffled_samples
	) {
		using namespace green;
		constexpr size_t maxsamples = 100000;
		size_t start = 0;
		if (mesh.n_vertices() > maxsamples) start = std::uniform_int_distribution<size_t>(0, mesh.n_vertices() - maxsamples)(rand);
		assert(mesh.has_vertex_normals());
		const float hist_irange = (NBins - 1) / maxval;
		std::array<float, NBins> hist{};
		float atot = 0;
		// note: this sampling is not area-uniform
		for (size_t i = start; i < std::min(start + maxsamples, mesh.n_vertices()); i++) {
			PolyMesh::VertexHandle v(shuffled_samples[i]);
			const float area = mesh.property(vertexAreasProperty, v);
			const float don = std::min(maxdon(mesh, v, normalPower), maxval);
			const auto bin = intptr_t(hist_irange * don);
			hist[bin] += area;
			atot += area;
		}
		const float iatot = 1.f / atot;
		for (auto &a : hist) {
			a *= iatot;
		}
		return hist;
	}

}

namespace green {

	ModelBase::ModelBase(const std::filesystem::path &fpath) {
		m_mesh.request_face_normals();
		m_mesh.request_vertex_normals();
		m_mesh.request_vertex_colors();
		// face status not really required here
		// NOTE face status currently breaks decimation
		//m_trimesh.request_face_status();
		// load custom "quality" property and vertex colors if they exist
		// FIXME what else do we need to ask for and preserve on export? texcoords?
		OpenMesh::IO::Options readOptions = OpenMesh::IO::Options::Custom | OpenMesh::IO::Options::VertexColor;

		std::cerr << "Loading model " << fpath.u8string() << std::endl;
		if (!std::filesystem::is_regular_file(fpath)) {
			throw std::runtime_error("file does not exist");
		}

		if (OpenMesh::IO::read_mesh(m_mesh, fpath, readOptions)) {
			std::cerr << "Loaded" << std::endl;
		} else {
			std::cerr << "Failed" << std::endl;
			throw std::runtime_error("failed to load model");
		}

		// TODO necessary? optional?
		// NOTE currently done by assimp too
		m_mesh.triangulate();
		
		// copy original vertex colors
		// (because we need to be able to overwrite the actual vertex colors during export)
		if (readOptions.check(OpenMesh::IO::Options::VertexColor)) {
			std::cout << "Found vertex colors" << std::endl;
			m_mesh.add_property(m_prop_vcolor_original);
			for (auto vIt = m_mesh.vertices_begin(); vIt != m_mesh.vertices_end(); ++vIt) {
				m_mesh.property(m_prop_vcolor_original, *vIt) = m_mesh.color(*vIt);
			}
		}

		// check for 'raw' saliency properties
		std::vector<std::pair<saliency_prop_t, model_saliency_data>> sprops;
		for (auto it = m_mesh.vprops_begin(); it != m_mesh.vprops_end(); ++it) {
			// TODO any way to do this other than dynamic cast?
			auto prop = dynamic_cast<OpenMesh::PropertyT<float> *>(*it);
			if (!prop) continue;
			if (prop->name().substr(0, 7) == "quality") {
				std::cout << "Found saliency property with internal name " << prop->name() << std::endl;
				auto ph0 = saliency_prop_t{int(it - m_mesh.vprops_begin())};
				if (&m_mesh.property(ph0) != prop) abort();
				model_saliency_data sd;
				sd.filename = fpath.filename().u8string();
				auto &pname = prop->name();
				if (auto i = pname.find('['); i != std::string::npos && pname.back() == ']') {
					// property name contains metadata
					sd.dispname = pname.substr(7, i - 7);
					auto meta = std::string_view(pname).substr(i + 1, pname.size() - i - 2);
					std::string_view metamisc = "";
					std::string_view metauparams = meta;
					if (auto j = meta.find(';'); j != std::string::npos) {
						// misc metadata before saliency user params
						metamisc = meta.substr(0, j);
						metauparams = meta.substr(j + 1);
					}
					sd.decimated = metamisc.find('D') != std::string::npos;
					sd.uparams_known = sd.uparams.parse(metauparams);
				} else {
					// no metadata
					sd.dispname = pname.substr(7);
				}
				if (auto i = sd.dispname.find('@'); i != std::string::npos) {
					// name has file-local unique index
					sd.dispname = sd.dispname.substr(0, i);
				}
				// re-persist by default
				sd.persistent = true;
				sd.should_export = true;
				sprops.emplace_back(ph0, std::move(sd));
			}
		}

		// move 'raw' saliency to unnamed properties
		// need 2 passes because we cant modify the properties while iterating them
		for (auto &p : sprops) {
			saliency_prop_t ph;
			m_mesh.add_property(ph);
			m_mesh.property(ph).data_vector() = std::move(m_mesh.property(p.first).data_vector());
			m_mesh.remove_property(p.first);
			auto sd = std::move(p.second);
			sd.prop_saliency = ph;
			m_saliency.push_back(std::move(sd));
		}

		// move original vids to unnamed property
		if (OpenMesh::VPropHandleT<int> pvid; m_mesh.get_property_handle(pvid, "original_vid")) {
			std::cout << "Found original vertex ids property" << std::endl;
			m_mesh.add_property(m_prop_vid_original);
			m_mesh.property(m_prop_vid_original).data_vector() = std::move(m_mesh.property(pvid).data_vector());
			m_mesh.remove_property(pvid);
		}

	}

	std::vector<model_saliency_data>::const_iterator ModelBase::find_saliency(std::string_view name) const {
		if (name.size() > 0 && name.front() == '@') {
			// parse index
			int i = 0;
			auto r = std::from_chars(&*name.begin() + 1, &*name.end(), i);
			if (r.ec != std::errc{}) return m_saliency.end();
			if (i < 0) i = int(m_saliency.size()) + i;
			if (i < 0) return m_saliency.begin();
			if (i >= m_saliency.size()) return m_saliency.end();
			return m_saliency.begin() + i;
		} else {
			// search by name
			for (auto it = m_saliency.begin(); it != m_saliency.end(); ++it) {
				if (it->dispname == name) return it;
			}
			return m_saliency.end();
		}
	}

	Model::Model(ModelBase &&base) : ModelBase(std::move(base)) {
		
		// calculate normals always
		std::cout << "Computing vertex normals" << std::endl;
		m_mesh.update_face_normals();
		m_mesh.update_vertex_normals();

		std::cout << "Computing bounding box" << std::endl;
		for (auto vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit) {
			m_bound_min = min(m_bound_min, om2glm(m_mesh.point(*vit)));
			m_bound_max = max(m_bound_max, om2glm(m_mesh.point(*vit)));
		}

		std::cout << "Computing vertex areas" << std::endl;
		m_prop_vertex_area = computeVertexAreas(m_mesh);

		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_mesh);

		std::cout << "Computing raw don curvature" << std::endl;
		// TODO what if we didnt want this? (to run and time with a difference curv measure)
		const auto time_curv_start = std::chrono::steady_clock::now();
		computeDoNMaxDiffs(m_mesh, m_prop_doncurv_raw, m_prop_vertex_area, 1);
		const auto time_curv_finish = std::chrono::steady_clock::now();
		std::cout << "Curvature took " << ((time_curv_finish - time_curv_start) / std::chrono::duration<double>(1.0)) << "s" << std::endl;

		// TODO move autocontrast to a function
		std::cout << "Computing auto contrast for saliency" << std::endl;
		{
			// experimental auto contrast
			// TODO run this after decimation too
			static constexpr int hist_bits = 8;
			std::minstd_rand rand{std::random_device{}()};
			std::vector<int> samples(m_mesh.n_vertices());
			std::iota(samples.begin(), samples.end(), 0);
			std::shuffle(samples.begin(), samples.end(), rand);
			auto eval_entropy = [&](float contrast) {
				//std::cout << "computing entropy for contrast=" << contrast << std::endl;
				// use upper bound of 1 (and saliency now always uses 0-1 binning too)
				// (determining the upper bound exactly would need a prepass to calculate the range before binning)
				auto hist = maxdon_histogram<(1 << hist_bits)>(m_mesh, m_prop_vertex_area, contrast, 1.f, rand, samples);
				float s = histogram_entropy(hist);
				//std::cout << "entropy=" << s << std::endl;
				return s;
			};
			auto eval_entropy_gradient = [&](float contrast, float delta) {
				float s0 = eval_entropy(contrast);
				float s1 = eval_entropy(contrast + delta);
				return std::pair{s0, (s1 - s0) / delta};
			};
			const float target_entropy = float(hist_bits) * 0.45f;
			float best_contrast = 1;
			float contrast_delta = 0.1f;
			for (int i = 0; i < 5; i++) {
				auto [s, dsdc] = eval_entropy_gradient(best_contrast, contrast_delta);
				best_contrast = best_contrast - (s - target_entropy) / dsdc;
				contrast_delta *= 0.7f;
			}
			std::cout << "Auto contrast: " << best_contrast << std::endl;
			m_auto_contrast = best_contrast;
		}
	}

	Model::Model(const std::filesystem::path &fpath) : Model(ModelBase(fpath)) {
		
	}

	void Model::save(const std::filesystem::path &fpath, const model_save_params &sparams0) {
		auto sparams = sparams0;
		if (sparams.color_mode == model_color_mode::saliency) {
			if (!sparams.prop_saliency.is_valid()) {
				sparams.color_mode = model_color_mode::none;
			} else {
				std::cout << "Colorizing saliency" << std::endl;
			}
		}
		if (sparams.color_mode == model_color_mode::saliency_comparison) {
			if (!sparams.prop_saliency.is_valid()) {
				sparams.color_mode = model_color_mode::none;
			} else if (!sparams.prop_saliency_baseline.is_valid()) {
				sparams.color_mode = model_color_mode::none;
			} else {
				std::cout << "Colorizing saliency comparison" << std::endl;
			}
		}
		if (sparams.color_mode == model_color_mode::doncurv) {
			if (!m_prop_doncurv_raw.is_valid()) {
				sparams.color_mode == model_color_mode::none;
			} else {
				std::cout << "Colorizing don curvature" << std::endl;
			}
		}
		// assign vertex colors
		if (sparams.color_mode != model_color_mode::none) {
			m_mesh.request_vertex_colors();
			for (auto vIt = m_mesh.vertices_begin(); vIt != m_mesh.vertices_end(); ++vIt) {
				PolyMesh::Color col;
				OpenMesh::Vec3f v;
				bool usev = true;
				switch (sparams.color_mode) {
				case model_color_mode::saliency:
				{
					const float s = m_mesh.property(sparams.prop_saliency, *vIt);
					mapScalarToColor(v, s, TransferFunction::ZBRUSH);
					break;
				}
				case model_color_mode::saliency_comparison:
				{
					const float s = m_mesh.property(sparams.prop_saliency, *vIt);
					const float b = m_mesh.property(sparams.prop_saliency_baseline, *vIt);
					mapScalarToColor(v, std::clamp((s - b) * sparams.error_scale * 0.5f + 0.5f, 0.f, 1.f), TransferFunction::ZBRUSH);
					break;
				}
				case model_color_mode::vcolor:
				{
					col = m_mesh.property(m_prop_vcolor_original, *vIt);
					usev = false;
					break;
				}
				case model_color_mode::doncurv:
				{
					// note: need much lower contrast for display than for saliency
					const float c = std::pow(m_mesh.property(m_prop_doncurv_raw, *vIt), m_auto_contrast * 0.2f);
					mapScalarToColor(v, c, TransferFunction::ZBRUSH);
					break;
				}
				default:
					break;
				}
				if (usev) {
					v = v * 255;
					col[0] = v[0];
					col[1] = v[1];
					col[2] = v[2];
				}
				m_mesh.set_color(*vIt, col);
			}
		}
		// copy export saliency to named properties
		std::vector<saliency_prop_t> exprops;
		for (auto &sd : m_saliency) {
			if (!sd.should_export || !sd.prop_saliency.is_valid()) continue;
			saliency_prop_t p;
			const auto name = sd.export_propname(exprops.size());
			std::cout << "Export saliency property with internal name " << name << std::endl;
			m_mesh.add_property(p, name);
			exprops.push_back(p);
			// persistent => openmesh will export with Options::Custom
			m_mesh.property(p).set_persistent(true);
			// note - copying saliency data, not exactly efficient
			m_mesh.property(p).data_vector() = m_mesh.property(sd.prop_saliency).data_vector();
		}
		// copy original vids to named property
		OpenMesh::VPropHandleT<int> exprop_vid;
		if (sparams.original_vids && m_prop_vid_original.is_valid()) {
			std::cout << "Export property for original vertex ids" << std::endl;
			m_mesh.add_property(exprop_vid, "original_vid");
			m_mesh.property(exprop_vid).set_persistent(true);
			m_mesh.property(exprop_vid).data_vector() = m_mesh.property(m_prop_vid_original).data_vector();
		}
		// export!
		OpenMesh::IO::Options opts{};
		if (exprops.size()) opts = opts | OpenMesh::IO::Options::Custom;
		if (sparams.color_mode != model_color_mode::none) opts = opts | OpenMesh::IO::Options::VertexColor;
		if (sparams.binary) opts = opts | OpenMesh::IO::Options::Binary;
		std::cerr << "Saving model " << fpath.u8string() << std::endl;
		auto res = OpenMesh::IO::write_mesh(m_mesh, fpath, opts);
		// remove temp named properties
		m_mesh.remove_property(exprop_vid);
		for (auto &p : exprops) {
			m_mesh.remove_property(p);
		}
		// done
		if (res) {
			std::cerr << "Saved" << std::endl;
		} else {
			std::cerr << "Could not save model " << fpath.u8string() << std::endl;
			throw std::runtime_error("failed to save model");
		}
	}

	Model Model::prepare_decimate(saliency_prop_t prop_saliency, const std::vector<model_saliency_data> &sdv) const {
		Model m;
		m.m_bound_min = m_bound_min;
		m.m_bound_max = m_bound_max;
		// copy vertex positions and connectivity
		// FIXME if source has face status etc, those props are broken in the copy! (openmesh bug)
		if (m_mesh.has_face_status()) std::abort();
		m.m_mesh.assign(m_mesh);
		// copy original vertex colors if present
		if (m_prop_vcolor_original.is_valid()) {
			m.m_mesh.add_property(m.m_prop_vcolor_original);
			m.m_mesh.property(m.m_prop_vcolor_original).data_vector() = m_mesh.property(m_prop_vcolor_original).data_vector();
		}
		// TODO not currently copying existing original vids
		// (would need a nice way of choosing when to keep or regen)
		// copy specified saliency properties
		for (auto &sd : sdv) {
			model_saliency_data sd2 = sd;
			sd2.prop_saliency.invalidate();
			sd2.prop_sampled.invalidate();
			sd2.decimated = true;
			if (sd.prop_saliency.is_valid()) {
				m.m_mesh.add_property(sd2.prop_saliency);
				if (sd.prop_saliency == prop_saliency) m.m_prop_sal_dec = sd2.prop_saliency;
				m.m_mesh.property(sd2.prop_saliency).data_vector() = m_mesh.property(sd.prop_saliency).data_vector();
				m.m_saliency.push_back(std::move(sd2));
			}
		}
		return m;
	}

	bool Model::decimate(const decimate_user_params &uparams, decimate_progress &progress) {
		// generate original vertex ids if not present
		// never present atm
		if (!m_prop_vid_original.is_valid()) {
			m_mesh.add_property(m_prop_vid_original);
			for (auto v : m_mesh.vertices()) {
				m_mesh.property(m_prop_vid_original, v) = v.idx();
			}
		}
		// decimate!
		std::cout << "Decimating model" << std::endl;
		decimate_mesh_params mparams;
		mparams.mesh = &m_mesh;
		mparams.prop_saliency = m_prop_sal_dec;
		if (!green::decimate(mparams, uparams, progress)) return false;
		// recompute normals
		std::cout << "Computing vertex normals" << std::endl;
		m_mesh.request_face_normals();
		m_mesh.request_vertex_normals();
		m_mesh.update_face_normals();
		m_mesh.update_vertex_normals();
		// recompute vertex areas and edge length
		std::cout << "Computing vertex areas" << std::endl;
		m_prop_vertex_area = computeVertexAreas(m_mesh);
		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_mesh);
		return true;
	}

	void Model::update_vao() {
		bool makevao = false;
		if (!m_vao) makevao = true, m_vao = cgu::gl_object::gen_vertex_array();
		if (!m_ibo) m_ibo = cgu::gl_object::gen_buffer();
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		const size_t ntris = m_mesh.n_faces();
		const size_t nverts = m_mesh.n_vertices();
		assert(ntris <= size_t(INT_MAX));
		assert(nverts <= size_t(INT_MAX));
		m_vao_ntris = ntris;
		m_vao_nverts = nverts;
		const size_t size_ibo = ntris * 3 * sizeof(GLuint);
		glBindVertexArray(m_vao);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_ibo, nullptr, GL_STATIC_DRAW);
		if (ntris) {
			auto pibo = reinterpret_cast<GLuint *>(
				glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, size_ibo, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto fit = m_mesh.faces_begin(); fit != m_mesh.faces_end(); ++fit) {
				auto hit = m_mesh.cfv_iter(*fit);
				// assume triangles
				*pibo++ = GLuint(hit++->idx());
				*pibo++ = GLuint(hit++->idx());
				*pibo++ = GLuint(hit++->idx());
			}
			glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		}
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 4, GL_FLOAT, false, sizeof(glm::vec4), 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vbos() {
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBufferData(GL_ARRAY_BUFFER, m_mesh.n_vertices() * sizeof(glm::vec3), m_mesh.points(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
		glBufferData(GL_ARRAY_BUFFER, m_mesh.n_vertices() * sizeof(glm::vec3), m_mesh.vertex_normals(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		glBufferData(GL_ARRAY_BUFFER, m_mesh.n_vertices() * sizeof(glm::vec4), nullptr, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	bool Model::update_color(const model_color_params &cparams, model_saliency_errors *err) {
		model_saliency_errors err2{};
		if (!err) err = &err2;
		if (cparams.color_mode == model_color_mode::saliency && !cparams.prop_saliency.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::vcolor && !m_prop_vcolor_original.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::saliency_comparison && !cparams.prop_saliency.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::saliency_comparison && !cparams.prop_saliency_baseline.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::doncurv && !m_prop_doncurv_raw.is_valid()) return false;
		if (!m_vbo_col) return false;
		const auto nverts = vao_nverts();
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		auto *data = reinterpret_cast<glm::vec4 *>(
			glMapBufferRange(GL_ARRAY_BUFFER, 0, nverts * sizeof(glm::vec4), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
		);
		if (cparams.color_mode == model_color_mode::saliency) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float s = m_mesh.property(cparams.prop_saliency, v);
				data[i] = glm::vec4(s, 0, 0, 0);
			}
		} else if (cparams.color_mode == model_color_mode::saliency_comparison) {
			err->min = 9001;
			err->max = -9001;
			float sse = 0;
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float b = m_mesh.property(cparams.prop_saliency_baseline, v);
				const float s = m_mesh.property(cparams.prop_saliency, v);
				const float e = s - b;
				err->min = std::min(err->min, e);
				err->max = std::max(err->max, e);
				sse += e * e;
				// TODO do scale on gpu
				data[i] = glm::vec4(e * cparams.error_scale, 0, 0, 0);
			}
			err->rms = sqrt(sse / nverts);
		} else if (cparams.color_mode == model_color_mode::vcolor) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const auto col = m_mesh.property(m_prop_vcolor_original, v);
				// need to gamma decode the color because the shader expects linear
				data[i] = glm::vec4(pow(glm::vec3(col[0], col[1], col[2]) / 255.f, glm::vec3(2.2f)), 1);
			}
		} else if (cparams.color_mode == model_color_mode::doncurv) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				// note: need much lower contrast for display than for saliency
				const float c = std::pow(m_mesh.property(m_prop_doncurv_raw, v), auto_contrast() * 0.2f);
				data[i] = glm::vec4(c, 0, 0, 0);
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		return true;
	}

	void Model::draw(GLenum polymode) const {
		if (!m_vao) return;
		glBindVertexArray(m_vao);
		if (polymode == GL_POINT) {
			glDrawArrays(GL_POINTS, 0, m_vao_nverts);
		} else {
			glPolygonMode(GL_FRONT_AND_BACK, polymode);
			glDrawElements(GL_TRIANGLES, m_vao_ntris * 3, GL_UNSIGNED_INT, 0);
		}
		glBindVertexArray(0);
	}

	void Model::draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params, GLenum polymode) const {

		static GLuint prog_smooth = 0;
		static GLuint prog_flat = 0;
		if (!params.shade_flat && !prog_smooth) {
			// note: needs custom depth env
			prog_smooth = cgu::make_shader_program(
				"green::model.smooth",
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{"#define SHADE_SMOOTH\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}
		if (params.shade_flat && !prog_flat) {
			// note: needs custom depth env
			prog_flat = cgu::make_shader_program(
				"green::model.flat",
				"330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{"#define SHADE_FLAT\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}

		GLuint prog = params.shade_flat ? prog_flat : prog_smooth;
		glUseProgram(prog);
		glUniform1f(glGetUniformLocation(prog, cgu::glsl_uniform_zfar_name), zfar);
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(modelview));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(projection));
		glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, value_ptr(params.color));
		glm::vec3 pos_bias = -bound_center();
		glUniform3fv(glGetUniformLocation(prog, "u_pos_bias"), 1, value_ptr(pos_bias));
		glUniform1f(glGetUniformLocation(prog, "u_shading"), params.shading);
		float bias = 0;
		// biases adjusted for log-depth
		if (polymode == GL_LINE) bias = -0.00001;
		if (polymode == GL_POINT) bias = -0.00002;
		glUniform1f(glGetUniformLocation(prog, "u_depth_bias"), bias);
		glUniform1i(glGetUniformLocation(prog, "u_entity_id"), params.entity_id);
		glUniform1i(glGetUniformLocation(prog, "u_vert_color_map"), params.vert_color_map);
		glUniform1i(glGetUniformLocation(prog, "u_show_samples"), params.show_samples);

		draw(polymode);
	}

}
