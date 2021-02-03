
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

	auto & mesh_io_mutex() {
		static std::mutex m;
		return m;
	}

	struct saliency_clipboard_data {
		green::model_saliency_data sd;
		std::vector<float> data;
	};

	auto & saliency_clipboard() {
		static saliency_clipboard_data d;
		return d;
	}

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
		const green::TriMesh &mesh,
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
			TriMesh::VertexHandle v(shuffled_samples[i]);
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
		m_trimesh.request_face_normals();
		m_trimesh.request_vertex_normals();
		m_trimesh.request_vertex_colors();
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

		if (OpenMesh::IO::read_mesh(m_trimesh, fpath, readOptions)) {
			std::cerr << "Loaded" << std::endl;
		} else {
			std::cerr << "Failed" << std::endl;
			throw std::runtime_error("failed to load model");
		}

		// TODO necessary? optional?
		// NOTE currently done by assimp too
		m_trimesh.triangulate();
		
		// copy original vertex colors
		// (because we need to be able to overwrite the actual vertex colors during export)
		if (readOptions.check(OpenMesh::IO::Options::VertexColor)) {
			std::cout << "Found vertex colors" << std::endl;
			m_trimesh.add_property(m_prop_vcolor_original);
			for (auto vIt = m_trimesh.vertices_begin(); vIt != m_trimesh.vertices_end(); ++vIt) {
				m_trimesh.property(m_prop_vcolor_original, *vIt) = m_trimesh.color(*vIt);
			}
		}

		// check for 'raw' saliency properties
		std::vector<std::pair<saliency_prop_t, model_saliency_data>> sprops;
		for (auto it = m_trimesh.vprops_begin(); it != m_trimesh.vprops_end(); ++it) {
			// TODO any way to do this other than dynamic cast?
			auto prop = dynamic_cast<OpenMesh::PropertyT<float> *>(*it);
			if (!prop) continue;
			if (prop->name().substr(0, 7) == "quality") {
				std::cout << "Found saliency property with internal name " << prop->name() << std::endl;
				auto ph0 = saliency_prop_t{int(it - m_trimesh.vprops_begin())};
				if (&m_trimesh.property(ph0) != prop) abort();
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
			m_trimesh.add_property(ph);
			m_trimesh.property(ph).data_vector() = std::move(m_trimesh.property(p.first).data_vector());
			m_trimesh.remove_property(p.first);
			auto sd = std::move(p.second);
			sd.prop_saliency = ph;
			m_saliency.push_back(std::move(sd));
		}

		// move original vids to unnamed property
		if (OpenMesh::VPropHandleT<int> pvid; m_trimesh.get_property_handle(pvid, "original_vid")) {
			std::cout << "Found original vertex ids property" << std::endl;
			m_trimesh.add_property(m_prop_vid_original);
			m_trimesh.property(m_prop_vid_original).data_vector() = std::move(m_trimesh.property(pvid).data_vector());
			m_trimesh.remove_property(pvid);
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
		m_trimesh.update_face_normals();
		m_trimesh.update_vertex_normals();

		std::cout << "Computing bounding box" << std::endl;
		for (auto vit = m_trimesh.vertices_begin(); vit != m_trimesh.vertices_end(); ++vit) {
			m_bound_min = min(m_bound_min, om2glm(m_trimesh.point(*vit)));
			m_bound_max = max(m_bound_max, om2glm(m_trimesh.point(*vit)));
		}

		std::cout << "Computing vertex areas" << std::endl;
		m_prop_vertex_area = computeVertexAreas(m_trimesh);

		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_trimesh);

		std::cout << "Computing raw don curvature" << std::endl;
		// TODO what if we didnt want this? (to run and time with a difference curv measure)
		const auto time_curv_start = std::chrono::steady_clock::now();
		// TODO do away with this histogram thing
		lce::Histogram curvhist;
		computeDoNMaxDiffs(m_trimesh, m_prop_doncurv_raw, curvhist, m_prop_vertex_area, 1);
		const auto time_curv_finish = std::chrono::steady_clock::now();
		std::cout << "Curvature took " << ((time_curv_finish - time_curv_start) / std::chrono::duration<double>(1.0)) << "s" << std::endl;
		m_don_raw_min = curvhist.getMin();
		m_don_raw_max = curvhist.getMax();

		// TODO move autocontrast to a function
		std::cout << "Computing auto contrast for saliency" << std::endl;
		{
			// experimental auto contrast
			// TODO run this after decimation too
			static constexpr int hist_bits = 8;
			std::minstd_rand rand{std::random_device{}()};
			std::vector<int> samples(m_trimesh.n_vertices());
			std::iota(samples.begin(), samples.end(), 0);
			std::shuffle(samples.begin(), samples.end(), rand);
			auto eval_entropy = [&](float contrast) {
				//std::cout << "computing entropy for contrast=" << contrast << std::endl;
				// use upper bound of 1 (and saliency now always uses 0-1 binning too)
				// (determining the upper bound exactly would need a prepass to calculate the range before binning)
				auto hist = maxdon_histogram<(1 << hist_bits)>(m_trimesh, m_prop_vertex_area, contrast, 1.f, rand, samples);
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
			m_trimesh.request_vertex_colors();
			for (auto vIt = m_trimesh.vertices_begin(); vIt != m_trimesh.vertices_end(); ++vIt) {
				TriMesh::Color col;
				OpenMesh::Vec3f v;
				bool usev = true;
				switch (sparams.color_mode) {
				case model_color_mode::saliency:
				{
					const float s = m_trimesh.property(sparams.prop_saliency, *vIt);
					mapScalarToColor(v, s, TransferFunction::ZBRUSH);
					break;
				}
				case model_color_mode::saliency_comparison:
				{
					const float s = m_trimesh.property(sparams.prop_saliency, *vIt);
					const float b = m_trimesh.property(sparams.prop_saliency_baseline, *vIt);
					mapScalarToColor(v, std::clamp((s - b) * sparams.error_scale * 0.5f + 0.5f, 0.f, 1.f), TransferFunction::ZBRUSH);
					break;
				}
				case model_color_mode::vcolor:
				{
					col = m_trimesh.property(m_prop_vcolor_original, *vIt);
					usev = false;
					break;
				}
				case model_color_mode::doncurv:
				{
					// note: need much lower contrast for display than for saliency
					const float c = std::pow(m_trimesh.property(m_prop_doncurv_raw, *vIt), m_auto_contrast * 0.2f);
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
				m_trimesh.set_color(*vIt, col);
			}
		}
		// copy export saliency to named properties
		std::vector<saliency_prop_t> exprops;
		for (auto &sd : m_saliency) {
			if (!sd.should_export || !sd.prop_saliency.is_valid()) continue;
			saliency_prop_t p;
			const auto name = sd.export_propname(exprops.size());
			std::cout << "Export saliency property with internal name " << name << std::endl;
			m_trimesh.add_property(p, name);
			exprops.push_back(p);
			// persistent => openmesh will export with Options::Custom
			m_trimesh.property(p).set_persistent(true);
			// note - copying saliency data, not exactly efficient
			m_trimesh.property(p).data_vector() = m_trimesh.property(sd.prop_saliency).data_vector();
		}
		// copy original vids to named property
		OpenMesh::VPropHandleT<int> exprop_vid;
		if (sparams.original_vids && m_prop_vid_original.is_valid()) {
			std::cout << "Export property for original vertex ids" << std::endl;
			m_trimesh.add_property(exprop_vid, "original_vid");
			m_trimesh.property(exprop_vid).set_persistent(true);
			m_trimesh.property(exprop_vid).data_vector() = m_trimesh.property(m_prop_vid_original).data_vector();
		}
		// export!
		OpenMesh::IO::Options opts{};
		if (exprops.size()) opts = opts | OpenMesh::IO::Options::Custom;
		if (sparams.color_mode != model_color_mode::none) opts = opts | OpenMesh::IO::Options::VertexColor;
		if (sparams.binary) opts = opts | OpenMesh::IO::Options::Binary;
		std::cerr << "Saving model " << fpath.u8string() << std::endl;
		auto res = OpenMesh::IO::write_mesh(m_trimesh, fpath, opts);
		// remove temp named properties
		m_trimesh.remove_property(exprop_vid);
		for (auto &p : exprops) {
			m_trimesh.remove_property(p);
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
		if (m_trimesh.has_face_status()) std::abort();
		m.m_trimesh.assign(m_trimesh);
		// copy original vertex colors if present
		if (m_prop_vcolor_original.is_valid()) {
			m.m_trimesh.add_property(m.m_prop_vcolor_original);
			m.m_trimesh.property(m.m_prop_vcolor_original).data_vector() = m_trimesh.property(m_prop_vcolor_original).data_vector();
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
				m.m_trimesh.add_property(sd2.prop_saliency);
				if (sd.prop_saliency == prop_saliency) m.m_prop_sal_dec = sd2.prop_saliency;
				m.m_trimesh.property(sd2.prop_saliency).data_vector() = m_trimesh.property(sd.prop_saliency).data_vector();
				m.m_saliency.push_back(std::move(sd2));
			}
		}
		return m;
	}

	bool Model::decimate(const decimate_user_params &uparams, decimate_progress &progress) {
		// generate original vertex ids if not present
		// never present atm
		if (!m_prop_vid_original.is_valid()) {
			m_trimesh.add_property(m_prop_vid_original);
			for (auto v : m_trimesh.vertices()) {
				m_trimesh.property(m_prop_vid_original, v) = v.idx();
			}
		}
		// decimate!
		std::cout << "Decimating model" << std::endl;
		decimate_mesh_params mparams;
		mparams.mesh = &m_trimesh;
		mparams.prop_saliency = m_prop_sal_dec;
		if (!green::decimate(mparams, uparams, progress)) return false;
		// recompute normals
		std::cout << "Computing vertex normals" << std::endl;
		m_trimesh.request_face_normals();
		m_trimesh.request_vertex_normals();
		m_trimesh.update_face_normals();
		m_trimesh.update_vertex_normals();
		// recompute vertex areas and edge length
		std::cout << "Computing vertex areas" << std::endl;
		m_prop_vertex_area = computeVertexAreas(m_trimesh);
		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_trimesh);
		return true;
	}

	void Model::update_vao() {
		bool makevao = false;
		if (!m_vao) makevao = true, m_vao = cgu::gl_object::gen_vertex_array();
		if (!m_ibo) m_ibo = cgu::gl_object::gen_buffer();
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		const size_t ntris = m_trimesh.n_faces();
		const size_t nverts = m_trimesh.n_vertices();
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
			for (auto fit = m_trimesh.faces_begin(); fit != m_trimesh.faces_end(); ++fit) {
				auto hit = m_trimesh.cfv_iter(*fit);
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
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec3), m_trimesh.points(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec3), m_trimesh.vertex_normals(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec4), nullptr, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
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

	ModelEntity::ModelEntity() {
		
	}

	void ModelEntity::load(const std::filesystem::path &fpath) {
		// lock in case of premature closure
		// mainly because the pending_load future's dtor will block otherwise
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async([=, lock=std::move(lock)]() {
			// apparently openmesh load is not threadsafe?
			// TODO check assimp too
			std::lock_guard iolock(mesh_io_mutex());
			return std::make_unique<Model>(fpath);
		});
	}

	void ModelEntity::save(const std::filesystem::path &fpath) {
		if (!m_model) return;
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		const auto &saliency = m_model->saliency();
		const auto sdcolor = m_saliency_export_index < saliency.size() ? saliency[m_saliency_export_index] : model_saliency_data{};
		const auto sdbase = m_saliency_baseline_index < saliency.size() ? saliency[m_saliency_baseline_index] : model_saliency_data{};
		m_fpath_save = fpath;
		m_pending_save = std::async([=, lock=std::move(lock)]() {
			// TODO is this also not threadsafe? idk
			// ... can it run parallel with load?
			std::lock_guard iolock(mesh_io_mutex());
			model_save_params sparams;
			sparams.prop_saliency = sdcolor.prop_saliency;
			sparams.prop_saliency_baseline = sdbase.prop_saliency;
			sparams.color_mode = m_exp_color_mode;
			sparams.original_vids = m_save_original_vids;
			sparams.binary = m_save_binary;
			m_model->save(fpath, sparams);
			return true;
		});
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m) {
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async(std::launch::deferred, [m=std::move(m)]() mutable {
			return std::make_unique<Model>(std::move(m));
		});
		// non-blocking wait will never succeed otherwise
		m_pending_load.wait();
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress) {
		load(fpath, std::move(m));
		m_decimated = true;
		m_dec_uparams = dec_uparams;
		m_dec_progress = dec_progress;
	}

	void ModelEntity::move_by(const glm::vec3 &d) {
		m_translation += d;
		if (d != glm::vec3(0)) invalidate_scene();
	}

	glm::mat4 ModelEntity::transform() const {
		if (!m_model) return glm::mat4(1);
		glm::mat4 transform(1);
		transform = glm::translate(transform, m_translation);
		transform *= glm::eulerAngleYXZ(m_rotation_euler_yxz.y, m_rotation_euler_yxz.x, m_rotation_euler_yxz.z);
		transform = glm::scale(transform, glm::vec3(m_scale));
		glm::mat3 basis(m_basis_vectors[m_basis_right].v, m_basis_vectors[m_basis_up].v, m_basis_vectors[m_basis_back].v);
		return transform * glm::mat4(transpose(basis));
	}

	void ModelEntity::update_vbo() {
		if (!m_model) return;
		if (!m_saliency_vbo_dirty) return;
		std::shared_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) return;
		auto &saliency_outputs = m_model->saliency();
		if (m_disp_color_mode == model_color_mode::saliency && m_saliency_index >= saliency_outputs.size()) return;
		if (m_disp_color_mode == model_color_mode::vcolor && !m_model->prop_vcolor_original().is_valid()) return;
		if (m_disp_color_mode == model_color_mode::saliency_comparison && m_saliency_index >= saliency_outputs.size()) return;
		if (m_disp_color_mode == model_color_mode::saliency_comparison && m_saliency_baseline_index >= saliency_outputs.size()) return;
		if (m_disp_color_mode == model_color_mode::doncurv && !m_model->prop_doncurv_raw().is_valid()) return;
		GLuint vbo_col = m_model->vbo_color();
		if (!vbo_col) return;
		const auto &mesh = m_model->trimesh();
		const auto nverts = m_model->vao_nverts();
		glBindBuffer(GL_ARRAY_BUFFER, vbo_col);
		auto *data = reinterpret_cast<glm::vec4 *>(
			glMapBufferRange(GL_ARRAY_BUFFER, 0, nverts * sizeof(glm::vec4), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
		);
		if (m_disp_color_mode == model_color_mode::saliency) {
			auto salprop = saliency_outputs[m_saliency_index].prop_saliency;
			auto sampledprop = saliency_outputs[m_saliency_index].prop_sampled;
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float s = mesh.property(salprop, v);
				const bool sampled = sampledprop.is_valid() && mesh.property(sampledprop, v);
				data[i] = glm::vec4(s, 0, 0, sampled);
			}
		} else if (m_disp_color_mode == model_color_mode::saliency_comparison) {
			auto baseprop = saliency_outputs[m_saliency_baseline_index].prop_saliency;
			auto salprop = saliency_outputs[m_saliency_index].prop_saliency;
			auto sampledprop = saliency_outputs[m_saliency_index].prop_sampled;
			m_saliency_errors.min = 9001;
			m_saliency_errors.max = -9001;
			float sse = 0;
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float b = mesh.property(baseprop, v);
				const float s = mesh.property(salprop, v);
				const bool sampled = sampledprop.is_valid() && mesh.property(sampledprop, v);
				const float e = s - b;
				m_saliency_errors.min = std::min(m_saliency_errors.min, e);
				m_saliency_errors.max = std::max(m_saliency_errors.max, e);
				sse += e * e;
				// TODO do scale on gpu
				data[i] = glm::vec4(e * m_saliency_error_scale, 0, 0, sampled);
			}
			m_saliency_errors.rms = sqrt(sse / nverts);
		} else if (m_disp_color_mode == model_color_mode::vcolor) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const auto col = mesh.property(m_model->prop_vcolor_original(), v);
				// need to gamma decode the color because the shader expects linear
				data[i] = glm::vec4(pow(glm::vec3(col[0], col[1], col[2]) / 255.f, glm::vec3(2.2f)), 1);
			}
		} else if (m_disp_color_mode == model_color_mode::doncurv) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				// note: need much lower contrast for display than for saliency
				const float c = std::pow(mesh.property(m_model->prop_doncurv_raw(), v), m_model->auto_contrast() * 0.2f);
				data[i] = glm::vec4(c, 0, 0, 0);
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		m_saliency_vbo_dirty = false;
		invalidate_scene();
	}

	void ModelEntity::draw_window_models(bool selected) {
		using namespace ImGui;
		if (Begin("Models")) {
			PushID(this);
			bool dirty = false;
			draw_select_header(selected);
			if (m_pending_load.valid()) {
				Text("Loading...");
			} else if (m_model) {
				// ok
			} else {
				Text("Failed to load model");
			}
			dirty |= Checkbox("Faces", &m_show_faces);
			SameLine();
			dirty |= Checkbox("Edges", &m_show_edges);
			SameLine();
			dirty |= Checkbox("Verts", &m_show_verts);
			SameLine();
			SetCursorPosX(GetCursorPosX() + std::max(0.f, GetContentRegionAvail().x - 20));
			PushStyleColor(ImGuiCol_Button, ImVec4{0.6f, 0.3f, 0.3f, 1});
			if (Button("X", {-1, 0}) || m_try_kill) {
				OpenPopup("##close");
				m_try_kill = false;
			}
			if (IsItemHovered()) SetTooltip("Close model");
			PopStyleColor();
			if (BeginPopupModal("##close", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoDecoration)) {
				Text("Close model \"%s\" ?", name().c_str());
				SetHoveredTooltip("%s", make_name_tooltip().c_str());
				if (Button("Close")) {
					std::unique_lock lock(m_modelmtx, std::defer_lock);
					if (lock.try_lock()) {
						m_dead = true;
						invalidate_scene();
					} else {
						spawn_locked_notification();
					}
					CloseCurrentPopup();
				}
				SameLine();
				if (Button("Cancel")) {
					CloseCurrentPopup();
				}
				EndPopup();
			}
			Separator();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
			if (dirty) invalidate_scene();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_selection() {
		using namespace ImGui;
		if (Begin("Selection")) {
			PushID(this);
			PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
			Selectable(name().c_str(), true, 0, {0, GetTextLineHeightWithSpacing()});
			PopStyleColor();
			SetHoveredTooltip("%s", make_name_tooltip().c_str());
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);

			Text("%zd vertices, %zd triangles", m_model->trimesh().n_vertices(), m_model->trimesh().n_faces());
			if (m_decimated) Text(
				"Decimated %d vertices (%.1f%%) in %.3fs",
				m_dec_progress.completed_collapses,
				100.f * m_dec_progress.completed_collapses / float(m_dec_uparams.targetverts + m_dec_progress.completed_collapses),
				m_dec_progress.elapsed_time / std::chrono::duration<double>(1.0)
			);

			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			if (m_pending_save.valid()) {
				if (m_pending_save.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
					try {
						m_save_ok = m_pending_save.get();
					} catch (std::exception &e) {
						std::cerr << "failed to save model: " << e.what() << std::endl;
						m_save_ok = false;
					} catch (...) {
						std::cerr << "failed to save model" << std::endl;
						m_save_ok = false;
					}
				} else {
					TextDisabled("Exporting %s...", m_fpath_save.filename().u8string().c_str());
				}
			} else if (m_save_ok) {
				TextDisabled("Export %s succeeded", m_fpath_save.filename().u8string().c_str());
			} else if (!m_fpath_save.empty()) {
				TextColored(badcol, "Export %s failed", m_fpath_save.filename().u8string().c_str());
			} else {
				TextDisabled("Not exported");
			}
			if (IsItemHovered() && !m_fpath_save.empty()) SetTooltip("%s", m_fpath_save.u8string().c_str());

			if (CollapsingHeader("Transform")) {
				auto pick_basis = [&](const char *label, int *basis) {
					SetNextItemWidth(GetTextLineHeight() * 3.5f);
					Combo(
						label, basis,
						[](void *data, int item, const char **out_text) {
							auto vs = reinterpret_cast<basis_vector *>(data);
							*out_text = vs[item].name;
							return true;
						}, m_basis_vectors, 6
					);
				};
				pick_basis("Right", &m_basis_right);
				SameLine();
				pick_basis("Up", &m_basis_up);
				SameLine();
				pick_basis("Back", &m_basis_back);
				if (Button("Reset##scale")) m_scale = m_model->unit_bound_scale() * 4;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderFloat("Scale", &m_scale, 0, 1000, "%.4f", 8);
				if (Button("Reset##translation")) m_translation = glm::vec3{0};
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				InputFloat3("Translation", value_ptr(m_translation));
				if (Button("Reset##roty")) m_rotation_euler_yxz.y = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Yaw", &m_rotation_euler_yxz.y, -180, 180);
				if (Button("Reset##rotx")) m_rotation_euler_yxz.x = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Pitch", &m_rotation_euler_yxz.x, -180, 180);
				if (Button("Reset##rotz")) m_rotation_euler_yxz.z = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Roll", &m_rotation_euler_yxz.z, -180, 180);
			}

			if (CollapsingHeader("Rendering")) {
				bool r = false;
				if (RadioButton("Smooth Shading", !m_shade_flat)) {
					m_shade_flat = false;
					r = true;
				}
				SameLine();
				if (RadioButton("Flat Shading", m_shade_flat)) {
					m_shade_flat = true;
					r = true;
				}
				r |= Checkbox("Color Faces", &m_color_faces);
				SetHoveredTooltip("Apply color mode to faces");
				SameLine();
				r |= Checkbox("Color Verts", &m_color_verts);
				SetHoveredTooltip("Apply color mode to vertices");
				r |= Checkbox("Cull Faces", &m_cull_faces);
				SameLine();
				r |= Checkbox("Cull Edges", &m_cull_edges);
				r |= SliderInt("Point Size", &m_vert_point_size, 1, 5);
				m_vert_point_size = std::max(1, m_vert_point_size);
				if (r) invalidate_scene();
			}

			Separator();

			//SetNextItemWidth(GetContentRegionAvail().x);
			int cur_color_mode = int(m_disp_color_mode);
			if (Combo("Color Mode", &cur_color_mode, "None\0Vertex Color\0Saliency\0Saliency Comparison\0Curvature (DoN)\0")) {
				m_disp_color_mode = model_color_mode(cur_color_mode);
				invalidate_saliency_vbo();
			}

			Separator();

			std::vector<model_saliency_data> empty_saliency_outputs;
			auto &saliency_outputs = m_model ? m_model->saliency() : empty_saliency_outputs;

			if (Combo(
				"Saliency", &m_saliency_index,
				[](void *data, int item, const char **out_text) {
					auto ds = reinterpret_cast<model_saliency_data *>(data);
					const auto &so = ds[item];
					// TODO better?
					static std::string str;
					str = std::string(so);
					*out_text = str.c_str();
					return true;
				}, 
				saliency_outputs.data(), saliency_outputs.size()
			)) {
				invalidate_saliency_vbo();
			}

			if (m_disp_color_mode == model_color_mode::saliency_comparison && m_saliency_baseline_index < saliency_outputs.size()) {
				Text("Baseline: %s", saliency_outputs[m_saliency_baseline_index].str().c_str());
				if (m_saliency_index < saliency_outputs.size()) {
					auto &err = m_saliency_errors;
					Text("Errors: min=%.3f, max=%.3f, rmse=%.3f", err.min, err.max, err.rms);
					if (SliderFloat("Error Scale", &m_saliency_error_scale, 1, 100, "%.3f", 2)) {
						invalidate_saliency_vbo();
					}
				}
			}

			if (Button("Paste")) {
				auto &clip = saliency_clipboard();
				std::unique_lock lock(m_modelmtx, std::defer_lock);
				if (!lock.try_lock()) {
					spawn_locked_notification();
				} else if (clip.data.size() == m_model->trimesh().n_vertices()) {
					model_saliency_data sd = std::move(clip.sd);
					m_model->trimesh().add_property(sd.prop_saliency);
					m_model->trimesh().property(sd.prop_saliency).data_vector() = std::move(clip.data);
					clip.sd = {};
					clip.data.clear();
					saliency_outputs.push_back(std::move(sd));
					invalidate_saliency_vbo();
					// references/iterators into saliency outputs are invalidated
				} else {
					OpenPopup("Paste Error##pasteerror");
				}
			}

			if (m_saliency_index < saliency_outputs.size()) {
				auto &salout = saliency_outputs[m_saliency_index];
				SameLine();
				if (Button("Copy")) {
					std::shared_lock lock(m_modelmtx, std::defer_lock);
					if (!lock.try_lock()) {
						spawn_locked_notification();
					} else {
						auto &clip = saliency_clipboard();
						clip.sd = salout;
						clip.sd.prop_saliency.reset();
						clip.data = m_model->trimesh().property(salout.prop_saliency).data_vector();
					}
				}
				SameLine();
				if (Button("Remove")) OpenPopup("##remove");
				SameLine();
				if (Button("Rename")) OpenPopup("Rename Saliency##rename");
				SameLine();
				if (Button("Baseline")) m_saliency_baseline_index = m_saliency_index;
				Checkbox("Persistent", &salout.persistent);
				SetHoveredTooltip("Persistent properties will be preserved when decimating\nand will be exported by default");
				Separator();
				if (CollapsingHeader("Saliency Parameters")) {
					if (salout.uparams_known) {
						draw_saliency_params(salout.uparams);
						if (Button("Reload")) ui_saliency_user_params() = salout.uparams;
					} else {
						TextDisabled("Parameters unknown");
					}
				}
				if (salout.filename.empty()) {
					if (CollapsingHeader("Saliency Progress")) {
						draw_saliency_progress(salout.progress);
					}
				}
				if (BeginPopupModal("Paste Error##pasteerror")) {
					auto &clip = saliency_clipboard();
					if (clip.data.size() == m_model->trimesh().n_vertices()) {
						CloseCurrentPopup();
					} else {
						Text("Can't paste saliency for %d vertices into model with %d vertices", clip.data.size(), m_model->trimesh().n_vertices());
						if (Button("OK")) CloseCurrentPopup();
					}
				}
				bool rename_open = true;
				SetNextWindowSize({300, 100}, ImGuiCond_Appearing);
				if (BeginPopupModal("Rename Saliency##rename", &rename_open)) {
					Text("New saliency property name");
					// TODO better?
					static char buf[1024]{};
					if (IsWindowAppearing()) {
						strncpy(buf, salout.dispname.c_str(), sizeof(buf));
						SetKeyboardFocusHere();
					}
					SetNextItemWidth(GetContentRegionAvail().x);
					const char *badchars = " \t@[]";
					if (InputText("", buf, sizeof(buf), ImGuiInputTextFlags_EnterReturnsTrue)) {
						if (std::string_view(buf).find_first_of(badchars) == std::string::npos) {
							salout.dispname = buf;
							CloseCurrentPopup();
						}
					}
					if (std::string_view(buf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of: %s", badchars);
					}
					EndPopup();
				}
				if (BeginPopupModal("##remove")) {
					Text("Remove saliency result?");
					if (Button("Remove")) {
						std::unique_lock lock(m_modelmtx, std::defer_lock);
						if (!lock.try_lock()) {
							spawn_locked_notification();
						} else {
							auto it = saliency_outputs.begin() + m_saliency_index;
							m_model->trimesh().remove_property(it->prop_saliency);
							saliency_outputs.erase(it);
							// references/iterators into saliency outputs are invalidated
							// (so this section should come last)
							if (m_saliency_index >= saliency_outputs.size()) {
								m_saliency_index = std::max(0, m_saliency_index - 1);
								invalidate_saliency_vbo();
							}
						}
						CloseCurrentPopup();
					}
					SameLine();
					if (Button("Cancel")) CloseCurrentPopup();
					EndPopup();
				}
			} else {
				SameLine();
				ButtonDisabled("Copy");
				SameLine();
				ButtonDisabled("Remove");
				SameLine();
				ButtonDisabled("Rename");
				SameLine();
				ButtonDisabled("Baseline");
			}
			Separator();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_export() {
		using namespace ImGui;
		if (!m_model) return;
		auto &saliency_outputs = m_model->saliency();
		if (m_try_export) {
			m_try_export = false;
			m_saliency_export_index = m_saliency_index;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				sd.should_export = sd.persistent;
			}
			OpenPopup("Export##export");
		}
		SetNextWindowSize({500, 400}, ImGuiCond_Appearing);
		bool export_window_open = true;
		if (BeginPopupModal("Export##export", &export_window_open)) {
			PushID(this);
			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
			Selectable(m_fpath_load.filename().u8string().c_str(), true, ImGuiSelectableFlags_DontClosePopups, {0, GetTextLineHeightWithSpacing()});
			PopStyleColor();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
			SetHoveredTooltip("%s", make_name_tooltip().c_str());
			// export path
			const auto &pathhint = m_fpath_save.empty() ? m_fpath_load : m_fpath_save;
			// TODO better?
			static char pathbuf[1024]{};
			if (IsWindowAppearing()) snprintf(pathbuf, sizeof(pathbuf), "%s", pathhint.u8string().c_str());
			auto &fpath_save_fut = ui_save_path(pathhint, Button("Browse"));
			if (fpath_save_fut.valid() && fpath_save_fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				auto s = fpath_save_fut.get().u8string();
				if (!s.empty()) snprintf(pathbuf, sizeof(pathbuf), "%s", s.c_str());
			}
			SameLine();
			InputText("Path", pathbuf, sizeof(pathbuf));
			// export color mode
			int cur_color_mode = int(m_exp_color_mode);
			if (Combo("Color Mode", &cur_color_mode, "None\0Vertex Color\0Saliency\0Saliency Comparison\0Curvature (DoN)\0")) {
				m_exp_color_mode = model_color_mode(cur_color_mode);
			}
			if (m_exp_color_mode == model_color_mode::saliency_comparison) {
				if (m_saliency_baseline_index < m_model->saliency().size()) {
					Text("Baseline: %s", m_model->saliency()[m_saliency_baseline_index].str().c_str());
				} else {
					TextColored(badcol, "Invalid baseline");
				}
			}
			Separator();
			Text("Saliency Properties");
			if (BeginChild("saliency", {0, -100}, false, ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
				for (int i = 0; i < saliency_outputs.size(); i++) {
					PushID(i);
					auto &sd = saliency_outputs[i];
					Checkbox("##shouldexport", &sd.should_export);
					SetHoveredTooltip("Export this saliency property?");
					SameLine();
					if (RadioButton("##colorize", i == m_saliency_export_index)) m_saliency_export_index = i;
					SetHoveredTooltip("Colorize this saliency property?");
					const char *badchars = " \t@[]";
					char propbuf[256]{};
					strncpy(propbuf, sd.expname.c_str(), sizeof(propbuf));
					SameLine();
					SetNextItemWidth(150);
					if (InputText("##expname", propbuf, sizeof(propbuf))) {
						if (std::string_view(propbuf).find_first_of(badchars) == std::string::npos) {
							sd.expname = propbuf;
						}
					}
					SetHoveredTooltip("Export property name.\nDefault when empty is the current display name.");
					SameLine();
					if (sd.expname.empty()) {
						auto p0 = GetCursorPos();
						SetCursorPosX(p0.x - 150);
						TextDisabled("%s", sd.dispname.c_str());
						SameLine();
						SetCursorPosX(p0.x);
					}
					if (Button("<")) sd.expname = sd.dispname;
					SetHoveredTooltip("Copy display name to export name");
					SameLine();
					Text("%s", sd.str().c_str());
					SetHoveredTooltip(sd.uparams_known ? sd.uparams.str().c_str() : "Parameters unknown");
					if (std::string_view(propbuf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of \"%s\"", badchars);
					}
					PopID();
				}
			}
			EndChild();
			Separator();
			// other options
			Checkbox("Binary", &m_save_binary);
			SetHoveredTooltip("Save file as binary if supported");
			SameLine();
			Checkbox("Original Vertex IDs", &m_save_original_vids);
			SetHoveredTooltip("Export extra property with original vertex ids from before decimation");
			TextDisabled("Supported formats are currently PLY and OBJ.\nColorization and property export are only supported for PLY.");
			// validate path and maybe export
			auto fpath = std::filesystem::u8path(pathbuf);
			auto stat = std::filesystem::status(fpath);
			bool cansave = !fpath.empty();
			const char *badchars = "\\/:*?\"<>|";
			if (fpath.is_relative()) {
				TextColored(badcol, "Path must be absolute");
				cansave = false;
			} else if (std::filesystem::is_directory(stat)) {
				TextColored(badcol, "Path is a directory");
				cansave = false;
			} else if (!std::filesystem::is_directory(fpath.parent_path())) {
				TextColored(badcol, "Directory does not exist");
				cansave = false;
			} else if (fpath.filename().u8string().find_first_of(badchars) != std::string::npos) {
				TextColored(badcol, "File name must not contain any of \"%s\"", badchars);
				cansave = false;
			} else if (std::filesystem::exists(stat)) {
				TextColored(badcol, "Path exists! Save will overwrite");
				//TextColored(badcol, u8"保存先が存在します！保存したら上書きします");
			}
			if (Button("Save") && cansave && !m_pending_save.valid()) {
				save(std::filesystem::absolute(fpath));
				CloseCurrentPopup();
			}
			if (IsItemHovered()) {
				if (m_pending_save.valid()) {
					SetTooltip("Export already in progress");
				} else if (!cansave) {
					SetTooltip("Can't save to this path");
				}
			}
			SameLine();
			if (Button("Cancel")) CloseCurrentPopup();
			PopID();
			EndPopup();
		}
	}

	void ModelEntity::draw_window_decimation(bool selected) {
		using namespace ImGui;
		if (m_decimated && ui_decimation_window_open() && (m_dec_progress.state < decimation_state::done)) {
			if (Begin("Decimation")) {
				PushID(this);
				draw_select_header(selected);
				draw_decimate_progress(m_dec_progress);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				PopID();
			}
		}
	}

	bool ModelEntity::draw_select_header(bool selected) {
		using namespace ImGui;
		bool r = false;
		PushStyleColor(ImGuiCol_Header, selected ? ImVec4{0.7f, 0.4f, 0.1f, 1} : GetStyle().Colors[ImGuiCol_Button]);
		// hack - selectable is always 'selected' in order to show highlight, it just changes colour
		if (Selectable(name().c_str(), true, 0, {0, GetTextLineHeightWithSpacing()})) {
			auto &sel = ui_selection();
			sel.select_entity = id();
			sel.select_vertex = -1;
			r = true;
		}
		PopStyleColor();
		SetHoveredTooltip("%s", make_name_tooltip().c_str());
		SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
		return r;
	}

	void ModelEntity::spawn_locked_notification() const {
		using namespace ImGui;
		ui_spawn([nametooltip=make_name_tooltip(), name=name()](bool *p_open) {
			if (*p_open) OpenPopup("Model in use");
			if (BeginPopupModal("Model in use", p_open, ImGuiWindowFlags_AlwaysAutoResize)) {
				Selectable(name.c_str(), true, ImGuiSelectableFlags_DontClosePopups, {0, GetTextLineHeightWithSpacing()});
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				SetHoveredTooltip("%s", nametooltip.c_str());
				Text("The operation cannot be performed because\nthe model is currently in use.");
				if (Button("  OK  ")) {
					*p_open = false;
					CloseCurrentPopup();
				}
				EndPopup();
			}
		});
	}

	std::string ModelEntity::make_name_tooltip() const {
		auto s = m_fpath_load.u8string();
		if (m_decimated) {
			s += "\n";
			s += m_dec_uparams.str();
		}
		return s;
	}

	void ModelEntity::invalidate_saliency_vbo() {
		m_saliency_vbo_dirty = true;
		invalidate_scene();
	}

	void ModelEntity::draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) {
		const auto xform0 = transform();
		if (m_pending_load.valid() && m_pending_load.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
			try {
				// is this the best place to do this?
				m_model = m_pending_load.get();
				// gl stuff has to happen on main thread
				m_model->update_vao();
				m_model->update_vbos();
				m_scale = m_model->unit_bound_scale() * 4;
				invalidate_saliency_vbo();
			} catch (std::exception &e) {
				std::cerr << "failed to load model: " << e.what() << std::endl;
			} catch (...) {
				std::cerr << "failed to load model" << std::endl;
			}
		}
		auto &sel = ui_selection();
		const bool selected = sel.select_entity == id();
		if (!dead() && selected) ui_select_model(this);
		if (!dead() && sel.hover_entity == id()) ui_hover_model(this);
		draw_window_models(selected);
		draw_window_decimation(selected);
		if (m_model) {
			if (selected) draw_window_selection();
			if (selected) draw_window_export();
			if (draw_scene) {
				update_vbo();
				// determine color map to apply in shader
				auto &saliency_outputs = m_model->saliency();
				int color_map = 0;
				const bool sal_valid = m_saliency_index < saliency_outputs.size();
				if (m_disp_color_mode == model_color_mode::saliency && sal_valid) color_map = 3;
				if (m_disp_color_mode == model_color_mode::saliency_comparison && sal_valid && m_saliency_baseline_index < saliency_outputs.size()) color_map = 4;
				if (m_disp_color_mode == model_color_mode::vcolor && m_model->prop_vcolor_original().is_valid()) color_map = 1;
				if (m_disp_color_mode == model_color_mode::doncurv) color_map = 3;
				// prepare to draw
				model_draw_params params;
				params.sel = sel;
				params.entity_id = id();
				auto set_cull_faces = [](bool b) {
					if (b) {
						glEnable(GL_CULL_FACE);
						glCullFace(GL_BACK);
					} else {
						glDisable(GL_CULL_FACE);
					}
				};
				// faces
				params.shade_flat = m_shade_flat;
				params.color = {0.6f, 0.6f, 0.5f, 1};
				params.vert_color_map = m_color_faces ? color_map : 0;
				set_cull_faces(m_cull_faces);
				glColorMaski(1, GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_faces) m_model->draw(view * transform(), proj, zfar, params, GL_FILL);
				// edges
				params.shade_flat = false;
				params.shading = 0;
				params.color = {0.03f, 0.03f, 0.03f, 1};
				params.vert_color_map = 0;
				set_cull_faces(m_cull_edges);
				glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_edges) m_model->draw(view * transform(), proj, zfar, params, GL_LINE);
				// verts
				params.color = {0.5f, 0, 0, 1};
				params.vert_color_map = m_color_verts ? color_map : 0;
				params.sel.hover_entity = -1;
				params.sel.select_entity = -1;
				//params.show_samples = sal_valid && (m_color_mode == color_mode::saliency || m_color_mode == color_mode::saliency_comparison);
				glDisable(GL_CULL_FACE);
				glPointSize(m_vert_point_size);
				glColorMaski(1, GL_TRUE, GL_TRUE, GL_FALSE, GL_FALSE);
				if (m_show_verts) m_model->draw(view * transform(), proj, zfar, params, GL_POINT);
			}
		}
		const auto xform1 = transform();
		if (xform0 != xform1) invalidate_scene();
	}

	std::future<saliency_result> ModelEntity::compute_saliency_async(const saliency_user_params &uparams0, saliency_progress &progress) {
		if (!m_model) return {};
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return {};
		}
		saliency_user_params uparams = uparams0;
		if (uparams.auto_contrast) uparams.normal_power = m_model->auto_contrast();
		saliency_mesh_params mparams;
		mparams.mesh = &m_model->trimesh();
		mparams.prop_vertex_area = m_model->prop_vertex_area();
		mparams.prop_doncurv_raw = m_model->prop_doncurv_raw();
		mparams.prop_edge_length = m_model->prop_edge_length();
		// create properties
		mparams.prop_saliency_levels.resize(uparams.levels);
		mparams.mesh->add_property(mparams.prop_curvature);
		mparams.mesh->add_property(mparams.prop_saliency);
		mparams.mesh->add_property(mparams.prop_sampled);
		for (int i = 0; i < uparams.levels; i++) {
			mparams.mesh->add_property(mparams.prop_saliency_levels[i]);
		}
		// cleanup will be run when the result is received from the future
		// can use shared lock for actual computation because only property contents are being used
		// HACK make_shared is an ugly workaround for std::function needing to be copyable
		lock.unlock();
		mparams.cleanup = [=, pprogress=&progress, slock=std::make_shared<std::shared_lock<std::shared_mutex>>(m_modelmtx)](bool r) mutable noexcept {
			// now need to upgrade to unique lock
			slock->unlock();
			std::unique_lock lock(m_modelmtx);
			if (!m_model) return;
			auto &saliency_outputs = m_model->saliency();
			// destroy temp properties
			mparams.mesh->remove_property(mparams.prop_curvature);
			for (int i = 0; i < uparams.levels; i++) {
				mparams.mesh->remove_property(mparams.prop_saliency_levels[i]);
			}
			if (r) {
				// if preview, remove any previous preview results
				saliency_outputs.erase(std::remove_if(saliency_outputs.begin(), saliency_outputs.end(), [](const auto &sd) { return sd.uparams.preview; }), saliency_outputs.end());
				// save user params, progress output and actual saliency mesh property
				model_saliency_data sd;
				sd.uparams = uparams;
				sd.progress = *pprogress;
				sd.prop_saliency = mparams.prop_saliency;
				sd.prop_sampled = mparams.prop_sampled;
				sd.uparams_known = true;
				saliency_outputs.push_back(sd);
				// give focus to this result
				m_saliency_index = saliency_outputs.size() - 1;
				invalidate_saliency_vbo();
			} else {
				// cancelled, destroy saliency property too
				mparams.mesh->remove_property(mparams.prop_saliency);
				mparams.mesh->remove_property(mparams.prop_sampled);
			}
		};
		return green::compute_saliency_async(mparams, uparams, progress);
	}

	std::unique_ptr<ModelEntity> ModelEntity::decimate_async(const decimate_user_params &uparams) {
		// this function isn't const because it needs to lock the mutex - mutable?
		std::shared_lock lock1(m_modelmtx, std::defer_lock);
		if (!lock1.try_lock()) {
			spawn_locked_notification();
			return {};
		}
		auto e = std::make_unique<ModelEntity>();
		std::unique_lock lock2(e->m_modelmtx);
		e->m_fpath_load = m_fpath_load;
		e->m_decimated = true;
		e->m_dec_uparams = uparams;
		e->m_basis_right = m_basis_right;
		e->m_basis_up = m_basis_up;
		e->m_basis_back = m_basis_back;
		e->m_rotation_euler_yxz = m_rotation_euler_yxz;
		e->m_cull_faces = m_cull_faces;
		e->m_cull_edges = m_cull_edges;
		e->m_pending_load = std::async([=, e=e.get(), lock1=std::move(lock1), lock2=std::move(lock2)]() mutable {
			// prepare mesh copy to decimate
			// this is slow enough on big models to be async
			auto &saliency_outputs = m_model->saliency();
			std::vector<model_saliency_data> sdv;
			saliency_prop_t srcprop;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				if (i == m_saliency_index || sd.persistent) {
					srcprop = sd.prop_saliency;
					sdv.push_back(sd);
				} else if (sd.persistent) {
					sdv.push_back(sd);
				}
			}
			auto m = m_model->prepare_decimate(srcprop, sdv);
			// unlock source model
			lock1.unlock();
			// now decimate
			if (!m.decimate(uparams, e->m_dec_progress)) throw std::runtime_error("decimation was cancelled");
			return std::make_unique<Model>(std::move(m));
		});
		return e;
	}

	ModelEntity::~ModelEntity() {

	}

}
