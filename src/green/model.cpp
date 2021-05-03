
#include "model.hpp"

#include <OpenMesh/Core/IO/MeshIO.hh>

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
#include <limits>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/euler_angles.hpp>

#include <cgu/shader.hpp>

#include <imgui.h>

#include "imguiex.hpp"
#include "main.hpp"
#include "curvature.hpp"
#include "neighborhood.hpp"
#include "meshutils.hpp"

#include "model.glsl.hpp"

namespace green {

	float autocontrast_target_entropy = log2f(MeshCache::ncurvbins) * 0.45f;

	ModelBase::ModelBase(const std::filesystem::path &fpath) {
		// face status not really required here
		// NOTE face status currently breaks decimation
		//m_trimesh.request_face_status();
		// TODO what else do we need to ask for and preserve on export?
		OpenMesh::IO::Options readOptions{OpenMesh::IO::OptionBits::Custom | OpenMesh::IO::OptionBits::Fuse};
		readOptions.vattribs |= OpenMesh::AttributeBits::Color;
		readOptions.hattribs |= OpenMesh::AttributeBits::Normal;
		readOptions.hattribs |= OpenMesh::AttributeBits::TexCoordAll;
		
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
		
		// TODO avoid triangulating
		// not triangulating breaks some stuff atm, probably in vertex area
		m_mesh.triangulate();

		if (readOptions.halfedge_has_normal() && m_mesh.has_halfedge_normals()) {
			std::cout << "Found halfedge normals" << std::endl;
		}

		if (readOptions.halfedge_has_texcoord2D() && m_mesh.has_halfedge_texcoords2D()) {
			std::cout << "Found halfedge texcoords2D" << std::endl;
		}

		// copy original vertex colors
		// (because we need to be able to overwrite the actual vertex colors during export)
		if (readOptions.vertex_has_color()) {
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
		const float surfarea = surfaceArea(m_mesh);

		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_mesh);

		// TODO how should we not compute unwanted curvature measures?

		{
			std::cout << "Computing raw don curvature" << std::endl;
			const auto time_curv_start = std::chrono::steady_clock::now();
			curvature_measure curv(m_mesh);
			compute_maxdon(m_mesh, 1, curv);
			const auto time_curv_finish = std::chrono::steady_clock::now();
			std::cout << "Don curvature took " << ((time_curv_finish - time_curv_start) / std::chrono::duration<double>(1.0)) << "s" << std::endl;
			// TODO time autocontrast?
			std::cout << "Computing auto contrast for don curvature" << std::endl;
			m_curv_don = autocontrast(m_mesh, curv, m_prop_vertex_area, true);
		}

		{
			std::cout << "Computing raw mean curvature" << std::endl;
			const auto time_curv_start = std::chrono::steady_clock::now();
			curvature_measure curv(m_mesh);
			compute_mean_curvature(m_mesh, 1, curv);
			const auto time_curv_finish = std::chrono::steady_clock::now();
			std::cout << "Mean curvature took " << ((time_curv_finish - time_curv_start) / std::chrono::duration<double>(1.0)) << "s" << std::endl;
			// TODO time autocontrast?
			std::cout << "Computing auto contrast for mean curvature" << std::endl;
			m_curv_mean = autocontrast(m_mesh, curv, m_prop_vertex_area, false);
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
		} else if (sparams.color_mode == model_color_mode::saliency_comparison) {
			if (!sparams.prop_saliency.is_valid()) {
				sparams.color_mode = model_color_mode::none;
			} else if (!sparams.prop_saliency_baseline.is_valid()) {
				sparams.color_mode = model_color_mode::none;
			} else {
				std::cout << "Colorizing saliency comparison" << std::endl;
			}
		} else if (sparams.color_mode == model_color_mode::doncurv) {
			if (!m_curv_don.prop_curv.is_valid()) {
				sparams.color_mode == model_color_mode::none;
			} else {
				std::cout << "Colorizing don curvature" << std::endl;
			}
		} else {
			std::cout << "Unsupported color mode for model export" << std::endl;
			sparams.color_mode == model_color_mode::none;
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
					const float c = std::pow(m_mesh.property(m_curv_don.prop_curv, *vIt), m_curv_don.contrast * 0.2f);
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
		if (exprops.size()) opts += OpenMesh::IO::OptionBits::Custom;
		if (sparams.binary) opts += OpenMesh::IO::OptionBits::Binary;
		opts.hattribs |= OpenMesh::AttributeBits::Normal;
		opts.hattribs |= OpenMesh::AttributeBits::TexCoordAll;
		if (sparams.color_mode != model_color_mode::none) opts.vattribs |= OpenMesh::AttributeBits::Color;
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
		m.m_mesh.assign(m_mesh, true);
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

	void Model::init_saliency_params(saliency_mesh_params &mparams, const saliency_user_params uparams) {
		mparams.mesh = &m_mesh;
		// set input properties
		mparams.curv = curv(uparams.curv_mode);
		mparams.prop_vertex_area = m_prop_vertex_area;
		mparams.prop_edge_length = m_prop_edge_length;
		// create output properties
		mparams.mesh->add_property(mparams.prop_saliency);
		mparams.mesh->add_property(mparams.prop_sampled);
		// create temp properties
		mparams.mesh->add_property(mparams.prop_curvature);
		mparams.prop_saliency_levels.resize(uparams.levels);
		for (auto &prop : mparams.prop_saliency_levels) {
			mparams.mesh->add_property(prop);
		}
	}

	void Model::cleanup_saliency_params(saliency_mesh_params &mparams, bool success) noexcept {
		assert(mparams.mesh == &m_mesh);
		// destroy temp properties
		mparams.mesh->remove_property(mparams.prop_curvature);
		for (auto &prop : mparams.prop_saliency_levels) {
			mparams.mesh->remove_property(prop);
		}
		if (!success) {
			// cancelled, also destroy output properties
			mparams.mesh->remove_property(mparams.prop_saliency);
			mparams.mesh->remove_property(mparams.prop_sampled);
		}
	}

	void Model::update_vao_verts() {
		bool makevao = false;
		if (!m_vao_verts) makevao = true, m_vao_verts = cgu::gl_object::gen_vertex_array();
		const auto nverts = m_mesh.n_vertices();
		assert(nverts <= std::numeric_limits<decltype(m_vao_nverts)>::max());
		m_vao_nverts = decltype(m_vao_nverts)(nverts);
		glBindVertexArray(m_vao_verts);
		// no index buffer
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 4, GL_FLOAT, false, sizeof(glm::vec4), 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vao_edges() {
		bool makevao = false;
		if (!m_vao_edges) makevao = true, m_vao_edges = cgu::gl_object::gen_vertex_array();
		const auto nedges = m_mesh.n_edges();
		assert(nedges <= std::numeric_limits<decltype(m_vao_nverts)>::max());
		m_vao_nedges = decltype(m_vao_nverts)(nedges);
		const size_t size_ibo = nedges * 2 * sizeof(GLuint);
		glBindVertexArray(m_vao_edges);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo_edges);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_ibo, nullptr, GL_STATIC_DRAW);
		if (nedges) {
			auto pibo = reinterpret_cast<GLuint *>(
				glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, size_ibo, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto eit = m_mesh.edges_begin(); eit != m_mesh.edges_end(); ++eit) {
				auto hh = m_mesh.s_halfedge_handle(*eit, 0);
				*pibo++ = m_mesh.from_vertex_handle(hh).idx();
				*pibo++ = m_mesh.to_vertex_handle(hh).idx();
			}
			glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		}
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 4, GL_FLOAT, false, sizeof(glm::vec4), 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vao_tris() {
		bool makevao = false;
		if (!m_vao_tris) makevao = true, m_vao_tris = cgu::gl_object::gen_vertex_array();
		// calculate triangle count (upper bound)
		const double he_per_face = double(m_mesh.n_halfedges()) / m_mesh.n_faces();
		const double fntris = std::ceil(std::max(1.0, std::ceil(he_per_face - 2)) * m_mesh.n_faces());
		assert(fntris <= double(std::numeric_limits<decltype(m_vao_nverts)>::max()));
		const auto max_ntris = decltype(m_vao_nverts)(fntris);
		// tris counted during triangulation
		m_vao_ntris = 0;
		const size_t size_ibo = size_t(max_ntris) * 3 * sizeof(GLuint);
		glBindVertexArray(m_vao_tris);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo_tris);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_ibo, nullptr, GL_STATIC_DRAW);
		if (max_ntris) {
			auto pibo = reinterpret_cast<GLuint *>(
				glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, size_ibo, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto fit = m_mesh.faces_begin(); fit != m_mesh.faces_end(); ++fit) {
				auto hit = m_mesh.cfh_iter(*fit);
				const auto hit0 = hit++;
				auto hit1 = hit++;
				// triangulate
				for (; hit.is_valid(); ++hit1, ++hit) {
					// record actual tri count
					assert(m_vao_ntris < max_ntris);
					m_vao_ntris++;
					*pibo++ = GLuint(hit0->idx());
					*pibo++ = GLuint(hit1->idx());
					*pibo++ = GLuint(hit->idx());
				}
			}
			glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		}
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_halfedges);
			glEnableVertexAttribArray(0);
			glVertexAttribIPointer(0, 1, GL_INT, sizeof(halfedge_vbo_props), (void *) offsetof(halfedge_vbo_props, vid));
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, false, sizeof(halfedge_vbo_props), (void *) offsetof(halfedge_vbo_props, norm));
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 2, GL_FLOAT, false, sizeof(halfedge_vbo_props), (void *) offsetof(halfedge_vbo_props, uv));
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vao_boundaries() {
		// like edges, but only boundary ones (next to a topological hole)
		bool makevao = false;
		if (!m_vao_boundaries) makevao = true, m_vao_boundaries = cgu::gl_object::gen_vertex_array();
		const auto nedges = m_mesh.n_edges();
		assert(nedges <= std::numeric_limits<decltype(m_vao_nboundaries)>::max());
		// boundary edges counted while filling ibo
		m_vao_nboundaries = 0;
		// TODO this will be generally well over-allocated atm
		const size_t size_ibo = nedges * 2 * sizeof(GLuint);
		glBindVertexArray(m_vao_boundaries);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo_boundaries);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_ibo, nullptr, GL_STATIC_DRAW);
		if (nedges) {
			auto pibo = reinterpret_cast<GLuint *>(
				glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, size_ibo, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto eit = m_mesh.edges_begin(); eit != m_mesh.edges_end(); ++eit) {
				auto hh = m_mesh.s_halfedge_handle(*eit, 0);
				auto ohh = m_mesh.opposite_halfedge_handle(hh);
				if (m_mesh.is_boundary(hh) || m_mesh.is_boundary(ohh)) {
					// only add to ibo if a boundary edge
					m_vao_nboundaries++;
					*pibo++ = m_mesh.from_vertex_handle(hh).idx();
					*pibo++ = m_mesh.to_vertex_handle(hh).idx();
				}
			}
			glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		}
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 4, GL_FLOAT, false, sizeof(glm::vec4), 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vaos() {
		if (!m_ibo_edges) m_ibo_edges = cgu::gl_object::gen_buffer();
		if (!m_ibo_tris) m_ibo_tris = cgu::gl_object::gen_buffer();
		if (!m_ibo_boundaries) m_ibo_boundaries = cgu::gl_object::gen_buffer();
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		if (!m_vbo_halfedges) m_vbo_halfedges = cgu::gl_object::gen_buffer();
		if (!m_tex_col) m_tex_col = cgu::gl_object::gen_texture();
		if (!m_tex_pos) m_tex_pos = cgu::gl_object::gen_texture();
		// actually create buffer objects for binding to textures
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		// setup buffer textures
		glBindTexture(GL_TEXTURE_BUFFER, m_tex_col);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, m_vbo_col);
		// TODO three-component formats can only be used in GL 4.0 or with ARB_texture_buffer_object_rgb32
		glBindTexture(GL_TEXTURE_BUFFER, m_tex_pos);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, m_vbo_pos);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
		update_vao_verts();
		update_vao_edges();
		update_vao_tris();
		update_vao_boundaries();
	}

	void Model::update_vbos() {
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		if (!m_vbo_halfedges) m_vbo_halfedges = cgu::gl_object::gen_buffer();
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBufferData(GL_ARRAY_BUFFER, m_mesh.n_vertices() * sizeof(glm::vec3), m_mesh.points(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		glBufferData(GL_ARRAY_BUFFER, m_mesh.n_vertices() * sizeof(glm::vec4), nullptr, GL_STATIC_DRAW);
		const size_t size_vbo_halfedges = m_mesh.n_halfedges() * sizeof(halfedge_vbo_props);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_halfedges);
		glBufferData(GL_ARRAY_BUFFER, size_vbo_halfedges, nullptr, GL_STATIC_DRAW);
		m_force_shade_flat = !m_mesh.has_halfedge_normals() && !m_mesh.has_vertex_normals();
		if (m_mesh.n_halfedges()) {
			auto pvbo = reinterpret_cast<halfedge_vbo_props *>(
				glMapBufferRange(GL_ARRAY_BUFFER, 0, size_vbo_halfedges, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto hit = m_mesh.halfedges_begin(); hit != m_mesh.halfedges_end(); ++hit) {
				auto vh = m_mesh.to_vertex_handle(*hit);
				pvbo->vid = vh.idx();
				if (m_mesh.has_halfedge_normals()) {
					pvbo->norm = om2glm(m_mesh.normal(*hit));
				} else if (m_mesh.has_vertex_normals()) {
					pvbo->norm = om2glm(m_mesh.normal(vh));
				}
				if (m_mesh.has_halfedge_texcoords2D()) {
					pvbo->uv = om2glm(m_mesh.texcoord2D(*hit));
				}
				pvbo++;
			}
			glUnmapBuffer(GL_ARRAY_BUFFER);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	bool Model::update_color(const model_color_params &cparams, model_saliency_errors *err) {
		model_saliency_errors err2{};
		if (!err) err = &err2;
		if (cparams.color_mode == model_color_mode::saliency && !cparams.prop_saliency.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::vcolor && !m_prop_vcolor_original.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::saliency_comparison && !cparams.prop_saliency.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::saliency_comparison && !cparams.prop_saliency_baseline.is_valid()) return false;
		if (cparams.color_mode == model_color_mode::doncurv && !m_curv_don.prop_curv.is_valid()) return false;
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
				const float c = std::pow(m_mesh.property(m_curv_don.prop_curv, v), m_curv_don.contrast * 0.2f);
				data[i] = glm::vec4(c, 0, 0, 0);
			}
		} else {
			// uvs/texture are handled in the shader directly
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		return true;
	}

	void Model::draw(GLenum polymode, bool boundaries) const {
		if (polymode == GL_POINT && !boundaries) {
			if (!m_vao_verts) return;
			glBindVertexArray(m_vao_verts);
			glDrawArrays(GL_POINTS, 0, m_vao_nverts);
		} else if (polymode == GL_LINE && !boundaries) {
			if (!m_vao_edges) return;
			glBindVertexArray(m_vao_edges);
			// must be uint (not int)
			glDrawElements(GL_LINES, m_vao_nedges * 2, GL_UNSIGNED_INT, 0);
		} else if (polymode == GL_LINE && boundaries) {
			if (!m_vao_boundaries) return;
			glBindVertexArray(m_vao_boundaries);
			// must be uint (not int)
			glDrawElements(GL_LINES, m_vao_nboundaries * 2, GL_UNSIGNED_INT, 0);
		} else if (polymode == GL_FILL && !boundaries) {
			if (!m_vao_tris) return;
			glBindVertexArray(m_vao_tris);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			// must be uint (not int)
			glDrawElements(GL_TRIANGLES, m_vao_ntris * 3, GL_UNSIGNED_INT, 0);
		}
		glBindVertexArray(0);
	}

	void Model::draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params) const {

		const GLenum polymode = params.polymode;
		const bool shade_flat = m_force_shade_flat || params.shade_flat;

		// note: need custom depth env
		static GLuint prog_verts = 0;
		static GLuint prog_edges = 0;
		static GLuint prog_tris_smooth = 0;
		static GLuint prog_tris_flat = 0;
		if (polymode == GL_POINT && !prog_verts) {
			prog_verts = cgu::make_shader_program(
				"green::model.verts",
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{"#define MODEL_VERTS\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}
		if (polymode == GL_LINE && !prog_edges) {
			prog_edges = cgu::make_shader_program(
				"green::model.edges",
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{"#define MODEL_EDGES\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}
		if (polymode == GL_FILL && !shade_flat && !prog_tris_smooth) {
			prog_tris_smooth = cgu::make_shader_program(
				"green::model.tris.smooth",
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{"#define MODEL_TRIS\n#define SHADE_SMOOTH\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}
		if (polymode == GL_FILL && shade_flat && !prog_tris_flat) {
			prog_tris_flat = cgu::make_shader_program(
				"green::model.tris.flat",
				"330 core",
				{GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER},
				{"#define MODEL_TRIS\n#define SHADE_FLAT\n", cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}

		GLuint prog = 0;
		if (polymode == GL_POINT) prog = prog_verts;
		if (polymode == GL_LINE) prog = prog_edges;
		if (polymode == GL_FILL) prog = shade_flat ? prog_tris_flat : prog_tris_smooth;
		if (!prog) return;

		glUseProgram(prog);
		
		if (polymode == GL_FILL) {
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_BUFFER, m_tex_pos);
			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_BUFFER, m_tex_col);
			// texunit 0 is for user-supplied texture
			glUniform1i(glGetUniformLocation(prog, "u_sampler_col"), 0);
			glUniform1i(glGetUniformLocation(prog, "u_buf_pos"), 1);
			glUniform1i(glGetUniformLocation(prog, "u_buf_col"), 2);
		}

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
		glUniform1i(glGetUniformLocation(prog, "u_color_map"), params.vert_color_map);
		glUniform1i(glGetUniformLocation(prog, "u_show_samples"), params.show_samples);

		draw(polymode, params.boundaries);

		if (polymode == GL_FILL) {
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_BUFFER, 0);
			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_BUFFER, 0);
		}
	}
}
