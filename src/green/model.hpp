
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

#include <string>
#include <string_view>
#include <filesystem>
#include <utility>
#include <memory>
#include <algorithm>
#include <future>
#include <shared_mutex>

#include <cgu/opengl.hpp>

// ensure assimp openmesh modules are registered via static init
#include "assimp_openmesh.hpp"

#include "meshutils.hpp"
#include "entity.hpp"
#include "decimate.hpp"

namespace green {

	enum class model_color_mode : unsigned char {
		none,
		vcolor,
		saliency,
		saliency_comparison, 
		doncurv,
		uv,
		checkerboard
	};

	struct model_draw_params {
		GLenum polymode = GL_FILL;
		entity_selection sel;
		glm::vec4 color{0.6f, 0.6f, 0.5f, 1};
		float shading = 0.9f;
		int entity_id = -1;
		int vert_color_map = 0;
		bool show_samples = false;
		bool shade_flat = false;
		bool boundaries = false;
	};

	struct model_color_params {
		saliency_prop_t prop_saliency{};
		saliency_prop_t prop_saliency_baseline{};
		model_color_mode color_mode = model_color_mode::none;
		float error_scale = 1;
	};

	struct model_save_params : model_color_params {
		bool original_vids = false;
		bool binary = true;
	};

	struct model_saliency_errors {
		float min = 0;
		float max = 0;
		float rms = 0;
	};

	struct model_saliency_data {
		std::string filename;
		std::string dispname;
		std::string expname;
		saliency_user_params uparams;
		saliency_progress progress;
		saliency_prop_t prop_saliency{};
		OpenMesh::VPropHandleT<unsigned char> prop_sampled{};
		bool decimated = false;
		bool persistent = false;
		bool uparams_known = false;
		bool should_export = false;

		std::string str() const {
			std::string s;
			if (decimated) {
				s += "<decimated> ";
			} else if (uparams.preview) {
				s += "<preview> ";
			}
			if (filename.empty()) {
				s += dispname.empty() ? uparams.str() : dispname;
			} else {
				s += filename;
				if (!dispname.empty()) {
					s += "?";
					s += dispname;
				}
			}
			return s;
		}

		explicit operator std::string() const {
			return str();
		}

		std::string export_propname(int i) const {
			std::string s = "quality";
			s += !expname.empty() ? expname : (!dispname.empty() ? dispname : "");
			s += '@';
			s += std::to_string(i);
			if (uparams_known) {
				s += '[';
				if (decimated) s += 'D';
				if (s.back() != '[') s += ';';
				s += uparams.str(true);
				s += ']';
			}
			return s;
		}
	};

	// this class exists to separate the minimal loading of a model from file
	// from any computations that we do once the model is loaded (vertex areas, auto contrast).
	// one good reason for this is to exclude only the initial load from timing.
	class ModelBase {
	protected:
		// FIXME openmesh types are not movable
		PolyMesh m_mesh;
		OpenMesh::VPropHandleT<PolyMesh::Color> m_prop_vcolor_original;
		OpenMesh::VPropHandleT<int> m_prop_vid_original;
		std::vector<model_saliency_data> m_saliency;

	public:
		ModelBase() {}

		ModelBase(const ModelBase &) = delete;
		ModelBase & operator=(const ModelBase &) = delete;

		// FIXME noexcept
		ModelBase(ModelBase &&) = default;
		ModelBase & operator=(ModelBase &&) = default;

		ModelBase(const std::filesystem::path &fpath);

		const PolyMesh & mesh() const {
			return m_mesh;
		}

		PolyMesh & mesh() {
			return m_mesh;
		}

		bool has_texcoords2d() const {
			return m_mesh.has_halfedge_texcoords2D();
		}

		std::vector<model_saliency_data> & saliency() {
			return m_saliency;
		}

		const std::vector<model_saliency_data> & saliency() const {
			return m_saliency;
		}

		// search by display name (first found) or @index (into current saliency properties).
		// note that for @index will only correlate with property names in files if there have been no removals.
		std::vector<model_saliency_data>::const_iterator find_saliency(std::string_view name) const;

		std::vector<model_saliency_data>::iterator find_saliency(std::string_view name) {
			const auto &cthis = *this;
			auto cit = cthis.find_saliency(name);
			return m_saliency.begin() + (cit - m_saliency.cbegin());
		}

		OpenMesh::VPropHandleT<PolyMesh::Color> prop_vcolor_original() const {
			return m_prop_vcolor_original;
		}

	};

	class Model : public ModelBase {
	private:
		OpenMesh::VPropHandleT<float> m_prop_vertex_area;
		OpenMesh::VPropHandleT<float> m_prop_doncurv_raw;
		OpenMesh::EPropHandleT<float> m_prop_edge_length;
		saliency_prop_t m_prop_sal_dec;

		glm::vec3 m_bound_min{9001e19f}, m_bound_max{-9001e19f};
		float m_auto_contrast = 1;

		// glsl texelFetch only takes int
		GLint m_vao_nverts = 0, m_vao_nedges = 0, m_vao_ntris = 0, m_vao_nboundaries = 0;
		cgu::gl_object m_vao_verts, m_vao_edges, m_vao_tris, m_vao_boundaries;
		cgu::gl_object m_ibo_edges, m_ibo_tris, m_ibo_boundaries;
		cgu::gl_object m_vbo_pos, m_vbo_col, m_vbo_halfedges;
		cgu::gl_object m_tex_pos, m_tex_col;
		bool m_force_shade_flat = false;

	public:
		Model() {}

		Model(const Model &) = delete;
		Model & operator=(const Model &) = delete;

		// FIXME noexcept
		Model(Model &&) = default;
		Model & operator=(Model &&) = default;

		Model(ModelBase &&base);

		Model(const std::filesystem::path &fpath);

		void save(const std::filesystem::path &fpath, const model_save_params &sparams);

		Model prepare_decimate(saliency_prop_t prop_saliency, const std::vector<model_saliency_data> &sdv) const;

		bool decimate(const decimate_user_params &uparams, decimate_progress &progress);

		OpenMesh::VPropHandleT<float> prop_vertex_area() const {
			return m_prop_vertex_area;
		}

		OpenMesh::VPropHandleT<float> prop_doncurv_raw() const {
			return m_prop_doncurv_raw;
		}

		OpenMesh::EPropHandleT<float> prop_edge_length() const {
			return m_prop_edge_length;
		}

		saliency_prop_t prop_sal_dec() const {
			return m_prop_sal_dec;
		}

		glm::vec3 bound_min() const {
			return m_bound_min;
		}

		glm::vec3 bound_max() const {
			return m_bound_max;
		}

		glm::vec3 bound_size() const {
			return m_bound_max - m_bound_min;
		}

		glm::vec3 bound_center() const {
			return m_bound_min + bound_size() * 0.5f;
		}

		float unit_bound_scale() const {
			glm::vec3 x = 1.f / bound_size();
			float xmin = x[0];
			for (float f : {x.x, x.y, x.z}) {
				xmin = std::min(f, xmin);
			}
			return xmin;
		}

		float auto_contrast() const {
			return m_auto_contrast;
		}

		// TODO does this still need exposing?
		GLuint vbo_color() const {
			return m_vbo_col;
		}

		GLuint vao_nverts() const {
			return m_vao_nverts;
		}

		void update_vaos();

		void update_vbos();

		bool update_color(const model_color_params &cparams, model_saliency_errors * = nullptr);

		void draw(GLenum polymode = GL_FILL, bool boundaries = false) const;

		void draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params) const;

	private:
		void update_vao_verts();
		void update_vao_edges();
		void update_vao_tris();
		void update_vao_boundaries();

		struct halfedge_vbo_props {
			// glsl texelFetch only takes int
			GLint vid = 0;
			glm::vec3 norm{};
			glm::vec2 uv{};
		};
	};

}

#endif
