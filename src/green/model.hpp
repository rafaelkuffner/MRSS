
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

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

	struct model_draw_params {
		entity_selection sel;
		glm::vec4 color{0.6f, 0.6f, 0.5f, 1};
		float shading = 0.9f;
		int entity_id = -1;
		int vert_color_map = 0;
		bool show_samples = false;
	};

	struct model_saliency_data {
		std::string filename;
		std::string propname = "quality";
		saliency_user_params uparams;
		saliency_progress progress;
		OpenMesh::VPropHandleT<float> prop_saliency{};
		OpenMesh::VPropHandleT<unsigned char> prop_sampled{};
		bool decimated = false;

		std::string str() const {
			// actual property name should be prefixed with 'quality'
			// this is then stripped for the display name
			std::string pdname = propname.substr(0, 7) == "quality" ? propname.substr(7) : propname;
			std::string s;
			if (decimated) {
				s += "<decimated> ";
			} else if (uparams.preview) {
				s += "<preview> ";
			}
			if (filename.empty()) {
				s += pdname.empty() ? uparams.str() : pdname;
			} else {
				s += filename;
				if (!pdname.empty()) {
					s += "?";
					s += pdname;
				}
			}
			return s;
		}

		explicit operator std::string() const {
			return str();
		}
	};

	class Model {
	private:
		Model(const Model &) = delete;
		Model & operator=(const Model &) = delete;

		TriMesh m_trimesh;
		OpenMesh::VPropHandleT<float> m_prop_vertex_area;
		OpenMesh::EPropHandleT<float> m_prop_edge_length;
		OpenMesh::VPropHandleT<TriMesh::Color> m_prop_vcolor_original;
		std::vector<model_saliency_data> m_original_saliency;

		glm::vec3 m_bound_min{9001e19f}, m_bound_max{-9001e19f};

		GLuint m_vao_ntris = 0, m_vao_nverts = 0;
		cgu::gl_object m_vao;
		cgu::gl_object m_ibo;
		cgu::gl_object m_vbo_pos, m_vbo_norm, m_vbo_col;

	public:
		Model() {}

		Model(Model &&other) = default;
		Model & operator=(Model &&other) = default;

		Model(const std::filesystem::path &fpath);

		void save(const std::filesystem::path &fpath, OpenMesh::VPropHandleT<float> prop_saliency, bool binary);

		Model prepare_decimate(const model_saliency_data &sd) const;

		bool decimate(const decimate_user_params &uparams, decimate_progress &progress);

		const TriMesh & trimesh() const {
			return m_trimesh;
		}

		TriMesh & trimesh() {
			return m_trimesh;
		}

		OpenMesh::VPropHandleT<float> prop_vertex_area() const {
			return m_prop_vertex_area;
		}

		OpenMesh::EPropHandleT<float> prop_edge_length() const {
			return m_prop_edge_length;
		}

		std::vector<model_saliency_data> & original_saliency() {
			return m_original_saliency;
		}

		const std::vector<model_saliency_data> & original_saliency() const {
			return m_original_saliency;
		}

		OpenMesh::VPropHandleT<TriMesh::Color> prop_vcolor_original() const {
			return m_prop_vcolor_original;
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

		GLuint vbo_color() const {
			return m_vbo_col;
		}

		GLuint vao_nverts() const {
			return m_vao_nverts;
		}

		void update_vao();

		void update_vbos();

		void draw(GLenum polymode = GL_FILL) const;

		void draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params, GLenum polymode = GL_FILL) const;

	};

	class ModelEntity : public Entity {
	private:
		// unique: adding/removing properties, vertices etc
		// shared: reading/writing contents of existing properties
		std::shared_mutex m_modelmtx;
		std::unique_ptr<Model> m_model;

		std::filesystem::path m_fpath_load;
		std::future<std::unique_ptr<Model>> m_pending_load;
		decimate_user_params m_dec_uparams;
		decimate_progress m_dec_progress;
		bool m_decimated = false;

		std::filesystem::path m_fpath_save;
		bool m_save_binary = true;
		bool m_save_ok = false;
		std::future<bool> m_pending_save;

		std::vector<model_saliency_data> m_saliency_outputs;
		int m_saliency_index = 0;
		int m_saliency_baseline_index = 0;
		bool m_saliency_vbo_dirty = false;

		struct saliency_errors {
			float min = 0;
			float max = 0;
			float rms = 0;
		};

		saliency_errors m_saliency_errors;
		float m_saliency_error_scale = 1;

		float m_scale = 1;
		glm::vec3 m_translation{0};
		// TODO store as quat, only use euler for editing?
		glm::vec3 m_rotation_euler_yxz{0};

		struct basis_vector {
			glm::vec3 v;
			const char *name;
		};

		basis_vector m_basis_vectors[6]{
			{{+1, 0, 0}, "+X"},
			{{-1, 0, 0}, "-X"},
			{{0, +1, 0}, "+Y"},
			{{0, -1, 0}, "-Y"},
			{{0, 0, +1}, "+Z"},
			{{0, 0, -1}, "-Z"}
		};

		int m_basis_right = 0, m_basis_up = 2, m_basis_back = 4;

		enum class color_mode : unsigned char {
			none, vcolor, saliency, saliency_comparison
		};

		int m_vert_point_size = 3;
		color_mode m_color_mode = color_mode::saliency;
		bool m_color_faces = true;
		bool m_color_verts = false;
		bool m_show_faces = true;
		bool m_show_edges = false;
		bool m_show_verts = false;
		bool m_cull_faces = true;
		bool m_cull_edges = true;
		
		bool m_try_export = false;
		bool m_try_kill = false;
		bool m_dead = false;

		void update_vbo();
		void draw_window_models(bool selected);
		void draw_window_selection();
		void draw_window_export();
		void draw_window_decimation(bool selected);
		bool draw_select_header(bool selected);
		void spawn_locked_notification() const;
		std::string make_name_tooltip() const;

	public:
		ModelEntity();

		const Model * model() const {
			return m_model.get();
		}

		Model * model() {
			return m_model.get();
		}

		void load(const std::filesystem::path &fpath);

		void save(const std::filesystem::path &fpath, bool binary);

		void load(const std::filesystem::path &fpath, Model m);

		void load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress);

		void try_export() {
			m_try_export = true;
		}

		bool show_entity_outline() const {
			return m_show_faces;
		}

		bool decimated() const {
			return m_decimated;
		}

		const model_saliency_data * selected_saliency() const {
			if (m_saliency_index >= m_saliency_outputs.size()) return nullptr;
			return &m_saliency_outputs[m_saliency_index];
		}

		virtual std::string name() const override {
			std::string s = m_fpath_load.filename().u8string();
			if (m_decimated) {
				s += " [dec ";
				s += m_dec_uparams.str();
				s += "]";
			}
			return s;
		}

		virtual void move_by(const glm::vec3 &d) override {
			m_translation += d;
		}

		virtual void try_kill() override {
			m_try_kill = true;
		}

		virtual bool dead() const override {
			return m_dead;
		}

		virtual glm::mat4 transform() const override;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar) override;

		virtual std::future<saliency_result> compute_saliency_async(const saliency_user_params &uparams, saliency_progress &progress) override;

		std::unique_ptr<ModelEntity> decimate_async(const decimate_user_params &uparams);

		virtual ~ModelEntity();
	};
}

#endif
