
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

#include <filesystem>
#include <utility>
#include <memory>
#include <algorithm>
#include <future>

#include <cgu/opengl.hpp>

// ensure assimp openmesh modules are registered via static init
#include "assimp_openmesh.hpp"

#include "meshutils.hpp"
#include "entity.hpp"

namespace green {

	struct model_draw_params {
		entity_selection sel;
		glm::vec4 color{0.6f, 0.6f, 0.5f, 1};
		glm::vec4 color_hover{1, 1, 0, 1};
		glm::vec4 color_select{1, 0.5f, 0, 1};
		float shading = 0.9f;
		int entity_id = -1;
		bool use_vert_color = false;
	};

	struct model_saliency_data {
		saliency_user_params uparams;
		saliency_progress progress;
		OpenMesh::VPropHandleT<float> prop_saliency;
	};

	class Model {
	private:
		Model(const Model &) = delete;
		Model & operator=(const Model &) = delete;

		TriMesh m_trimesh;
		OpenMesh::VPropHandleT<float> m_prop_vertex_area;
		OpenMesh::EPropHandleT<float> m_prop_edge_length;

		glm::vec3 m_bound_min{9001e19f}, m_bound_max{-9001e19f};

		GLuint m_vao_ntris = 0, m_vao_nverts = 0;
		cgu::gl_object m_vao;
		cgu::gl_object m_ibo;
		cgu::gl_object m_vbo_pos, m_vbo_norm;

	public:
		Model(const std::filesystem::path &fpath);

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

		void update_vao();

		void update_vbos();

		void draw(GLenum polymode = GL_FILL) const;

		void draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params, GLenum polymode = GL_FILL) const;

	};

	class ModelEntity : public Entity {
	private:
		std::filesystem::path m_fpath;
		std::unique_ptr<Model> m_model;
		std::future<std::unique_ptr<Model>> m_pending;
		std::vector<model_saliency_data> m_saliency_outputs;
		int m_saliency_index = 0;

		float m_scale = 1;
		glm::vec3 m_translation{0};

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

		int basis_right = 0, basis_up = 2, basis_back = 4;

		bool m_show_faces = true;
		bool m_show_edges = false;
		bool m_show_verts = false;
		bool m_show_saliency = true;
		bool m_dead = false;

	public:
		ModelEntity();

		const Model * model() const {
			return m_model.get();
		}

		Model * model() {
			return m_model.get();
		}

		void load(const std::filesystem::path &fpath);

		virtual std::string name() const override {
			// TODO unicode...
			return m_fpath.filename().string();
		}

		virtual void move_by(const glm::vec3 &d) override {
			m_translation += d;
		}

		virtual bool dead() const override {
			return m_dead;
		}

		virtual glm::mat4 transform() const override;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar) override;

		virtual std::future<saliency_result> compute_saliency_async(const saliency_user_params &uparams, saliency_progress &progress) override;

		virtual ~ModelEntity();
	};
}

#endif
