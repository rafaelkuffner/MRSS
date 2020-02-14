
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

#include <filesystem>
#include <utility>
#include <memory>
#include <algorithm>
#include <future>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <cgu/opengl.hpp>

// ensure assimp openmesh modules are registered via static init
#include "assimp_openmesh.hpp"

#include "entity.hpp"

namespace green {

	inline OpenMesh::Vec3f glm2om(const glm::vec3 &v) {
		OpenMesh::Vec3f vv;
		vv[0] = v.x;
		vv[1] = v.y;
		vv[2] = v.z;
		return vv;
	}

	inline glm::vec3 om2glm(const OpenMesh::Vec3f &v) {
		glm::vec3 vv;
		vv.x = v[0];
		vv.y = v[1];
		vv.z = v[2];
		return vv;
	}

	struct DefaultMeshTraitsLCE
	{
		// OpenMesh default
		typedef OpenMesh::Vec3f  Point;

		// OpenMesh default
		typedef OpenMesh::Vec3f  Normal;

		// OpenMesh default
		typedef float  TexCoord1D;

		// use higher precision for texture coordinates
		typedef OpenMesh::Vec2d  TexCoord2D;

		// OpenMesh default
		typedef OpenMesh::Vec3f  TexCoord3D;

		// OpenMesh default
		typedef int TextureIndex;

		// OpenMesh default
		typedef OpenMesh::Vec3uc Color;

#ifndef DOXY_IGNORE_THIS
		VertexTraits    {};
		HalfedgeTraits  {};
		EdgeTraits      {};
		FaceTraits      {};
#endif

		VertexAttributes(0);
		HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
		EdgeAttributes(0);
		FaceAttributes(0);
	};

	struct TriMesh : public OpenMesh::PolyMesh_ArrayKernelT<DefaultMeshTraitsLCE>
	{
	};

	struct model_draw_params {
		selection sel;
		glm::vec4 color{0.6f, 0.6f, 0.5f, 1};
		glm::vec4 color_hover{1, 1, 0, 1};
		glm::vec4 color_select{1, 0.5f, 0, 1};
		float shading = 0.9f;
		int entity_id = -1;
		bool use_vert_color = false;
	};

	class Model {
	private:
		Model(const Model &) = delete;
		Model & operator=(const Model &) = delete;

		TriMesh m_trimesh;

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

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, selection &sel) override;

		virtual std::future<bool> compute_saliency_async(const saliency_params &params, saliency_progress &progress) override;

		virtual ~ModelEntity();
	};
}

#endif
