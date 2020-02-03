
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

#include <filesystem>
#include <utility>
#include <algorithm>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

// ensure assimp openmesh modules are registered via static init
#include "assimp_openmesh.hpp"

#include <cgu/opengl.hpp>

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

	class Model {
	private:
		Model(const Model &) = delete;
		Model & operator=(const Model &) = delete;

		TriMesh m_trimesh;

		glm::vec3 m_bound_min{9001e19f}, m_bound_max{-9001e19f};

		GLuint m_vao_ntris = 0;
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

		void draw() const;

		void draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar) const;

	};

}

#endif
