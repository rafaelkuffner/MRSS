
#pragma once

#ifndef GREEN_MODEL_HPP
#define GREEN_MODEL_HPP

#include <filesystem>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

// ensure assimp openmesh modules are registered via static init
#include "assimp_openmesh.hpp"

#include <cgu/opengl.hpp>

namespace green {

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

		void update_vao();

		void update_vbos();

		void draw() const;

	};

}

#endif
