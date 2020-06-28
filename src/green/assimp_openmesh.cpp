
#include <cassert>
#include <climits>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/IO/IOManager.hh>

#include "assimp_openmesh.hpp"

namespace stdfs = std::filesystem;

namespace {

	inline OpenMesh::Vec3f ai2om(const aiVector3D &v) {
		OpenMesh::Vec3f vv;
		vv[0] = v.x;
		vv[1] = v.y;
		vv[2] = v.z;
		return vv;
	}

	inline OpenMesh::Vec4f ai2om(const aiColor4D &v) {
		OpenMesh::Vec4f vv;
		vv[0] = v.r;
		vv[1] = v.g;
		vv[2] = v.b;
		vv[3] = v.a;
		return vv;
	}

}

namespace green {

	AssimpReader::AssimpReader() {
		// commented - dont pollute cli output
		//omerr() << "[Assimp] : registering openmesh reader module" << std::endl;
		OpenMesh::IO::IOManager().register_module(this);
	}

	bool AssimpReader::read(const std::filesystem::path &_filename, OpenMesh::IO::BaseImporter &_bi, OpenMesh::IO::Options &_opt) {
		omerr() << "[Assimp] : loading " << _filename.u8string() << std::endl;
		// NOTE assimp appears to support utf8 filenames; see assimp-5.0.1\code\Common\DefaultIOSystem.cpp
		Assimp::Importer importer;
		//importer.SetProgressHandler
		// we need to fuse colocated vertices so the topology works
		// TODO what is the best way of achieving this?
		// what is the runtime cost of using assimp's join identical?
		// remove components that could prevent colocated vertices from being identical for joining purposes
		// NOTE we want to load vertex colors
		importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS, aiComponent_TEXCOORDS | aiComponent_NORMALS | aiComponent_TANGENTS_AND_BITANGENTS);
		// scene is owned by the importer
		const aiScene *scene = importer.ReadFile(
			_filename.u8string(),
			aiProcess_PreTransformVertices
			| aiProcess_Triangulate
			| aiProcess_JoinIdenticalVertices
			| aiProcess_RemoveComponent
		);
		if (!scene) {
			omerr() << "[Assimp] : failed to load " << _filename.u8string() << ": " << importer.GetErrorString() << std::endl;
			return false;
		}
		OpenMesh::IO::BaseImporter::VHandles vhandles;
		// combine all meshes into one
		for (unsigned i = 0; i < scene->mNumMeshes; i++) {
			const aiMesh *mesh = scene->mMeshes[i];
			// add vertices for this mesh
			// TODO individual mesh transforms
			OpenMesh::VertexHandle vh0;
			for (unsigned j = 0; j < mesh->mNumVertices; j++) {
				auto vh = _bi.add_vertex(ai2om(mesh->mVertices[j]));
				if (j == 0) vh0 = vh;
				if (mesh->HasNormals()) {
					_bi.set_normal(vh, ai2om(mesh->mNormals[j]));
					_opt.set(OpenMesh::IO::Options::VertexNormal);
				}
				if (mesh->HasTangentsAndBitangents()) {
					// openmesh doesnt care?
				}
				if (mesh->HasTextureCoords(0)) {
					_bi.set_texcoord(vh, ai2om(mesh->mTextureCoords[0][j]));
					_opt.set(OpenMesh::IO::Options::VertexTexCoord);
				}
				if (mesh->HasVertexColors(0)) {
					_bi.set_color(vh, ai2om(mesh->mColors[0][j]));
					_opt.set(OpenMesh::IO::Options::VertexColor);
				}
			}
			// add face indices for this mesh
			// TODO join same-position vertices (and average their normals or something)
			for (unsigned j = 0; j < mesh->mNumFaces; j++) {
				vhandles.clear();
				const aiFace &face = mesh->mFaces[j];
				for (unsigned k = 0; k < face.mNumIndices; k++) {
					const unsigned uvi = face.mIndices[k];
					// TODO proper overflow check
					vhandles.push_back(OpenMesh::VertexHandle{vh0.idx() + int(uvi)});
				}
				_bi.add_face(vhandles);
			}
		}
		return true;
	}

	bool AssimpReader::read(std::istream &_is, OpenMesh::IO::BaseImporter &_bi, OpenMesh::IO::Options &_opt) {
		// TODO assimp doesnt expose a good way to implement this
		// reading entire file into buffer to use load from memory is undesirable
		// best way might be via a named pipe / fifo
		throw std::logic_error("not implemented");
	}

	bool AssimpReader::can_u_read(const std::filesystem::path &_filename) const {
		// assume utf8
		auto ext = _filename.extension();
		for (const char *e : {".om", ".OM", ".ply", ".PLY", ".obj", ".OBJ", ".off", ".OFF", ".stl", ".STL"}) {
			// skip formats that openmesh can already load
			if (ext == e) return false;
		}
		// TODO check if assimp knows the extension?
		// this function is seemingly supposed to inspect the file contents
		// but we don't have a way to do that with assimp
		// we are supposed to check if we can open the file, see PLYReader.cc
		std::ifstream ifs(_filename);
		return ifs.is_open();
	}

	AssimpReaderInit::AssimpReaderInit() {
		static AssimpReader reader;
	}

}
