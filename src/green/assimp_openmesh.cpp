
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

#define OMLOG_SOURCE AssimpReader

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

	bool AssimpReader::read(const std::filesystem::path &_filename, OpenMesh::IO::BaseImporter &_bi) {
		OMLOG_INFO << "loading " << _filename.u8string();
		// NOTE assimp appears to support utf8 filenames; see assimp-5.0.1\code\Common\DefaultIOSystem.cpp
		Assimp::Importer importer;
		//importer.SetProgressHandler();
		// we dont care about tangents atm (openmesh doesnt know about them)
		importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS, aiComponent_TANGENTS_AND_BITANGENTS);
		// scene is owned by the importer
		// join identical is slow but necessary, else assimp will create new vertices for each face
		const aiScene *scene = importer.ReadFile(
			_filename.u8string(),
			aiProcess_PreTransformVertices
			| aiProcess_JoinIdenticalVertices
			| aiProcess_RemoveComponent
		);
		if (!scene) {
			OMLOG_ERROR << "failed to load " << _filename.u8string() << ": " << importer.GetErrorString();
			return false;
		}
		OpenMesh::IO::BaseImporter::VHandles vhandles;
		std::vector<OpenMesh::Vec3f> normals;
		std::vector<OpenMesh::Vec2f> texcoords2d;
		std::vector<OpenMesh::Vec3f> texcoords3d;
		// combine all meshes into one
		// TODO assign face group properties when supported
		for (unsigned i = 0; i < scene->mNumMeshes; i++) {
			const aiMesh *mesh = scene->mMeshes[i];
			if (mesh->HasNormals()) _bi.request_hattribs(OpenMesh::AttributeBits::Normal);
			// TODO maybe load other texcoord sets as custom attribs?
			if (mesh->HasTextureCoords(0) && mesh->mNumUVComponents[0] == 2) _bi.request_hattribs(OpenMesh::AttributeBits::TexCoord2D);
			if (mesh->HasTextureCoords(0) && mesh->mNumUVComponents[0] == 3) _bi.request_hattribs(OpenMesh::AttributeBits::TexCoord3D);
			// TODO maybe load other vertex color sets as custom attribs?
			if (mesh->HasVertexColors(0)) _bi.request_hattribs(OpenMesh::AttributeBits::Color);
			if (mesh->HasTangentsAndBitangents()) {
				// openmesh doesnt care
			}
			// add vertices for this mesh
			// TODO individual mesh transforms
			OpenMesh::VertexHandle vh0;
			for (unsigned j = 0; j < mesh->mNumVertices; j++) {
				auto vh = _bi.add_vertex(ai2om(mesh->mVertices[j]));
				if (j == 0) vh0 = vh;
				if (mesh->HasVertexColors(0)) {
					_bi.set_color(vh, ai2om(mesh->mColors[0][j]));
				}
			}
			// add faces for this mesh and their halfedge properties
			for (unsigned j = 0; j < mesh->mNumFaces; j++) {
				vhandles.clear();
				normals.clear();
				texcoords2d.clear();
				texcoords3d.clear();
				const aiFace &face = mesh->mFaces[j];
				for (unsigned k = 0; k < face.mNumIndices; k++) {
					const unsigned uvi = face.mIndices[k];
					// TODO proper overflow check
					vhandles.push_back(OpenMesh::VertexHandle{vh0.idx() + int(uvi)});
					if (mesh->HasNormals()) {
						normals.push_back(ai2om(mesh->mNormals[uvi]));
					}
					if (mesh->HasTextureCoords(0) && mesh->mNumUVComponents[0] == 2) {
						texcoords2d.push_back(OpenMesh::vector_cast<OpenMesh::Vec2f>(ai2om(mesh->mTextureCoords[0][uvi])));
					}
					if (mesh->HasTextureCoords(0) && mesh->mNumUVComponents[0] == 3) {
						texcoords3d.push_back(ai2om(mesh->mTextureCoords[0][uvi]));
					}
				}
				OpenMesh::FaceHandle fh = _bi.add_face(vhandles);
				if (fh.is_valid() && vhandles.size()) {
					_bi.add_face_normals(fh, vhandles[0], normals);
					_bi.add_face_texcoords(fh, vhandles[0], texcoords2d);
					_bi.add_face_texcoords(fh, vhandles[0], texcoords3d);
				}
			}
		}
		return true;
	}

	bool AssimpReader::read(std::istream &_is, OpenMesh::IO::BaseImporter &_bi) {
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
