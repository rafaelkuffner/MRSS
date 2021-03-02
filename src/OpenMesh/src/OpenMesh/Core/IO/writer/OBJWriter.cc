/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */




 //== INCLUDES =================================================================


 //STL
#include <fstream>
#include <limits>
#include <map>

// OpenMesh
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>

#define OMLOG_SOURCE OBJWriter

//=== NAMESPACES ==============================================================


namespace OpenMesh {
	namespace IO {


		//=== INSTANCIATE =============================================================


		// register the OBJLoader singleton with MeshLoader
		_OBJWriter_  __OBJWriterinstance;
		_OBJWriter_ &OBJWriter() { return __OBJWriterinstance; }


		//=== IMPLEMENTATION ==========================================================


		_OBJWriter_::_OBJWriter_() { IOManager().register_module(this); }


		//-----------------------------------------------------------------------------


		bool
			_OBJWriter_::
			write(const std::filesystem::path &_filename, BaseExporter &_be, Options _opt, std::streamsize _precision) const
		{
			std::fstream out(_filename.c_str(), std::ios_base::out);

			if (!out)
			{
				OMLOG_ERROR << "cannot open file " << _filename.u8string();
				return false;
			}

			// Set precision on output stream. The default is set via IOManager and passed through to all writers.
			out.precision(_precision);

			// Set fixed output to avoid problems with programs not reading scientific notation correctly
			out << std::fixed;

			{
#if defined(WIN32)
				//std::string::size_type dot = _filename.find_last_of("\\/");
#else
				//std::string::size_type dot = _filename.rfind("/");
#endif

				path_ = _filename.parent_path();
				objName_ = _filename.filename();

				//if (dot == std::string::npos){
				//  path_ = "./";
				//  objName_ = _filename;
				//}else{
				//  path_ = _filename.substr(0,dot+1);
				//  objName_ = _filename.substr(dot+1);
				//}

				//remove the file extension
				//dot = objName_.find_last_of(".");

				//if(dot != std::string::npos)
				  //objName_ = objName_.substr(0,dot);
				objName_.replace_extension("");
			}

			bool result = write(out, _be, _opt, _precision);

			out.close();
			return result;
		}

		//-----------------------------------------------------------------------------

		size_t _OBJWriter_::getMaterial(OpenMesh::Vec3f _color) const
		{
			for (size_t i = 0; i < material_.size(); i++)
				if (material_[i] == _color)
					return i;

			//not found add new material
			material_.push_back(_color);
			return material_.size() - 1;
		}

		//-----------------------------------------------------------------------------

		size_t _OBJWriter_::getMaterial(OpenMesh::Vec4f _color) const
		{
			for (size_t i = 0; i < materialA_.size(); i++)
				if (materialA_[i] == _color)
					return i;

			//not found add new material
			materialA_.push_back(_color);
			return materialA_.size() - 1;
		}

		//-----------------------------------------------------------------------------

		bool
			_OBJWriter_::
			writeMaterial(std::ostream &_out, BaseExporter &_be, Options _opt) const
		{
			OpenMesh::Vec3f c;
			OpenMesh::Vec4f cA;

			material_.clear();
			materialA_.clear();

			//iterate over faces
			for (size_t i = 0, nF = _be.n_faces(); i < nF; ++i)
			{
				//color with alpha
				if (_opt.color_has_alpha()) {
					cA = color_cast<OpenMesh::Vec4f> (_be.colorA(FaceHandle(int(i))));
					getMaterial(cA);
				} else {
					//and without alpha
					c = color_cast<OpenMesh::Vec3f> (_be.color(FaceHandle(int(i))));
					getMaterial(c);
				}
			}

			//write the materials
			if (_opt.color_has_alpha())
				for (size_t i = 0; i < materialA_.size(); i++) {
					_out << "newmtl " << "mat" << i << '\n';
					_out << "Ka 0.5000 0.5000 0.5000" << '\n';
					_out << "Kd " << materialA_[i][0] << ' ' << materialA_[i][1] << ' ' << materialA_[i][2] << '\n';
					_out << "Tr " << materialA_[i][3] << '\n';
					_out << "illum 1" << '\n';
				} else
					for (size_t i = 0; i < material_.size(); i++) {
						_out << "newmtl " << "mat" << i << '\n';
						_out << "Ka 0.5000 0.5000 0.5000" << '\n';
						_out << "Kd " << material_[i][0] << ' ' << material_[i][1] << ' ' << material_[i][2] << '\n';
						_out << "illum 1" << '\n';
					}

				return true;
		}

		//-----------------------------------------------------------------------------


		bool
			_OBJWriter_::
			write(std::ostream &_out, BaseExporter &_be, Options _opt, std::streamsize _precision) const
		{
			unsigned int idx;
			size_t i, j, nV, nF;
			Vec3f v, n;
			Vec2f t;
			VertexHandle vh;
			std::vector<VertexHandle> vhandles;
			bool useMatrial = false;
			OpenMesh::Vec3f c;
			OpenMesh::Vec4f cA;

			OMLOG_INFO << "write file";

			_out.precision(_precision);

			// No binary mode for OBJ
			if (_opt.check(OptionBits::Binary)) {
				OMLOG_WARNING << "Binary mode not supported by OBJ Writer, falling back to standard";
			}

			// check for unsupported writer features
			if (_opt.face_has_normal()) {
				OMLOG_WARNING << "FaceNormal not supported by OBJ Writer";
				_opt.fattribs &= ~AttributeBits::Normal;
			}

			// check for unsupported writer features
			if (_opt.vertex_has_color()) {
				// TODO implement OBJ vertex color export
				OMLOG_WARNING << "VertexColor not supported by OBJ Writer";
				_opt.vattribs &= ~AttributeBits::Color;
			}

			//create material file if needed
			if (_opt.face_has_color()) {

				auto matFile = path_ / objName_;
				matFile.replace_extension(".mat");
				// shouldnt it be .mtl ? - ben

				std::fstream matStream(matFile, std::ios_base::out);

				if (!matStream)
				{
					OMLOG_ERROR << "[OBJWriter] : cannot write material file " << matFile.u8string();

				} else {
					useMatrial = writeMaterial(matStream, _be, _opt);

					matStream.close();
				}
			}

			// header
			_out << "# " << _be.n_vertices() << " vertices, ";
			_out << _be.n_faces() << " faces" << '\n';

			// material file
			if (useMatrial && _opt.face_has_color())
				_out << "mtllib " << objName_ << ".mat" << '\n';

			// TODO export 3d texcoords?
			std::map<Vec2f, int> texMap;
			//collect Texturevertices from halfedges
			if (_opt.halfedge_has_texcoord2D())
			{
				std::vector<Vec2f> texCoords;
				//add all texCoords to map
				unsigned int num = _be.get_face_texcoords(texCoords);
				for (unsigned int i = 0; i < num; ++i)
				{
					texMap[texCoords[i]] = i;
				}
			}

			//collect Texture coordinates from vertices
			if (_opt.vertex_has_texcoord2D())
			{
				for (size_t i = 0, nV = _be.n_vertices(); i < nV; ++i)
				{
					vh = VertexHandle(static_cast<int>(i));
					t = _be.texcoord(vh);
					texMap[t] = static_cast<int>(i);
				}
			}

			// assign each texcoord in the map its id
			// and write the vt entries
			if (_opt.vertex_has_texcoord2D() || _opt.halfedge_has_texcoord2D())
			{
				int texCount = 0;
				for (std::map<Vec2f, int>::iterator it = texMap.begin(); it != texMap.end(); ++it)
				{
					_out << "vt " << it->first[0] << " " << it->first[1] << '\n';
					it->second = ++texCount;
				}
			}

			// vertex data (point, normals, texcoords)
			for (i = 0, nV = _be.n_vertices(); i < nV; ++i)
			{
				vh = VertexHandle(int(i));
				v = _be.point(vh);
				n = _be.normal(vh);
				t = _be.texcoord(vh);

				_out << "v " << v[0] << " " << v[1] << " " << v[2] << '\n';

				if (_opt.vertex_has_normal())
					_out << "vn " << n[0] << " " << n[1] << " " << n[2] << '\n';
			}

			size_t lastMat = std::numeric_limits<std::size_t>::max();

			// we do not want to write seperators if we only write vertex indices
			bool onlyVertices = !_opt.vertex_has_normal()
				&& !_opt.vertex_has_texcoord()
				&& !_opt.halfedge_has_normal()
				&& !_opt.halfedge_has_texcoord();

			// faces (indices starting at 1 not 0)
			for (i = 0, nF = _be.n_faces(); i < nF; ++i)
			{

				if (useMatrial && _opt.face_has_color()) {
					size_t material = std::numeric_limits<std::size_t>::max();

					//color with alpha
					if (_opt.color_has_alpha()) {
						cA = color_cast<OpenMesh::Vec4f> (_be.colorA(FaceHandle(int(i))));
						material = getMaterial(cA);
					} else {
						//and without alpha
						c = color_cast<OpenMesh::Vec3f> (_be.color(FaceHandle(int(i))));
						material = getMaterial(c);
					}

					// if we are ina a new material block, specify in the file which material to use
					if (lastMat != material) {
						_out << "usemtl mat" << material << '\n';
						lastMat = material;
					}
				}

				_out << "f";

				_be.get_vhandles(FaceHandle(int(i)), vhandles);

				for (j = 0; j < vhandles.size(); ++j)
				{

					// Write vertex index
					idx = vhandles[j].idx() + 1;
					_out << " " << idx;

					if (!onlyVertices) {
						// write separator
						_out << "/";

						//write texCoords index from halfedge
						if (_opt.halfedge_has_texcoord())
						{
							_out << texMap[_be.texcoord(_be.getHeh(FaceHandle(int(i)), vhandles[j]))];
						}

						else
						{
							// write vertex texture coordinate index
							if (_opt.vertex_has_texcoord())
								_out << texMap[_be.texcoord(vhandles[j])];
						}

						// FIXME export halfedge normals

						// write vertex normal index
						if (_opt.vertex_has_normal()) {
							// write separator
							_out << "/";
							_out << idx;
						}
					}
				}

				_out << '\n';
			}

			material_.clear();
			materialA_.clear();

			return true;
		}


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
