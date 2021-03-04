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


		bool _OBJWriter_::write(const std::filesystem::path &_filename, BaseExporter &_be) const
		{
			std::fstream out(_filename.c_str(), std::ios_base::out);

			if (!out)
			{
				OMLOG_ERROR << "cannot open file " << _filename.u8string();
				return false;
			}

			// Set precision on output stream. The default is set via IOManager and passed through to all writers.
			out.precision(9);

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

			bool result = write(out, _be);

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

		bool _OBJWriter_::writeMaterial(std::ostream &_out, BaseExporter &_be) const
		{
			OpenMesh::Vec3f c;
			OpenMesh::Vec4f cA;

			material_.clear();
			materialA_.clear();

			//iterate over faces
			for (size_t i = 0, nF = _be.n_faces(); i < nF; ++i)
			{
				//color with alpha
				if (_be.file_options().color_has_alpha()) {
					cA = color_cast<OpenMesh::Vec4f> (_be.colorA(FaceHandle(int(i))));
					getMaterial(cA);
				} else {
					//and without alpha
					c = color_cast<OpenMesh::Vec3f> (_be.color(FaceHandle(int(i))));
					getMaterial(c);
				}
			}

			//write the materials
			if (_be.file_options().color_has_alpha())
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


		bool _OBJWriter_::write(std::ostream &_out, BaseExporter &_be) const
		{
			// TODO face color to material needs a re-write

			if (!_out.good()) {
				OMLOG_ERROR << "cannot write to stream";
				return false;
			}

			OMLOG_INFO << "write file";

			_out.precision(9);

			// No binary mode for OBJ
			if (_be.file_options().check(OptionBits::Binary)) {
				OMLOG_WARNING << "Binary mode not supported, falling back to ascii";
			}

			// check for unsupported writer features
			if (_be.file_options().face_has_normal()) {
				OMLOG_WARNING << "Face normals not supported";
			}

			// cant export halfedge and vertex normals
			if (_be.file_options().vertex_has_normal() && _be.file_options().halfedge_has_normal()) {
				OMLOG_WARNING << "Prioritizing halfedge normals over vertex normals";
			}

			// cant export halfedge and vertex colors
			// not supporting proper halfedge color export (unlikely to be needed)
			if (_be.file_options().halfedge_has_color()) {
				if (_be.file_options().vertex_has_color()) {
					OMLOG_WARNING << "Prioritizing vertex colors over halfedge colors";
				} else {
					OMLOG_WARNING << "Truncating halfedge colors to vertex colors";
				}
			}

			bool do_texcoord2D = _be.file_options().halfedge_has_texcoord2D() || _be.file_options().vertex_has_texcoord2D();
			bool do_texcoord3D = _be.file_options().halfedge_has_texcoord3D() || _be.file_options().vertex_has_texcoord3D();
			bool do_normal = _be.file_options().halfedge_has_normal() || _be.file_options().vertex_has_normal();
			bool do_color = _be.file_options().halfedge_has_color() || _be.file_options().vertex_has_color();

			if (do_texcoord2D && do_texcoord3D) {
				OMLOG_WARNING << "Cannot export both 2D and 3D texture coords, prioritizing 2D";
				do_texcoord3D = false;
			}

			if (do_color && _be.file_options().color_has_alpha()) {
				OMLOG_WARNING << "Alpha color component not supported";
			}

			// header
			_out << "# " << _be.n_vertices() << " vertices, ";
			_out << _be.n_faces() << " faces" << '\n';

			// material file
			// _out << "mtllib " << objName_ << ".mtl" << '\n';

			// compressed halfedge properties
			std::vector<int> prop_tc_idx, prop_n_idx;

			if (do_texcoord2D && _be.file_options().halfedge_has_texcoord2D()) {
				// compress and write 2D texcoords
				OMLOG_DEBUG << "compressing halfedge texcoords2D";
				_be.compress_halfedge_properties<Vec2f>(
					prop_tc_idx, 0,
					[](const BaseExporter &be, VertexHandle vh, HalfedgeHandle hh) {
						return be.texcoord2D(hh);
					},
					[&](int, Vec2f &&tc) {
						_out << "vt " << tc[0] << " " << tc[1] << '\n';
					}
				);
			} else if (do_texcoord2D) {
				// write vertex 2D texcoords
				for (int i = 0; i < _be.n_vertices(); ++i) {
					Vec2f tc = _be.texcoord2D(VertexHandle(i));
					_out << "vt " << tc[0] << " " << tc[1] << '\n';
				}
			}

			if (do_texcoord3D && _be.file_options().halfedge_has_texcoord3D()) {
				// compress and write 3D texcoords
				OMLOG_DEBUG << "compressing halfedge texcoords3D";
				_be.compress_halfedge_properties<Vec3f>(
					prop_tc_idx, 0,
					[](const BaseExporter &be, VertexHandle vh, HalfedgeHandle hh) {
						return be.texcoord3D(hh);
					},
					[&](int, Vec3f &&tc) {
						_out << "vt " << tc[0] << " " << tc[1] << " " << tc[2] << '\n';
					}
				);
			} else if (do_texcoord3D) {
				// write vertex 3D texcoords
				for (int i = 0; i < _be.n_vertices(); ++i) {
					Vec3f tc = _be.texcoord3D(VertexHandle(i));
					_out << "vt " << tc[0] << " " << tc[1] << " " << tc[2] << '\n';
				}
			}

			if (_be.file_options().halfedge_has_normal()) {
				// compress and write normals
				OMLOG_DEBUG << "compressing halfedge normals";
				_be.compress_halfedge_properties<Vec3f>(
					prop_n_idx, 0,
					[](const BaseExporter &be, VertexHandle vh, HalfedgeHandle hh) {
						return be.normal(hh);
					},
					[&](int, Vec3f &&n) {
						_out << "vn " << n[0] << " " << n[1] << " " << n[2] << '\n';
					}
				);
			} else if (_be.file_options().vertex_has_normal()) {
				// write vertex normals
				for (int i = 0; i < _be.n_vertices(); ++i) {
					Vec3f n = _be.normal(VertexHandle(i));
					_out << "vn " << n[0] << " " << n[1] << " " << n[2] << '\n';
				}
			}

			// vertex points
			for (int i = 0; i < _be.n_vertices(); ++i) {
				VertexHandle vh{i};
				Vec3f p = _be.point(vh);
				_out << "v " << p[0] << " " << p[1] << " " << p[2];
				if (_be.file_options().vertex_has_color() || _be.file_options().halfedge_has_color()) {
					// alpha not supported, always float
					auto cf = _be.colorf(vh);
					_out << " " << cf[0] << " " << cf[1] << " " << cf[2];
				}
				_out << '\n';
			}

			// we do not want to write seperators if we only write vertex indices
			bool onlyVertices = !do_normal && !do_texcoord2D && !do_texcoord3D;

			std::vector<VertexHandle> vhandles;
			std::vector<HalfedgeHandle> hhandles;

			// faces (indices starting at 1 not 0)
			for (int i = 0; i < _be.n_faces(); ++i) {

				_out << "f";

				_be.face_vertex_handles(FaceHandle(int(i)), vhandles);
				_be.face_halfedge_handles(FaceHandle(int(i)), hhandles);

				for (int j = 0; j < vhandles.size(); ++j) {

					VertexHandle vh = vhandles[j];
					HalfedgeHandle hh = hhandles[j];

					// Write vertex index
					_out << " " << (vh.idx() + 1);

					if (!onlyVertices) {
						_out << "/";

						if (_be.file_options().halfedge_has_texcoord()) {
							// write texcoord index from halfedge
							_out << (prop_tc_idx[hh.idx()] + 1);
						} else if (_be.file_options().vertex_has_texcoord()) {
							// write texcoord index from vertex
							_out << (vh.idx() + 1);
						}

						if (_be.file_options().halfedge_has_normal()) {
							// write normal index from halfedge
							_out << "/";
							_out << (prop_n_idx[hh.idx()] + 1);
						} else if (_be.file_options().vertex_has_normal()) {
							// write normal index from vertex
							_out << "/";
							_out << (vh.idx() + 1);
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
