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


 // OpenMesh
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
// STL
#if defined(OM_CC_MIPS)
#  include <ctype.h>
/// \bug Workaround for STLPORT 4.6: isspace seems not to be in namespace std!
#elif defined(_STLPORT_VERSION) && (_STLPORT_VERSION==0x460)
#  include <cctype>
#else
using std::isspace;
#endif

#ifndef WIN32
#endif

#include <fstream>

#define OMLOG_SOURCE OBJReader

//=== NAMESPACES ==============================================================


namespace OpenMesh {
	namespace IO {


		//=== INSTANCIATE =============================================================


		_OBJReader_  __OBJReaderInstance;
		_OBJReader_ &OBJReader() { return __OBJReaderInstance; }


		//=== IMPLEMENTATION ==========================================================

		//-----------------------------------------------------------------------------

		void trimString(std::string &_string) {
			// Trim Both leading and trailing spaces

			size_t start = _string.find_first_not_of(" \t\r\n");
			size_t end = _string.find_last_not_of(" \t\r\n");

			if ((std::string::npos == start) || (std::string::npos == end))
				_string = "";
			else
				_string = _string.substr(start, end - start + 1);
		}

		//-----------------------------------------------------------------------------

		// remove duplicated indices from one face
		void remove_duplicated_vertices(BaseImporter::VHandles &_indices)
		{
			BaseImporter::VHandles::iterator endIter = _indices.end();
			for (BaseImporter::VHandles::iterator iter = _indices.begin(); iter != endIter; ++iter)
				endIter = std::remove(iter + 1, endIter, *(iter));

			_indices.erase(endIter, _indices.end());
		}

		//-----------------------------------------------------------------------------

		_OBJReader_::
			_OBJReader_()
		{
			IOManager().register_module(this);
		}


		//-----------------------------------------------------------------------------


		bool _OBJReader_::read(const std::filesystem::path &_filename, BaseImporter &_bi)
		{
			std::fstream in(_filename, std::ios_base::in);

			if (!in.is_open() || !in.good())
			{
				OMLOG_ERROR << "cannot not open file " << _filename.u8string();
				return false;
			}

			{
#if defined(WIN32)
				//std::string::size_type dot = _filename.find_last_of("\\/");
#else
				//std::string::size_type dot = _filename.rfind("/");
#endif
				path_ = _filename.parent_path();
				//path_ = (dot == std::string::npos)
				//  ? "./"
				//  : std::string(_filename.substr(0,dot+1));
			}

			bool result = read(in, _bi);

			in.close();
			return result;
		}

		//-----------------------------------------------------------------------------

		bool
			_OBJReader_::
			read_material(std::fstream &_in)
		{
			std::string line;
			std::string keyWrd;
			std::string textureName;

			std::stringstream  stream;

			std::string key;
			Material    mat;
			float       f1, f2, f3;
			bool        indef = false;
			int         textureId = 1;


			materials_.clear();
			mat.cleanup();

			while (_in && !_in.eof())
			{
				std::getline(_in, line);
				if (_in.bad()) {
					OMLOG_WARNING << "Could not read material file properly!";
					return false;
				}

				if (line.empty())
					continue;

				stream.str(line);
				stream.clear();

				stream >> keyWrd;

				if ((isspace(line[0]) && line[0] != '\t') || line[0] == '#')
				{
					if (indef && !key.empty() && mat.is_valid())
					{
						materials_[key] = mat;
						mat.cleanup();
					}
				}

				else if (keyWrd == "newmtl") // begin new material definition
				{
					stream >> key;
					indef = true;
				}

				else if (keyWrd == "Kd") // diffuse color
				{
					stream >> f1; stream >> f2; stream >> f3;

					if (!stream.fail())
						mat.set_Kd(f1, f2, f3);
				}

				else if (keyWrd == "Ka") // ambient color
				{
					stream >> f1; stream >> f2; stream >> f3;

					if (!stream.fail())
						mat.set_Ka(f1, f2, f3);
				}

				else if (keyWrd == "Ks") // specular color
				{
					stream >> f1; stream >> f2; stream >> f3;

					if (!stream.fail())
						mat.set_Ks(f1, f2, f3);
				}
#if 0
				else if (keyWrd == "illum") // diffuse/specular shading model
				{
					; // just skip this
				}

				else if (keyWrd == "Ns") // Shininess [0..200]
				{
					; // just skip this
				}

				else if (keyWrd == "map_") // map images
				{
					// map_Ks, specular map
					// map_Ka, ambient map
					// map_Bump, bump map
					// map_d,  opacity map
					; // just skip this
				}
#endif
				else if (keyWrd == "map_Kd") {
					// Get the rest of the line, removing leading or trailing spaces
					// This will define the filename of the texture
					std::getline(stream, textureName);
					trimString(textureName);
					if (! textureName.empty())
						mat.set_map_Kd(textureName, textureId++);
				} else if (keyWrd == "Tr") // transparency value
				{
					stream >> f1;

					if (!stream.fail())
						mat.set_Tr(f1);
				} else if (keyWrd == "d") // transparency value
				{
					stream >> f1;

					if (!stream.fail())
						mat.set_Tr(f1);
				}

				if (_in && indef && mat.is_valid() && !key.empty())
					materials_[key] = mat;
			}
			return true;
		}
		//-----------------------------------------------------------------------------

		bool
			_OBJReader_::
			read_vertices(std::istream &_in, BaseImporter &_bi,
				std::vector<Vec3f> &normals,
				std::vector<Vec3f> &colors,
				std::vector<Vec3f> &texcoords3d,
				std::vector<Vec2f> &texcoords,
				std::vector<VertexHandle> &vertexHandles
			)
		{
			float x, y, z, u, v, w;
			float r, g, b;

			std::string line;
			std::string keyWrd;

			std::stringstream stream;

			while (_in && !_in.eof())
			{
				std::getline(_in, line);
				if (_in.bad()) {
					OMLOG_WARNING << "Could not read file properly!";
					return false;
				}

				// Trim Both leading and trailing spaces
				trimString(line);

				// comment
				if (line.size() == 0 || line[0] == '#' || isspace(line[0])) {
					continue;
				}

				stream.str(line);
				stream.clear();

				stream >> keyWrd;

				// vertex
				if (keyWrd == "v")
				{
					stream >> x; stream >> y; stream >> z;

					if (!stream.fail())
					{
						vertexHandles.push_back(_bi.add_vertex(OpenMesh::Vec3f(x, y, z)));
						stream >> r; stream >> g; stream >> b;

						if (!stream.fail())
						{
							if (!!_bi.request_vattribs(AttributeBits::Color)) {
								colors.push_back(OpenMesh::Vec3f(r, g, b));
							}
						}
					}
				}

				// texture coord
				else if (keyWrd == "vt")
				{
					stream >> u; stream >> v;

					if (!stream.fail()) {

						// try to read the w component as it is optional
						stream >> w;

						if (stream.fail()) {
							// 2d
							if (!!_bi.request_hattribs(AttributeBits::TexCoord2D)) {
								texcoords.push_back(OpenMesh::Vec2f(u, v));
							}
						} else {
							// 3d
							if (!!_bi.request_hattribs(AttributeBits::TexCoord3D)) {
								texcoords3d.push_back(OpenMesh::Vec3f(u, v, w));
							}
						}

					} else {
						OMLOG_ERROR << "Only single 2D or 3D texture coordinate per vertex allowed!";
						return false;
					}
				}

				// color per vertex
				else if (keyWrd == "vc")
				{
					stream >> r; stream >> g; stream >> b;

					if (!stream.fail()) {
						if (!!_bi.request_vattribs(AttributeBits::Color)) {
							colors.push_back(OpenMesh::Vec3f(r, g, b));
						}
					}
				}

				// normal
				else if (keyWrd == "vn")
				{
					stream >> x; stream >> y; stream >> z;

					if (!stream.fail()) {
						if (!!_bi.request_hattribs(AttributeBits::Normal)) {
							normals.push_back(OpenMesh::Vec3f(x, y, z));
						}
					}
				}
			}

			if (texcoords.size() && texcoords3d.size()) {
				OMLOG_ERROR << "file has inconsistent 2D/3D texcoords, they will be broken";
			}

			return true;
		}

		//-----------------------------------------------------------------------------

		bool _OBJReader_::read(std::istream &_in, BaseImporter &_bi)
		{
			std::string line;
			std::string keyWrd;

			std::vector<Vec3f>        normals;
			std::vector<Vec3f>        colors;
			std::vector<Vec3f>        texcoords3d;
			std::vector<Vec2f>        texcoords;
			std::vector<VertexHandle> vertexHandles;

			BaseImporter::VHandles    vhandles;
			std::vector<Vec3f>        face_texcoords3d;
			std::vector<Vec2f>        face_texcoords;
			std::vector<Vec3f>        face_normals;

			std::string               matname;

			std::stringstream         stream, lineData, tmp;

			// pass 1: read vertices
			if (!read_vertices(_in, _bi,
				normals, colors, texcoords3d, texcoords,
				vertexHandles)) {
				return false;
			}

			// reset stream for second pass
			_in.clear();
			_in.seekg(0, std::ios::beg);

			int nCurrentPositions = 0,
				nCurrentTexcoords = 0,
				nCurrentNormals = 0;

			// pass 2: read faces
			while (_in && !_in.eof())
			{
				std::getline(_in, line);
				if (_in.bad()) {
					OMLOG_WARNING << "Could not read file properly!";
					return false;
				}

				// Trim Both leading and trailing spaces
				trimString(line);

				// comment
				if (line.size() == 0 || line[0] == '#' || isspace(line[0])) {
					continue;
				}

				stream.str(line);
				stream.clear();

				stream >> keyWrd;

				// material file
				if (keyWrd == "mtllib")
				{
					std::string matFile;

					// Get the rest of the line, removing leading or trailing spaces
					// This will define the filename of the texture
					std::getline(stream, matFile);
					trimString(matFile);

					auto matPath = path_ / std::filesystem::u8path(matFile);

					//omlog() << "Load material file " << matFile << std::endl;

					std::fstream matStream(matPath, std::ios_base::in);

					if (matStream) {

						if (!read_material(matStream)) matStream.close();

					} else
						OMLOG_WARNING << "Material file '" << matFile << "' not found!";

					//omlog() << "  " << materials_.size() << " materials loaded.\n";

					for (MaterialList::iterator material = materials_.begin(); material != materials_.end(); ++material)
					{
						// Save the texture information in a property
						if ((*material).second.has_map_Kd())
							_bi.add_texture_information((*material).second.map_Kd_index(), (*material).second.map_Kd());
					}

				}

				// usemtl
				else if (keyWrd == "usemtl")
				{
					stream >> matname;
					if (materials_.find(matname) == materials_.end())
					{
						OMLOG_WARNING << "Material '" << matname << "' not defined in material file.";
						matname = "";
					}
				}

				// track current number of parsed vertex attributes,
				// to allow for OBJs negative indices
				else if (keyWrd == "v")
				{
					++nCurrentPositions;
				} else if (keyWrd == "vt")
				{
					++nCurrentTexcoords;
				} else if (keyWrd == "vn")
				{
					++nCurrentNormals;
				}

				// faces
				else if (keyWrd == "f")
				{
					int component(0), nV(0);
					int value;

					vhandles.clear();
					face_texcoords.clear();
					face_texcoords3d.clear();
					face_normals.clear();

					// read full line after detecting a face
					std::string faceLine;
					std::getline(stream, faceLine);
					lineData.str(faceLine);
					lineData.clear();

					FaceHandle fh;
					BaseImporter::VHandles faceVertices;

					// work on the line until nothing left to read
					while (!lineData.eof())
					{
						// read one block from the line ( vertex/texCoord/normal )
						std::string vertex;
						lineData >> vertex;

						do {

							//get the component (vertex/texCoord/normal)
							size_t found = vertex.find("/");

							// parts are seperated by '/' So if no '/' found its the last component
							if (found != std::string::npos) {

								// read the index value
								tmp.str(vertex.substr(0, found));
								tmp.clear();

								// If we get an empty string this property is undefined in the file
								if (vertex.substr(0, found).empty()) {
									// Switch to next field
									vertex = vertex.substr(found + 1);

									// Now we are at the next component
									++component;

									// Skip further processing of this component
									continue;
								}

								// Read current value
								tmp >> value;

								// remove the read part from the string
								vertex = vertex.substr(found + 1);

							} else {

								// last component of the vertex, read it.
								tmp.str(vertex);
								tmp.clear();
								tmp >> value;

								// Clear vertex after finished reading the line
								vertex = "";

								// Nothing to read here ( garbage at end of line )
								if (tmp.fail()) {
									continue;
								}
							}

							// store the component ( each component is referenced by the index here! )
							switch (component)
							{
							case 0: // vertex
								if (value < 0) {
									// Calculation of index :
									// -1 is the last vertex in the list
									// As obj counts from 1 and not zero add +1
									value = nCurrentPositions + value + 1;
								}
								// Obj counts from 1 and not zero .. array counts from zero therefore -1
								vhandles.push_back(VertexHandle(value - 1));
								faceVertices.push_back(VertexHandle(value - 1));
								if (!!_bi.want_vattribs(AttributeBits::Color)) {
									if ((unsigned int) (value - 1) < colors.size()) {
										_bi.set_color(vhandles.back(), colors[value - 1]);
									} else {
										OMLOG_ERROR << "Error setting vertex color at " << value;
									}
								}
								break;

							case 1: // texture coord
								if (value < 0) {
									// Calculation of index :
									// -1 is the last vertex in the list
									// As obj counts from 1 and not zero add +1
									value = nCurrentTexcoords + value + 1;
								}
								assert(!vhandles.empty());

								// Obj counts from 1 and not zero .. array counts from zero therefore -1
								if (!!_bi.want_hattribs(AttributeBits::TexCoord2D)) {
									if (!texcoords.empty() && (unsigned int) (value - 1) < texcoords.size()) {
										face_texcoords.push_back(texcoords[value - 1]);
									} else {
										OMLOG_ERROR << "Error setting texture coordinates at " << value;
									}
								}

								if (!!_bi.want_hattribs(AttributeBits::TexCoord3D)) {
									if (!texcoords3d.empty() && (unsigned int) (value - 1) < texcoords3d.size()) {
										face_texcoords3d.push_back(texcoords3d[value - 1]);
									} else {
										OMLOG_ERROR << "Error setting texture coordinates at " << value;
									}
								}

								break;

							case 2: // normal
								if (value < 0) {
									// Calculation of index :
									// -1 is the last vertex in the list
									// As obj counts from 1 and not zero add +1
									value = nCurrentNormals + value + 1;
								}

								// Obj counts from 1 and not zero .. array counts from zero therefore -1
								if (!!_bi.want_hattribs(AttributeBits::Normal)) {
									if (!normals.empty() && (unsigned int) (value - 1) < normals.size()) {
										face_normals.push_back(normals[value - 1]);
									} else {
										OMLOG_ERROR << "Error setting vertex normal at " << value;
									}
								}

								break;
							}

							// Prepare for reading next component
							++component;

							// Read until line does not contain any other info
						} while (!vertex.empty());

						component = 0;
						nV++;

					}

					// note that add_face can possibly triangulate the faces, which is why we have to
					// store the current number of faces first
					size_t n_faces = _bi.n_faces();
					remove_duplicated_vertices(faceVertices);

					//A minimum of three vertices are required.
					if (faceVertices.size() > 2)
						fh = _bi.add_face(faceVertices);

					if (!vhandles.empty() && fh.is_valid())
					{
						_bi.add_face_texcoords(fh, vhandles[0], face_texcoords);
						_bi.add_face_texcoords(fh, vhandles[0], face_texcoords3d);
						_bi.add_face_normals(fh, vhandles[0], face_normals);
					}

					std::vector<FaceHandle> newfaces;
					for (size_t i = 0; i < _bi.n_faces() - n_faces; ++i) {
						newfaces.push_back(FaceHandle(int(n_faces + i)));
					}

					if (!matname.empty())
					{
						Material &mat = materials_[matname];

						if (mat.has_Kd()) {
							// apply mtl color as face color
							Vec3uc fc = color_cast<Vec3uc, Vec3f>(mat.Kd());
							if (!!_bi.request_fattribs(AttributeBits::Color)) {
								for (auto it = newfaces.begin(); it != newfaces.end(); ++it) {
									_bi.set_color(*it, fc);
								}
							}
						}

						if (mat.has_map_Kd()) {
							// Set the texture index in the face index property
							if (!!_bi.request_fattribs(AttributeBits::TextureIndex)) {
								for (auto it = newfaces.begin(); it != newfaces.end(); ++it) {
									_bi.set_face_texindex(*it, mat.map_Kd_index());
								}
							}
						} else {
							// If we don't have the info, set it to no texture
							if (!!_bi.request_fattribs(AttributeBits::TextureIndex)) {
								for (auto it = newfaces.begin(); it != newfaces.end(); ++it) {
									_bi.set_face_texindex(*it, 0);
								}
							}
						}

					} else {
						// Set the texture index to zero as we don't have any information
						if (!!_bi.request_fattribs(AttributeBits::TextureIndex)) {
							for (auto it = newfaces.begin(); it != newfaces.end(); ++it) {
								_bi.set_face_texindex(*it, 0);
							}
						}
					}

				}

			}

			// If we do not have any faces,
			// assume this is a point cloud and read the normals and colors directly
			if (_bi.n_faces() == 0)
			{
				int i = 0;

				// add normal per vertex
				if (normals.size() == _bi.n_vertices()) {
					if (!!_bi.request_vattribs(AttributeBits::Normal)) {
						for (auto it = vertexHandles.begin(); it != vertexHandles.end(); ++it, i++) {
							_bi.set_normal(*it, normals[i]);
						}
					}
				}

				// add color per vertex
				i = 0;
				if (colors.size() >= _bi.n_vertices()) {
					if (!!_bi.request_vattribs(AttributeBits::Color)) {
						for (auto it = vertexHandles.begin(); it != vertexHandles.end(); ++it, i++)
							_bi.set_color(*it, colors[i]);
					}
				}

			}

			return true;
		}


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
