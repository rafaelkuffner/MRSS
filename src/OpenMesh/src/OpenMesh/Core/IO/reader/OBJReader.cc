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

#include <iostream>

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

#include <iobuffer/iobuffer.hpp>

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

		// WHAT IS THE PURPOSE OF THIS? why did the original OpenMesh parser use it?
		// remove duplicated indices from one face
		void remove_duplicated_vertices(BaseImporter::VHandles &_indices)
		{
			BaseImporter::VHandles::iterator endIter = _indices.end();
			for (BaseImporter::VHandles::iterator iter = _indices.begin(); iter != endIter; ++iter)
				endIter = std::remove(iter + 1, endIter, *(iter));

			_indices.erase(endIter, _indices.end());
		}

		class Material
		{
		public:

			Material() = default;

			Material(std::string name) : name_(std::move(name)) {}

			bool is_valid(void) const
			{
				return name_.size() && (Kd_is_set_ || Ka_is_set_ || Ks_is_set_ || Tr_is_set_ || map_Kd_is_set_);
			}

			const std::string & name() const {
				return name_;
			}

			bool has_Kd(void) { return Kd_is_set_; }
			bool has_Ka(void) { return Ka_is_set_; }
			bool has_Ks(void) { return Ks_is_set_; }
			bool has_Tr(void) { return Tr_is_set_; }
			bool has_map_Kd(void) { return map_Kd_is_set_; }

			void set_Kd(float r, float g, float b)
			{
				Kd_ = Vec3f(r, g, b); Kd_is_set_ = true;
			}

			void set_Ka(float r, float g, float b)
			{
				Ka_ = Vec3f(r, g, b); Ka_is_set_ = true;
			}

			void set_Ks(float r, float g, float b)
			{
				Ks_ = Vec3f(r, g, b); Ks_is_set_ = true;
			}

			void set_Tr(float t)
			{
				Tr_ = t;            Tr_is_set_ = true;
			}

			void set_map_Kd(std::string _name, int _index_Kd)
			{
				map_Kd_ = _name, index_Kd_ = _index_Kd; map_Kd_is_set_ = true;
			};

			const Vec3f &Kd(void) const { return Kd_; }
			const Vec3f &Ka(void) const { return Ka_; }
			const Vec3f &Ks(void) const { return Ks_; }
			float  Tr(void) const { return Tr_; }
			const std::string &map_Kd(void) { return map_Kd_; }
			const int &map_Kd_index(void) { return index_Kd_; }

		private:
			std::string name_;

			Vec3f Kd_; // diffuse
			Vec3f Ka_; // ambient
			Vec3f Ks_; // specular
			float Tr_; // transperency

			int index_Kd_ = 0;
			std::string map_Kd_;  // Texture

			bool Kd_is_set_ = false;
			bool Ka_is_set_ = false;
			bool Ks_is_set_ = false;
			bool Tr_is_set_ = false;
			bool map_Kd_is_set_ = false;

		};

		class MTLParser {
		private:
			iob::text_reader m_in;

		public:
			MTLParser() = default;

			MTLParser(const MTLParser &) = delete;
			MTLParser & operator=(const MTLParser &) = delete;

			MTLParser(iob::iobuffer *buf_) :
				m_in{buf_}
			{}

			std::vector<Material> parse() {
				std::vector<Material> mtls;
				Material mtl;
				// TODO this is for the OpenMesh texture index property and should be per-OBJ not per-MTL
				int textureId = 1;
				while (m_in.good()) {
					// skip leading whitespace (including blank lines)
					m_in.get_while_ws();
					// skip comments
					if (m_in.has_any('#')) {
						m_in.get_line();
						continue;
					}
					std::string keyword{m_in.get_until_ws()};
					if (keyword == "newmtl") {
						if (mtl.is_valid()) mtls.push_back(std::move(mtl));
						m_in.get_while_ws();
						mtl = Material(std::string(m_in.get_until_ws()));
					} else if (keyword == "Kd") {
						float r = 0, g = 0, b = 0;
						auto [ec, nc] = m_in.get_floats(std::tie(r, g, b), ' ');
						if (ec == std::errc{}) {
							mtl.set_Kd(r, g, b);
						} else {
							OMLOG_ERROR << "failed to parse Kd for " << mtl.name();
						}
					} else if (keyword == "Ka") {
						float r = 0, g = 0, b = 0;
						auto [ec, nc] = m_in.get_floats(std::tie(r, g, b), ' ');
						if (ec == std::errc{}) {
							mtl.set_Ka(r, g, b);
						} else {
							OMLOG_ERROR << "failed to parse Ka for " << mtl.name();
						}
					} else if (keyword == "Ks") {
						float r = 0, g = 0, b = 0;
						auto [ec, nc] = m_in.get_floats(std::tie(r, g, b), ' ');
						if (ec == std::errc{}) {
							mtl.set_Ks(r, g, b);
						} else {
							OMLOG_ERROR << "failed to parse Ks for " << mtl.name();
						}
					} else if (keyword == "Tr" || keyword == "d") {
						float t = 0;
						auto [ec, nc] = m_in.get_floats(std::tie(t), ' ');
						if (ec == std::errc{}) {
							mtl.set_Tr(t);
						} else {
							OMLOG_ERROR << "failed to parse " << keyword << " for " << mtl.name();
						}
					} else if (keyword == "map_Kd") {
						std::string texname{iob::strip_any(m_in.get_until_any("\r\n"))};
						if (!texname.empty()) mtl.set_map_Kd(texname, textureId++);
					}
					// clean up end of line
					m_in.skip_line();
				}
				if (mtl.is_valid()) mtls.push_back(std::move(mtl));
				return mtls;
			}
		};

		class OBJParser {
		private:
			iob::text_reader m_in;
			BaseImporter *m_bi = nullptr;
			std::filesystem::path m_mtldir;
			std::map<std::string, Material> m_mtls;
			Material m_cur_mtl;

			// index tracking
			int m_n_v = 0;
			int m_n_vt = 0;
			int m_n_vn = 0;

			// mesh data
			std::vector<Vec3f>        m_normals;
			std::vector<Vec2f>        m_texcoords2d;
			std::vector<Vec3f>        m_texcoords3d;

			// face data (stored here to avoid reallocating)
			BaseImporter::VHandles    m_face_vhandles;
			std::vector<Vec3f>        m_face_normals;
			std::vector<Vec2f>        m_face_texcoords2d;
			std::vector<Vec3f>        m_face_texcoords3d;

		public:
			OBJParser() = default;

			OBJParser(const OBJParser &) = delete;
			OBJParser & operator=(const OBJParser &) = delete;

			OBJParser(iob::iobuffer *buf_, BaseImporter *bi_, std::filesystem::path mtldir_) :
				m_in{buf_}, m_bi{bi_}, m_mtldir{std::move(mtldir_)}
			{}

			bool parse() {
				// TODO report errors with line numbers?
				while (m_in.good()) {
					// skip leading whitespace (including blank lines)
					m_in.get_while_ws();
					// skip comments
					if (m_in.has_any('#')) {
						// TODO maybe parse vertex/face count hints?
						m_in.get_line();
						continue;
					}
					std::string keyword{m_in.get_until_ws()};
					// TODO other keywords like o and s
					if (keyword == "v") {
						parse_v();
						m_n_v++;
					} else if (keyword == "vt") {
						parse_vt();
						m_n_vt++;
					} else if (keyword == "vn") {
						parse_vn();
						m_n_vn++;
					} else if (keyword == "g") {
						parse_g();
					} else if (keyword == "f") {
						parse_f();
					} else if (keyword == "mtllib") {
						parse_mtllib();
					} else if (keyword == "usemtl") {
						parse_usemtl();
					}
					// clean up end of line
					m_in.skip_line();
				}
				if (m_bi->n_faces() == 0) treat_as_point_cloud();
				// TODO errors?
				return true;
			}

		private:
			void parse_v() {
				m_n_v++;
				float x = 0, y = 0, z = 0;
				auto [ec, nc] = m_in.get_floats(std::tie(x, y, z), ' ');
				if (nc == 3) {
					// try parse color
					auto vh = m_bi->add_vertex(Vec3f(x, y, z));
					float r = 0, g = 0, b = 0;
					auto [ec, nc] = m_in.get_floats(std::tie(r, g, b), ' ');
					if (nc == 3) {
						if (!!m_bi->request_vattribs(AttributeBits::Color)) {
							m_bi->set_color(vh, Vec3f(r, g, b));
						}
					}
				} else {
					OMLOG_ERROR << "failed to parse vertex position";
				}
			}

			void parse_vt() {
				m_n_vt++;
				float u = 0, v = 0, w = 0;
				auto [ec, nc] = m_in.get_floats(std::tie(u, v, w), ' ');
				if (nc == 3) {
					// 3d
					if (!!m_bi->request_hattribs(AttributeBits::TexCoord3D)) {
						m_texcoords3d.push_back(Vec3f(u, v, w));
					}
				} else if (nc == 2) {
					// 2d
					if (!!m_bi->request_hattribs(AttributeBits::TexCoord2D)) {
						m_texcoords2d.push_back(Vec2f(u, v));
					}
				} else {
					OMLOG_ERROR << "failed to parse texcoord";
				}
			}

			void parse_vn() {
				m_n_vn++;
				float x = 0, y = 0, z = 0;
				auto [ec, nc] = m_in.get_floats(std::tie(x, y, z), ' ');
				if (ec == std::errc{}) {
					if (!!m_bi->request_hattribs(AttributeBits::Normal)) {
						m_normals.push_back(Vec3f(x, y, z));
					}
				} else {
					OMLOG_ERROR << "failed to parse normal";
				}
			}

			void parse_g() {
				// TODO face groups
			}

			void parse_f() {
				m_face_vhandles.clear();
				m_face_normals.clear();
				m_face_texcoords2d.clear();
				m_face_texcoords3d.clear();
				m_in.get_while_any(' ');
				std::errc ec{};
				while (ec == std::errc{} && !m_in.has_any("\r\n")) {
					// parse v/vt/vn
					int v = 0, vt = 0, vn = 0;
					ec = m_in.get_int(v);
					if (ec == std::errc{} && m_in.get_any('/').size()) {
						if (m_in.get_any('/').size()) {
							ec = m_in.get_int(vn);
						} else {
							ec = m_in.get_int(vt);
							if (ec == std::errc{} && m_in.get_any('/').size()) {
								ec = m_in.get_int(vn);
							}
						}
					}
					if (ec != std::errc{}) {
						// check eof in case line was not terminated
						if (m_in.eof()) break;
						OMLOG_ERROR << "failed to parse face";
						return;
					}
					m_in.get_while_any(' ');
					// resolve negative indices and convert from 1-based to 0-based
					v = v < 0 ? (m_n_v + v) : (v - 1);
					vt = vt < 0 ? (m_n_vt + vt) : (vt - 1);
					vn = vn < 0 ? (m_n_vn + vn) : (vn - 1);
					// retrieve property values for this face
					if (0 <= v && v < int(m_bi->n_vertices())) {
						VertexHandle vh{v};
						Vec3f n{0};
						Vec2f tc2{0};
						Vec3f tc3{0};
						if (0 <= vn && vn < int(m_normals.size())) n = m_normals[vn];
						if (0 <= vt && vt < int(m_texcoords2d.size())) tc2 = m_texcoords2d[vt];
						if (0 <= vt && vt < int(m_texcoords3d.size())) tc3 = m_texcoords3d[vt];
						m_face_vhandles.push_back(vh);
						m_face_normals.push_back(n);
						m_face_texcoords2d.push_back(tc2);
						m_face_texcoords3d.push_back(tc3);
					} else {
						OMLOG_ERROR << "face vertex id out of range";
						// TODO possibly vertex defined after face, delay face parsing?
						// (this is why the original OpenMesh parser pessimistically did two passes)
					}
				}
				assert(ec == std::errc{});
				if (m_face_vhandles.size() > 2) {
					// note that add_face can possibly triangulate the faces, which is why we have to
					// store the current number of faces first
					size_t n_faces = m_bi->n_faces();
					FaceHandle fh = m_bi->add_face(m_face_vhandles);
					if (fh.is_valid()) {
						m_bi->add_face_normals(fh, m_face_vhandles[0], m_face_normals);
						m_bi->add_face_texcoords(fh, m_face_vhandles[0], m_face_texcoords2d);
						m_bi->add_face_texcoords(fh, m_face_vhandles[0], m_face_texcoords3d);
					}
					// set mtl color & texindex on new faces
					if (m_cur_mtl.has_Kd()) m_bi->request_fattribs(AttributeBits::Color);
					if (m_cur_mtl.has_map_Kd()) m_bi->request_fattribs(AttributeBits::TextureIndex);
					for (size_t i = n_faces; i < m_bi->n_faces(); ++i) {
						FaceHandle fh{int(i)};
						if (m_cur_mtl.is_valid()) {
							if (m_cur_mtl.has_Kd()) {
								// apply mtl color as face color
								m_bi->set_color(fh, m_cur_mtl.Kd());
							}
							if (m_cur_mtl.has_map_Kd()) {
								// Set the texture index in the face index property
								m_bi->set_face_texindex(fh, m_cur_mtl.map_Kd_index());
							} else {
								// If we don't have the info, set it to no texture
								m_bi->set_face_texindex(fh, 0);
							}
						} else {
							// Set the texture index to zero as we don't have any information
							m_bi->set_face_texindex(fh, 0);
						}
					}
				} else {
					OMLOG_ERROR << "face has < 3 vertices, skipping";
				}
			}

			void parse_mtllib() {
				const auto fpath = m_mtldir / std::filesystem::u8path(iob::strip_any(m_in.get_until_any("\r\n")));
				iob::file_buffer fbuf{fpath, iob::file_buffer::read};
				if (fbuf.is_open()) {
					MTLParser p{&fbuf};
					auto mtls = p.parse();
					for (auto &mtl : mtls) {
						m_mtls[mtl.name()] = std::move(mtl);
					}
				} else {
					OMLOG_ERROR << "failed to open mtl file " << fpath.u8string();
				}
			}

			void parse_usemtl() {
				m_in.get_while_ws();
				std::string mtlname{m_in.get_until_ws()};
				auto it = m_mtls.find(mtlname);
				if (it != m_mtls.end()) {
					m_cur_mtl = it->second;
				} else {
					m_cur_mtl = {};
					OMLOG_ERROR << "material " << mtlname << " not defined";
				}
			}

			void treat_as_point_cloud() {
				// add normal per vertex
				if (m_normals.size() == m_bi->n_vertices()) {
					if (!!m_bi->request_vattribs(AttributeBits::Normal)) {
						for (int i = 0; i < int(m_bi->n_vertices()); i++) {
							m_bi->set_normal(VertexHandle(i), m_normals[i]);
						}
					}
				}
			}
		};

		//-----------------------------------------------------------------------------

		_OBJReader_::
			_OBJReader_()
		{
			IOManager().register_module(this);
		}


		//-----------------------------------------------------------------------------


		bool _OBJReader_::read(const std::filesystem::path &_filename, BaseImporter &_bi) {
			iob::file_buffer fbuf{_filename, iob::file_buffer::read};

			if (!fbuf.is_open()) {
				OMLOG_ERROR << "Failed to open file " << _filename.u8string();
				return false;
			}

			OBJParser p{&fbuf, &_bi, _filename.parent_path()};
			return p.parse();
		}

		//-----------------------------------------------------------------------------

		bool _OBJReader_::read(std::istream &_in, BaseImporter &_bi) {
			iob::stream_buffer sbuf{_in.rdbuf()};
			OBJParser p{&sbuf, &_bi, ""};
			return p.parse();
		}


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
