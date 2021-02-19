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




 //=============================================================================
 //
 //  Implements an importer module for arbitrary OpenMesh meshes
 //
 //=============================================================================


#ifndef __IMPORTERT_HH__
#define __IMPORTERT_HH__


//=== INCLUDES ================================================================


#include <OpenMesh/Core/IO/importer/BaseImporter.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Core/System/omstream.hh>

#define OMLOG_SOURCE Importer

//== NAMESPACES ===============================================================


namespace OpenMesh {
	namespace IO {


		//=== IMPLEMENTATION ==========================================================


		/**
		 *  This class template provides an importer module for OpenMesh meshes.
		 */
		template <class Mesh>
		class ImporterT : public BaseImporter
		{
		public:

			typedef typename Mesh::Point       Point;
			typedef typename Mesh::Normal      Normal;
			typedef typename Mesh::Color       Color;
			typedef typename Mesh::TexCoord2D  TexCoord2D;
			typedef typename Mesh::TexCoord3D  TexCoord3D;
			typedef std::vector<VertexHandle>  VHandles;


			ImporterT(Mesh &_mesh, Options _useropts) : BaseImporter(_useropts), mesh_(_mesh) {}

			virtual void make_vattribs_impl(AttributeBits attribs) override
			{
				OMLOG_INFO << "making vertex attribs " << to_string(attribs);
				mesh_.request_vattribs(attribs);
			}

			virtual void make_eattribs_impl(AttributeBits attribs) override
			{
				OMLOG_INFO << "making edge attribs " << to_string(attribs);
				mesh_.request_eattribs(attribs);
			}

			virtual void make_hattribs_impl(AttributeBits attribs) override
			{
				OMLOG_INFO << "making halfedge attribs " << to_string(attribs);
				mesh_.request_hattribs(attribs);
			}

			virtual void make_fattribs_impl(AttributeBits attribs) override
			{
				OMLOG_INFO << "making face attribs " << to_string(attribs);
				mesh_.request_fattribs(attribs);
			}

			virtual VertexHandle add_vertex(const Vec3f &_point) override
			{
				return mesh_.add_vertex(vector_cast<Point>(_point));
			}

			virtual VertexHandle add_vertex() override
			{
				return mesh_.new_vertex();
			}

			virtual HalfedgeHandle add_edge(VertexHandle _vh0, VertexHandle _vh1) override
			{
				return mesh_.new_edge(_vh0, _vh1);
			}

			virtual FaceHandle add_face(const VHandles &_indices) override
			{
				FaceHandle fh;

				if (_indices.size() > 2)
				{
					VHandles::const_iterator it, it2, end(_indices.end());


					// test for valid vertex indices
					for (it = _indices.begin(); it != end; ++it)
						if (! mesh_.is_valid_handle(*it))
						{
							OMLOG_WARNING << "Face contains invalid vertex index";
							return fh;
						}


					// don't allow double vertices
					for (it = _indices.begin(); it != end; ++it)
						for (it2 = it + 1; it2 != end; ++it2)
							if (*it == *it2)
							{
								OMLOG_WARNING << "Face has equal vertices";
								return fh;
							}

					// try to add face
					fh = mesh_.add_face(_indices);

					// separate non-manifold faces and mark them
					if (!fh.is_valid())
					{
						VHandles vhandles(_indices.size());

						// double vertices
						for (unsigned int j = 0; j < _indices.size(); ++j)
						{
							// DO STORE p, reference may not work since vertex array
							// may be relocated after adding a new vertex !
							Point p = mesh_.point(_indices[j]);
							vhandles[j] = mesh_.add_vertex(p);

							// Mark vertices of failed face as non-manifold
							if (mesh_.has_vertex_status()) {
								mesh_.status(vhandles[j]).set_fixed_nonmanifold(true);
							}
						}

						// add face
						fh = mesh_.add_face(vhandles);

						// Mark failed face as non-manifold
						if (mesh_.has_face_status())
							mesh_.status(fh).set_fixed_nonmanifold(true);

						// Mark edges of failed face as non-two-manifold
						if (mesh_.has_edge_status()) {
							typename Mesh::FaceEdgeIter fe_it = mesh_.fe_iter(fh);
							for (; fe_it.is_valid(); ++fe_it) {
								mesh_.status(*fe_it).set_fixed_nonmanifold(true);
							}
						}
					}

					synthesize_face_hattribs(fh);
				}
				return fh;
			}

			virtual FaceHandle add_face(HalfedgeHandle _heh) override
			{
				auto fh = mesh_.new_face();
				mesh_.set_halfedge_handle(fh, _heh);
				return fh;
			}

			// vertex attributes

			virtual void set_point(VertexHandle _vh, const Vec3f &_point) override
			{
				mesh_.set_point(_vh, vector_cast<Point>(_point));
			}

			virtual void set_halfedge(VertexHandle _vh, HalfedgeHandle _heh) override
			{
				mesh_.set_halfedge_handle(_vh, _heh);
			}

			virtual void set_normal(VertexHandle _vh, const Vec3f &_normal) override
			{
				if (fileopts_.vertex_has_normal()) {
					mesh_.set_normal(_vh, vector_cast<Normal>(_normal));
				} else if (synthopts_.halfedge_has_normal()) {
					if (!vprop_fallback_normal_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex normals";
						mesh_.add_property(vprop_fallback_normal_);
					}
					mesh_.property(vprop_fallback_normal_, _vh) = vector_cast<Normal>(_normal);
				}
			}

			virtual void set_color(VertexHandle _vh, const Vec4uc &_color) override
			{
				if (fileopts_.vertex_has_color()) {
					mesh_.set_color(_vh, color_cast<Color>(_color));
				} else if (synthopts_.halfedge_has_color()) {
					if (!vprop_fallback_color_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex colors";
						mesh_.add_property(vprop_fallback_color_);
					}
					mesh_.property(vprop_fallback_color_, _vh) = color_cast<Color>(_color);
				}
			}

			virtual void set_color(VertexHandle _vh, const Vec3uc &_color) override
			{
				if (fileopts_.vertex_has_color()) {
					mesh_.set_color(_vh, color_cast<Color>(_color));
				} else if (synthopts_.halfedge_has_color()) {
					if (!vprop_fallback_color_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex colors";
						mesh_.add_property(vprop_fallback_color_);
					}
					mesh_.property(vprop_fallback_color_, _vh) = color_cast<Color>(_color);
				}
			}

			virtual void set_color(VertexHandle _vh, const Vec4f &_color) override
			{
				if (fileopts_.vertex_has_color()) {
					mesh_.set_color(_vh, color_cast<Color>(_color));
				} else if (synthopts_.halfedge_has_color()) {
					if (!vprop_fallback_color_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex colors";
						mesh_.add_property(vprop_fallback_color_);
					}
					mesh_.property(vprop_fallback_color_, _vh) = color_cast<Color>(_color);
				}
			}

			virtual void set_color(VertexHandle _vh, const Vec3f &_color) override
			{
				if (fileopts_.vertex_has_color()) {
					mesh_.set_color(_vh, color_cast<Color>(_color));
				} else if (synthopts_.halfedge_has_color()) {
					if (!vprop_fallback_color_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex colors";
						mesh_.add_property(vprop_fallback_color_);
					}
					mesh_.property(vprop_fallback_color_, _vh) = color_cast<Color>(_color);
				}
			}

			virtual void set_texcoord(VertexHandle _vh, const Vec2f &_texcoord) override
			{
				if (fileopts_.vertex_has_texcoord2D()) {
					mesh_.set_texcoord2D(_vh, vector_cast<TexCoord2D>(_texcoord));
				} else if (synthopts_.halfedge_has_texcoord2D()) {
					if (!vprop_fallback_texcoord2D_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex texcoords2D";
						mesh_.add_property(vprop_fallback_texcoord2D_);
					}
					mesh_.property(vprop_fallback_texcoord2D_, _vh) = vector_cast<TexCoord2D>(_texcoord);
				}
			}

			virtual void set_texcoord(VertexHandle _vh, const Vec3f &_texcoord) override
			{
				if (fileopts_.vertex_has_texcoord3D()) {
					mesh_.set_texcoord3D(_vh, vector_cast<TexCoord3D>(_texcoord));
				} else if (synthopts_.halfedge_has_texcoord3D()) {
					if (!vprop_fallback_texcoord3D_.is_valid()) {
						OMLOG_DEBUG << "making fallback vertex texcoords3d";
						mesh_.add_property(vprop_fallback_texcoord3D_);
					}
					mesh_.property(vprop_fallback_texcoord3D_, _vh) = vector_cast<TexCoord3D>(_texcoord);
				}
			}

			virtual void set_status(VertexHandle _vh, const OpenMesh::Attributes::StatusInfo &_status) override
			{
				if (!mesh_.has_vertex_status())
					mesh_.request_vertex_status();
				mesh_.status(_vh) = _status;
			}

			virtual void set_next(HalfedgeHandle _heh, HalfedgeHandle _next) override
			{
				mesh_.set_next_halfedge_handle(_heh, _next);
			}

			virtual void set_face(HalfedgeHandle _heh, FaceHandle _fh) override
			{
				mesh_.set_face_handle(_heh, _fh);
			}


			virtual void set_texcoord(HalfedgeHandle _heh, const Vec2f &_texcoord) override
			{
				if (fileopts_.halfedge_has_texcoord2D())
					mesh_.set_texcoord2D(_heh, vector_cast<TexCoord2D>(_texcoord));
				if (synthopts_.vertex_has_texcoord2D())
					mesh_.set_texcoord2D(mesh_.to_vertex_handle(_heh), vector_cast<TexCoord2D>(_texcoord));
			}

			virtual void set_texcoord(HalfedgeHandle _heh, const Vec3f &_texcoord) override
			{
				if (fileopts_.halfedge_has_texcoord3D())
					mesh_.set_texcoord3D(_heh, vector_cast<TexCoord3D>(_texcoord));
				if (synthopts_.vertex_has_texcoord3D())
					mesh_.set_texcoord3D(mesh_.to_vertex_handle(_heh), vector_cast<TexCoord3D>(_texcoord));
			}

			virtual void set_normal(HalfedgeHandle _heh, const Vec3f &_normal) override {
				if (fileopts_.halfedge_has_normal())
					mesh_.set_normal(_heh, vector_cast<Normal>(_normal));
				if (synthopts_.vertex_has_normal())
					mesh_.set_normal(mesh_.to_vertex_handle(_heh), vector_cast<Normal>(_normal));
			}

			// TODO halfedge color?

			virtual void set_status(HalfedgeHandle _heh, const OpenMesh::Attributes::StatusInfo &_status) override
			{
				if (!mesh_.has_halfedge_status())
					mesh_.request_halfedge_status();
				mesh_.status(_heh) = _status;
			}

			// edge attributes

			virtual void set_color(EdgeHandle _eh, const Vec4uc &_color) override
			{
				if (fileopts_.edge_has_color())
					mesh_.set_color(_eh, color_cast<Color>(_color));
			}

			virtual void set_color(EdgeHandle _eh, const Vec3uc &_color) override
			{
				if (fileopts_.edge_has_color())
					mesh_.set_color(_eh, color_cast<Color>(_color));
			}

			virtual void set_color(EdgeHandle _eh, const Vec4f &_color) override
			{
				if (fileopts_.edge_has_color())
					mesh_.set_color(_eh, color_cast<Color>(_color));
			}

			virtual void set_color(EdgeHandle _eh, const Vec3f &_color) override
			{
				if (fileopts_.edge_has_color())
					mesh_.set_color(_eh, color_cast<Color>(_color));
			}

			virtual void set_status(EdgeHandle _eh, const OpenMesh::Attributes::StatusInfo &_status) override
			{
				if (!mesh_.has_edge_status())
					mesh_.request_edge_status();
				mesh_.status(_eh) = _status;
			}

			// face attributes

			virtual void set_normal(FaceHandle _fh, const Vec3f &_normal) override
			{
				if (fileopts_.face_has_normal())
					mesh_.set_normal(_fh, vector_cast<Normal>(_normal));
			}

			virtual void set_color(FaceHandle _fh, const Vec3uc &_color) override
			{
				if (fileopts_.face_has_color())
					mesh_.set_color(_fh, color_cast<Color>(_color));
			}

			virtual void set_color(FaceHandle _fh, const Vec4uc &_color) override
			{
				if (fileopts_.face_has_color())
					mesh_.set_color(_fh, color_cast<Color>(_color));
			}

			virtual void set_color(FaceHandle _fh, const Vec3f &_color) override
			{
				if (fileopts_.face_has_color())
					mesh_.set_color(_fh, color_cast<Color>(_color));
			}

			virtual void set_color(FaceHandle _fh, const Vec4f &_color) override
			{
				if (fileopts_.face_has_color())
					mesh_.set_color(_fh, color_cast<Color>(_color));
			}

			virtual void set_status(FaceHandle _fh, const OpenMesh::Attributes::StatusInfo &_status) override
			{
				if (!mesh_.has_face_status())
					mesh_.request_face_status();
				mesh_.status(_fh) = _status;
			}

			virtual void add_face_texcoords(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec2f> &_face_texcoords) override
			{
				if (!want_hattribs(AttributeBits::TexCoord2D)) return;
				
				// get first halfedge handle
				HalfedgeHandle cur_heh = mesh_.halfedge_handle(_fh);
				HalfedgeHandle end_heh = mesh_.prev_halfedge_handle(cur_heh);

				// find start heh
				while (mesh_.to_vertex_handle(cur_heh) != _vh && cur_heh != end_heh)
					cur_heh = mesh_.next_halfedge_handle(cur_heh);

				for (unsigned int i = 0; i < _face_texcoords.size(); ++i)
				{
					set_texcoord(cur_heh, _face_texcoords[i]);
					cur_heh = mesh_.next_halfedge_handle(cur_heh);
				}
			}

			virtual void add_face_texcoords(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec3f> &_face_texcoords) override
			{
				if (!want_hattribs(AttributeBits::TexCoord3D)) return;
				
				// get first halfedge handle
				HalfedgeHandle cur_heh = mesh_.halfedge_handle(_fh);
				HalfedgeHandle end_heh = mesh_.prev_halfedge_handle(cur_heh);

				// find start heh
				while (mesh_.to_vertex_handle(cur_heh) != _vh && cur_heh != end_heh)
					cur_heh = mesh_.next_halfedge_handle(cur_heh);

				for (unsigned int i = 0; i < _face_texcoords.size(); ++i)
				{
					set_texcoord(cur_heh, _face_texcoords[i]);
					cur_heh = mesh_.next_halfedge_handle(cur_heh);
				}
			}

			virtual void add_face_normals(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec3f> &_face_normals) override
			{
				if (!want_hattribs(AttributeBits::Normal)) return;
				
				// get first halfedge handle
				HalfedgeHandle cur_heh = mesh_.halfedge_handle(_fh);
				HalfedgeHandle end_heh = mesh_.prev_halfedge_handle(cur_heh);

				// find start heh
				while (mesh_.to_vertex_handle(cur_heh) != _vh && cur_heh != end_heh)
					cur_heh = mesh_.next_halfedge_handle(cur_heh);

				for (unsigned int i = 0; i < _face_normals.size(); ++i)
				{
					set_normal(cur_heh, _face_normals[i]);
					cur_heh = mesh_.next_halfedge_handle(cur_heh);
				}
			}

			void synthesize_face_hattribs(FaceHandle fh)
			{
				if (!fh.is_valid() || !synthopts_.hattribs) return;
				// get vertex props to use for synthesis
				// note: could theoretically have valid fallback props but synthesis disabled
				// note: vertex attribs could exist without the importer being allowed to touch them
				auto vprop_normal = synthopts_.halfedge_has_normal()
					? (fileopts_.vertex_has_normal() ? mesh_.vertex_normals_pph() : vprop_fallback_normal_)
					: VPropHandleT<Normal>{};
				auto vprop_color = synthopts_.halfedge_has_color()
					? (fileopts_.vertex_has_color() ? mesh_.vertex_colors_pph() : vprop_fallback_color_)
					: VPropHandleT<Color>{};
				auto vprop_texcoord2D = synthopts_.halfedge_has_texcoord2D()
					? (fileopts_.vertex_has_texcoord2D() ? mesh_.vertex_texcoords2D_pph() : vprop_fallback_texcoord2D_)
					: VPropHandleT<TexCoord2D>{};
				auto vprop_texcoord3D = synthopts_.halfedge_has_texcoord3D()
					? (fileopts_.vertex_has_texcoord3D() ? mesh_.vertex_texcoords3D_pph() : vprop_fallback_texcoord3D_)
					: VPropHandleT<TexCoord3D>{};
				// synthesize props for all halfedges of the face
				for (auto hit = mesh_.fh_iter(fh); hit.is_valid(); ++hit) {
					HalfedgeHandle hh = *hit;
					VertexHandle vh = mesh_.to_vertex_handle(hh);
					if (vprop_normal.is_valid()) mesh_.set_normal(hh, mesh_.property(vprop_normal, vh));
					if (vprop_color.is_valid()) mesh_.set_color(hh, mesh_.property(vprop_color, vh));
					if (vprop_texcoord2D.is_valid()) mesh_.set_texcoord2D(hh, mesh_.property(vprop_texcoord2D, vh));
					if (vprop_texcoord3D.is_valid()) mesh_.set_texcoord3D(hh, mesh_.property(vprop_texcoord3D, vh));
				}
			}

			virtual void set_face_texindex(FaceHandle _fh, int _texId) override
			{
				if (fileopts_.face_has_texindex()) {
					mesh_.set_texture_index(_fh, _texId);
				}
			}

			virtual void add_texture_information(int _id, std::string _name) override
			{
				// note: dont check fileopts for texture index

				OpenMesh::MPropHandleT< std::map< int, std::string > > property;

				if (!mesh_.get_property_handle(property, "TextureMapping")) {
					mesh_.add_property(property, "TextureMapping");
				}

				if (mesh_.property(property).find(_id) == mesh_.property(property).end())
					mesh_.property(property)[_id] = _name;
			}

			// low-level access to mesh

			virtual BaseKernel *kernel() override { return &mesh_; }

			virtual bool is_triangle_mesh() const override
			{
				return Mesh::is_triangles();
			}

			virtual void reserve(unsigned int nV, unsigned int nE, unsigned int nF) override
			{
				mesh_.reserve(nV, nE, nF);
			}

			// query number of faces, vertices, normals, texcoords
			virtual size_t n_vertices()  const override { return mesh_.n_vertices(); }
			virtual size_t n_faces()     const override { return mesh_.n_faces(); }
			virtual size_t n_edges()     const override { return mesh_.n_edges(); }


			virtual void prepare() override { }


			virtual void finish() override
			{
				// remove temp properties
				mesh_.remove_property(vprop_fallback_color_);
				mesh_.remove_property(vprop_fallback_normal_);
				mesh_.remove_property(vprop_fallback_texcoord2D_);
				mesh_.remove_property(vprop_fallback_texcoord3D_);
			}


		private:
			Mesh &mesh_;

			// fallback properties for vertex -> halfedge attribute conversion
			VPropHandleT<Normal> vprop_fallback_normal_;
			VPropHandleT<Color> vprop_fallback_color_;
			VPropHandleT<TexCoord2D> vprop_fallback_texcoord2D_;
			VPropHandleT<TexCoord3D> vprop_fallback_texcoord3D_;
		};


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================

#undef OMLOG_SOURCE
