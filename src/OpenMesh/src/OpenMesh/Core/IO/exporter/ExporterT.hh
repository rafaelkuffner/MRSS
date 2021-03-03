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
 //  Implements an exporter module for arbitrary OpenMesh meshes
 //
 //=============================================================================


#ifndef __EXPORTERT_HH__
#define __EXPORTERT_HH__


//=== INCLUDES ================================================================

// C++
#include <vector>

// OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/IO/exporter/BaseExporter.hh>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

//=== NAMESPACES ==============================================================

namespace OpenMesh {
	namespace IO {


		//=== EXPORTER CLASS ==========================================================

		/**
		 *  This class template provides an exporter module for OpenMesh meshes.
		 */
		template <class Mesh>
		class ExporterT : public BaseExporter
		{
		public:

			// Constructor
			ExporterT(const Mesh &_mesh, const Options &_useropts) : BaseExporter(_useropts), mesh_(_mesh)
			{
				// mask user options by available mesh attribs
				fileopts_ = _useropts;
				fileopts_.vattribs &= _mesh.vattribs();
				fileopts_.hattribs &= _mesh.hattribs();
				fileopts_.eattribs &= _mesh.eattribs();
				fileopts_.fattribs &= _mesh.fattribs();
			}

			virtual Vec3f point(VertexHandle _vh) const override
			{
				return mesh_.point(_vh);
			}

			virtual Vec3f normal(VertexHandle _vh) const override
			{
				Vec3f r{0};
				if (fileopts_.vertex_has_normal()) r = vector_cast<Vec3f>(mesh_.normal(_vh));
				else if (fileopts_.halfedge_has_normal()) r = vector_cast<Vec3f>(mesh_.normal(mesh_.halfedge_handle(_vh)));
				return r;
			}

			virtual Vec4uc colorA(VertexHandle _vh) const override
			{
				Vec4uc r{0};
				if (fileopts_.vertex_has_color()) r = color_cast<Vec4uc>(mesh_.color(_vh));
				else if (fileopts_.halfedge_has_color()) r = color_cast<Vec4uc>(mesh_.color(mesh_.halfedge_handle(_vh)));
				return r;
			}

			virtual Vec4f colorAf(VertexHandle _vh) const override
			{
				Vec4f r{0};
				if (fileopts_.vertex_has_color()) r = color_cast<Vec4f>(mesh_.color(_vh));
				else if (fileopts_.halfedge_has_color()) r = color_cast<Vec4f>(mesh_.color(mesh_.halfedge_handle(_vh)));
				return r;
			}

			virtual Vec2f texcoord2D(VertexHandle _vh) const override
			{
				Vec2f r{0};
				if (fileopts_.vertex_has_texcoord2D()) r = vector_cast<Vec2f>(mesh_.texcoord2D(_vh));
				else if (fileopts_.halfedge_has_texcoord2D()) r = vector_cast<Vec2f>(mesh_.texcoord2D(mesh_.halfedge_handle(_vh)));
				return r;
			}

			virtual Vec3f texcoord3D(VertexHandle _vh) const  override
			{
				Vec3f r{0};
				if (fileopts_.vertex_has_texcoord3D()) r = vector_cast<Vec3f>(mesh_.texcoord3D(_vh));
				else if (fileopts_.halfedge_has_texcoord3D()) r = vector_cast<Vec3f>(mesh_.texcoord3D(mesh_.halfedge_handle(_vh)));
				return r;
			}

			virtual StatusInfo status(VertexHandle _vh) const override
			{
				StatusInfo r{};
				if (fileopts_.vertex_has_status()) r = mesh_.status(_vh);
				return r;
			}

			virtual Vec3f normal(HalfedgeHandle _heh) const override
			{
				Vec3f r{0};
				if (fileopts_.halfedge_has_normal()) r = vector_cast<Vec3f>(mesh_.normal(_heh));
				else if (fileopts_.vertex_has_normal()) r = vector_cast<Vec3f>(mesh_.normal(mesh_.to_vertex_handle(_heh)));
				return r;
			}

			virtual Vec4uc colorA(HalfedgeHandle _heh) const override
			{
				Vec4uc r{0};
				if (fileopts_.halfedge_has_color()) r = color_cast<Vec4uc>(mesh_.color(_heh));
				else if (fileopts_.vertex_has_color()) r = color_cast<Vec4uc>(mesh_.color(mesh_.to_vertex_handle(_heh)));
				return r;
			}

			virtual Vec4f colorAf(HalfedgeHandle _heh) const override
			{
				Vec4f r{0};
				if (fileopts_.halfedge_has_color()) r = color_cast<Vec4f>(mesh_.color(_heh));
				else if (fileopts_.vertex_has_color()) r = color_cast<Vec4f>(mesh_.color(mesh_.to_vertex_handle(_heh)));
				return r;
			}

			virtual Vec2f texcoord2D(HalfedgeHandle _heh) const override
			{
				Vec2f r{0};
				if (fileopts_.halfedge_has_texcoord2D()) r = vector_cast<Vec2f>(mesh_.texcoord2D(_heh));
				else if (fileopts_.vertex_has_texcoord2D()) r = vector_cast<Vec2f>(mesh_.texcoord2D(mesh_.to_vertex_handle(_heh)));
				return r;
			}

			virtual Vec3f texcoord3D(HalfedgeHandle _heh) const override
			{
				Vec3f r{0};
				if (fileopts_.halfedge_has_texcoord3D()) r = vector_cast<Vec3f>(mesh_.texcoord3D(_heh));
				else if (fileopts_.vertex_has_texcoord3D()) r = vector_cast<Vec3f>(mesh_.texcoord3D(mesh_.to_vertex_handle(_heh)));
				return r;
			}

			virtual StatusInfo status(HalfedgeHandle _heh) const override
			{
				StatusInfo r{};
				if (fileopts_.halfedge_has_status()) r = mesh_.status(_heh);
				return r;
			}

			virtual Vec4uc colorA(EdgeHandle _eh) const override
			{
				Vec4uc r{0};
				if (fileopts_.edge_has_color()) r = color_cast<Vec4uc>(mesh_.color(_eh));
				return r;
			}

			virtual Vec4f colorAf(EdgeHandle _eh) const override
			{
				Vec4f r{0};
				if (fileopts_.edge_has_color()) r = color_cast<Vec4f>(mesh_.color(_eh));
				return r;
			}

			virtual StatusInfo status(EdgeHandle _eh) const override
			{
				StatusInfo r{};
				if (fileopts_.edge_has_status()) r = mesh_.status(_eh);
				return r;
			}

			virtual Vec3f normal(FaceHandle _fh) const override
			{
				Vec3f r{0};
				if (fileopts_.face_has_normal()) r = vector_cast<Vec3f>(mesh_.normal(_fh));
				return r;
			}

			virtual Vec4uc colorA(FaceHandle _fh) const override
			{
				Vec4uc r{0};
				if (fileopts_.face_has_color()) r = color_cast<Vec4uc>(mesh_.color(_fh));
				return r;
			}

			virtual Vec4f colorAf(FaceHandle _fh) const override
			{
				Vec4f r{0};
				if (fileopts_.face_has_color()) r = color_cast<Vec4f>(mesh_.color(_fh));
				return r;
			}

			virtual int texindex(FaceHandle _fh) const override
			{
				int r = 0;
				if (fileopts_.face_has_texindex()) r = mesh_.texture_index(_fh);
				return r;
			}

			virtual StatusInfo status(FaceHandle _fh) const override
			{
				StatusInfo r{};
				if (fileopts_.face_has_status()) r = mesh_.status(_fh);
				return r;
			}

			virtual HalfedgeHandle halfedge_handle(FaceHandle _fh) const override
			{
				return mesh_.halfedge_handle(_fh);
			}

			virtual HalfedgeHandle halfedge_handle(VertexHandle _vh) const override
			{
				// note: outgoing
				return mesh_.halfedge_handle(_vh);
			}

			virtual HalfedgeHandle next_halfedge_handle(HalfedgeHandle _heh) const override
			{
				return mesh_.next_halfedge_handle(_heh);
			}

			virtual HalfedgeHandle opposite_halfedge_handle(HalfedgeHandle _heh) const override
			{
				return mesh_.opposite_halfedge_handle(_heh);
			}

			virtual VertexHandle to_vertex_handle(HalfedgeHandle _heh) const override
			{
				return mesh_.to_vertex_handle(_heh);
			}

			virtual FaceHandle face_handle(HalfedgeHandle _heh) const override
			{
				return mesh_.face_handle(_heh);
			}

			virtual size_t face_vertex_handles(FaceHandle _fh, std::vector<VertexHandle> &_vhandles) const override
			{
				_vhandles.clear();
				for (auto it = mesh_.cfv_begin(_fh); it.is_valid(); ++it) {
					_vhandles.push_back(*it);
				}
				return _vhandles.size();
			}

			virtual size_t face_halfedge_handles(FaceHandle _fh, std::vector<HalfedgeHandle> &_hhandles) const override
			{
				_hhandles.clear();
				for (auto it = mesh_.cfh_begin(_fh); it.is_valid(); ++it) {
					_hhandles.push_back(*it);
				}
				return _hhandles.size();
			}

			virtual size_t vertex_halfedge_handles(VertexHandle _vh, std::vector<HalfedgeHandle> &_hhandles) const override
			{
				// note: outgoing
				_hhandles.clear();
				for (auto it = mesh_.cvoh_begin(_vh); it.is_valid(); ++it) {
					_hhandles.push_back(*it);
				}
				return _hhandles.size();
			}

			virtual const BaseKernel *kernel() override { return &mesh_; }

			// number of faces, vertices, edges
			size_t n_vertices()  const override { return mesh_.n_vertices(); }
			size_t n_faces()     const override { return mesh_.n_faces(); }
			size_t n_edges()     const override { return mesh_.n_edges(); }

			// mesh information
			bool is_triangle_mesh() const override
			{
				return Mesh::is_triangles();
			}

		private:

			const Mesh &mesh_;
		};


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
