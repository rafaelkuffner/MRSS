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
 //  Implements the baseclass for MeshWriter exporter modules
 //
 //=============================================================================


#ifndef __BASEEXPORTER_HH__
#define __BASEEXPORTER_HH__


//=== INCLUDES ================================================================


// STL
#include <vector>
#include <type_traits>
#include <utility>

// OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/BaseKernel.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>

//=== NAMESPACES ==============================================================


namespace OpenMesh {
	namespace IO {


		//=== EXPORTER ================================================================


		/**
		   Base class for exporter modules.
		   The exporter modules provide an interface between the writer modules and
		   the target data structure.
		*/

		class OPENMESHDLLEXPORT BaseExporter
		{
		protected:
			Options useropts_;
			Options fileopts_;

		public:
			using StatusInfo = OpenMesh::Attributes::StatusInfo;

			virtual ~BaseExporter() { }

			BaseExporter(Options _useropts) : useropts_(_useropts) {}

			// user options (wanted/allowed by user)
			const Options & user_options() const
			{
				return useropts_;
			}

			// options exported to file (subset of user options)
			const Options & file_options() const
			{
				return fileopts_;
			}

			// set file option flags
			void set_file_options(OptionBits opts) {
				fileopts_.flags |= opts;
			}

			// unset file option flags
			void unset_file_options(OptionBits opts) {
				fileopts_.flags &= ~opts;
			}

			// vertex properties
			virtual Vec3f point(VertexHandle _vh) const = 0;
			virtual Vec3f normal(VertexHandle _vh) const = 0;
			virtual Vec4uc colorA(VertexHandle _vh) const = 0;
			virtual Vec4f colorAf(VertexHandle _vh) const = 0;
			virtual Vec2f texcoord2D(VertexHandle _vh) const = 0;
			virtual Vec3f texcoord3D(VertexHandle _vh) const = 0;
			virtual StatusInfo status(VertexHandle _vh) const = 0;
			Vec3uc color(VertexHandle _vh) const { return color_cast<Vec3uc>(colorA(_vh)); }
			Vec3f colorf(VertexHandle _vh) const { return color_cast<Vec3f>(colorAf(_vh)); }

			// halfedge properties
			virtual Vec3f normal(HalfedgeHandle _heh) const = 0;
			virtual Vec4uc colorA(HalfedgeHandle _heh) const = 0;
			virtual Vec4f colorAf(HalfedgeHandle _heh) const = 0;
			virtual Vec2f texcoord2D(HalfedgeHandle _heh) const = 0;
			virtual Vec3f texcoord3D(HalfedgeHandle _heh) const = 0;
			virtual StatusInfo status(HalfedgeHandle _heh) const = 0;
			Vec3uc color(HalfedgeHandle _heh) const { return color_cast<Vec3uc>(colorA(_heh)); }
			Vec3f colorf(HalfedgeHandle _heh) const { return color_cast<Vec3f>(colorAf(_heh)); }

			// edge properties
			virtual Vec4uc colorA(EdgeHandle _eh) const = 0;
			virtual Vec4f colorAf(EdgeHandle _eh) const = 0;
			virtual StatusInfo status(EdgeHandle _eh) const = 0;
			Vec3uc color(EdgeHandle _eh) const { return color_cast<Vec3uc>(colorA(_eh)); }
			Vec3f colorf(EdgeHandle _eh) const { return color_cast<Vec3f>(colorAf(_eh)); }

			// face properties
			virtual Vec3f normal(FaceHandle _fh) const = 0;
			virtual Vec4uc colorA(FaceHandle _fh) const = 0;
			virtual Vec4f colorAf(FaceHandle _fh) const = 0;
			virtual int texindex(FaceHandle _fh) const = 0;
			virtual StatusInfo status(FaceHandle _fh) const = 0;
			Vec3uc color(FaceHandle _fh) const { return color_cast<Vec3uc>(colorA(_fh)); }
			Vec3f colorf(FaceHandle _fh) const { return color_cast<Vec3f>(colorAf(_fh)); }

			// connectivity
			virtual HalfedgeHandle halfedge_handle(FaceHandle _fh) const = 0;
			// outgoing
			virtual HalfedgeHandle halfedge_handle(VertexHandle _vh) const = 0;
			virtual HalfedgeHandle next_halfedge_handle(HalfedgeHandle _heh) const = 0;
			virtual HalfedgeHandle opposite_halfedge_handle(HalfedgeHandle _heh) const = 0;
			virtual VertexHandle to_vertex_handle(HalfedgeHandle _heh) const = 0;
			virtual FaceHandle face_handle(HalfedgeHandle _heh) const = 0;
			virtual size_t face_vertex_handles(FaceHandle _fh, std::vector<VertexHandle> &_vhandles) const = 0;
			virtual size_t face_halfedge_handles(FaceHandle _fh, std::vector<HalfedgeHandle> &_hhandles) const = 0;
			// outgoing
			virtual size_t vertex_halfedge_handles(VertexHandle _vh, std::vector<HalfedgeHandle> &_hhandles) const = 0;
			
			// get reference to base kernel
			virtual const BaseKernel * kernel() const { return 0; }

			// number of faces, vertices, edges
			virtual size_t n_vertices()   const = 0;
			virtual size_t n_faces()      const = 0;
			virtual size_t n_edges()      const = 0;

			// mesh information
			virtual bool is_triangle_mesh()     const { return false; }

			template <
				typename T,
				typename PackFn,
				typename OutFn,
				typename = std::enable_if_t<std::is_invocable_r_v<T, PackFn, const BaseExporter &, HalfedgeHandle>>,
				typename = std::enable_if_t<std::is_invocable_v<OutFn, int, T &&>>
			>
			void compress_halfedge_properties(std::vector<int> &indices, const int idx0, PackFn &&packfn, OutFn &&outfn) {
				const BaseKernel &mesh = *kernel();
				indices.resize(n_edges() * 2);
				int idx_next = idx0;
				std::vector<HalfedgeHandle> hhandles;
				std::vector<T> cache;
				for (int iv = 0; iv < n_vertices(); ++iv) {
					// search incoming halfedges at each vertex for duplicate values
					VertexHandle vh{iv};
					cache.clear();
					vertex_halfedge_handles(vh, hhandles);
					for (auto ohh : hhandles) {
						auto hh = opposite_halfedge_handle(ohh);
						T t = packfn(static_cast<const BaseExporter &>(*this), hh);
						int idx1 = -1;
						for (int j = 0; j < int(cache.size()); ++j) {
							if (t == cache[j]) {
								// use index of existing value
								idx1 = idx_next - int(cache.size()) + j;
								break;
							}
						}
						if (idx1 < 0) {
							// add a new value
							idx1 = idx_next++;
							cache.push_back(t);
							outfn(idx1, std::move(t));
						}
						indices[hh.idx()] = idx1;
					}
				}
			}
		};


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
