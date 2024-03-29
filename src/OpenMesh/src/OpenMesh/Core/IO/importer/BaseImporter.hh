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
 //  Implements the baseclass for IOManager importer modules
 //
 //=============================================================================


#ifndef __BASEIMPORTER_HH__
#define __BASEIMPORTER_HH__


//=== INCLUDES ================================================================


// STL
#include <vector>

// OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/BaseKernel.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>

#define OMLOG_SOURCE Importer

//== NAMESPACES ===============================================================


namespace OpenMesh {
	namespace IO {


		//=== IMPLEMENTATION ==========================================================


		/**  Base class for importer modules. Importer modules provide an
		 *   interface between the loader modules and the target data
		 *   structure. This is basically a wrapper providing virtual versions
		 *   for the required mesh functions.
		 */
		class OPENMESHDLLEXPORT BaseImporter
		{
		protected:
			Options useropts_{};
			Options fileopts_{};
			Options synthopts_{};

		public:

			// base class needs virtual destructor
			virtual ~BaseImporter() {}

			BaseImporter(Options _useropts) : useropts_(_useropts) {}

			// user options (wanted/allowed by user)
			const Options & user_options() const
			{
				return useropts_;
			}

			// options loaded directly from file (subset of user options)
			const Options & file_options() const
			{
				return fileopts_;
			}

			// options loaded from file including those synthesized by the importer (subset of user options)
			const Options final_options() const
			{
				return fileopts_ | synthopts_;
			}

			// set file option flags
			void set_file_options(OptionBits opts) {
				// not exposing write access to full options to maintain control of attribs with request_* functions
				fileopts_.flags |= opts;
			}

			// unset file option flags
			void unset_file_options(OptionBits opts) {
				fileopts_.flags &= ~opts;
			}

			// returns settable subset of queried attributes.
			AttributeBits want_vattribs(AttributeBits attribs) const
			{
				return attribs & (fileopts_.vattribs | synthopts_.hattribs);
			}

			// returns settable subset of queried attributes.
			AttributeBits want_hattribs(AttributeBits attribs) const
			{
				return attribs & (fileopts_.hattribs | synthopts_.vattribs);
			}

			// returns settable subset of queried attributes.
			AttributeBits want_eattribs(AttributeBits attribs) const
			{
				return attribs & fileopts_.eattribs;
			}

			// returns settable subset of queried attributes.
			AttributeBits want_fattribs(AttributeBits attribs) const
			{
				return attribs & fileopts_.fattribs;
			}

			// for reader to request vertex attribs, falling back to halfedge attribs, checked against user options.
			// if fallback to halfedge attribs happens, vertex values for those attribs are collected in temporary
			// properties and applied to halfedges when faces are created.
			// returns settable subset of requested attributes.
			AttributeBits request_vattribs(AttributeBits attribs)
			{
				// request vertex attribs masked by user opts
				const auto va = attribs & useropts_.vattribs & ~fileopts_.vattribs;
				if (!!va) make_vattribs_impl(va & ~synthopts_.vattribs);
				fileopts_.vattribs |= va;
				synthopts_.vattribs &= ~va;
				// and also halfedge attribs for synthesis
				const auto ha = attribs & useropts_.hattribs & ~fileopts_.hattribs & ~synthopts_.hattribs;
				if (!!ha) {
					OMLOG_INFO << "synthesizing halfedge attribs " << to_string(ha);
					make_hattribs_impl(ha);
					synthopts_.hattribs |= ha;
				}
				return attribs & (fileopts_.vattribs | synthopts_.hattribs);
			}

			// for reader to request halfedge attribs, falling back to vertex attribs, checked against user options.
			// if fallback to vertex attribs happens, halfedge values for those attribs are applied to vertices instead.
			// returns settable subset of requested attributes.
			AttributeBits request_hattribs(AttributeBits attribs)
			{
				// request halfedge attribs masked by user opts
				const auto ha = attribs & useropts_.hattribs & ~fileopts_.hattribs;
				if (!!ha) make_hattribs_impl(ha & ~synthopts_.hattribs);
				fileopts_.hattribs |= ha;
				synthopts_.hattribs &= ~ha;
				// and also vertex attribs for synthesis
				const auto va = attribs & useropts_.vattribs & ~fileopts_.vattribs & ~synthopts_.vattribs;
				if (!!va) {
					OMLOG_INFO << "synthesizing vertex attribs " << to_string(va);
					make_vattribs_impl(va);
					synthopts_.vattribs |= va;
				}
				return attribs & (fileopts_.hattribs | synthopts_.vattribs);
			}

			// for reader to request edge attribs checked against user options.
			// returns settable subset of requested attributes.
			AttributeBits request_eattribs(AttributeBits attribs)
			{
				AttributeBits r = useropts_.eattribs & attribs & ~fileopts_.eattribs;
				if (!!r) make_eattribs_impl(r);
				return attribs & (fileopts_.eattribs |= r);
			}

			// for reader to request face attribs checked against user options.
			// returns settable subset of requested attributes.
			AttributeBits request_fattribs(AttributeBits attribs)
			{
				AttributeBits r = useropts_.fattribs & attribs & ~fileopts_.fattribs;
				if (!!r) make_fattribs_impl(r);
				return attribs & (fileopts_.hattribs |= r);
			}

			// add a vertex with coordinate \c _point
			virtual VertexHandle add_vertex(const Vec3f &_point) = 0;

			// add a vertex without coordinate. Use set_point to set the position deferred
			virtual VertexHandle add_vertex() = 0;

			// add an edge. Use set_next, set_vertex and set_face to set corresponding entities for halfedges
			virtual HalfedgeHandle add_edge(VertexHandle _vh0, VertexHandle _vh1) = 0;

			// add a face with indices _indices refering to vertices
			typedef std::vector<VertexHandle> VHandles;
			virtual FaceHandle add_face(const VHandles &_indices) = 0;

			// add a face with incident halfedge
			virtual FaceHandle add_face(HalfedgeHandle _heh) = 0;

			// add texture coordinates per face, _vh references the first texcoord
			virtual void add_face_texcoords(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec2f> &_face_texcoords) = 0;

			// add texture 3d coordinates per face, _vh references the first texcoord
			virtual void add_face_texcoords(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec3f> &_face_texcoords) = 0;

			// add normals per face, _vh references the first normal
			virtual void add_face_normals(FaceHandle _fh, VertexHandle _vh, const std::vector<Vec3f> &_face_normals) = 0;

			// Set the texture index for a face
			virtual void set_face_texindex(FaceHandle _fh, int _texId) = 0;

			// Set coordinate of the given vertex. Use this function, if you created a vertex without coordinate
			virtual void set_point(VertexHandle _vh, const Vec3f &_point) = 0;

			// Set outgoing halfedge for the given vertex.
			virtual void set_halfedge(VertexHandle _vh, HalfedgeHandle _heh) = 0;

			// set vertex normal
			virtual void set_normal(VertexHandle _vh, const Vec3f &_normal) = 0;

			// set vertex color
			virtual void set_color(VertexHandle _vh, const Vec3uc &_color) = 0;

			// set vertex color
			virtual void set_color(VertexHandle _vh, const Vec4uc &_color) = 0;

			// set vertex color
			virtual void set_color(VertexHandle _vh, const Vec3f &_color) = 0;

			// set vertex color
			virtual void set_color(VertexHandle _vh, const Vec4f &_color) = 0;

			// set vertex texture coordinate
			virtual void set_texcoord(VertexHandle _vh, const Vec2f &_texcoord) = 0;

			// set vertex status
			virtual void set_status(VertexHandle _vh, const OpenMesh::Attributes::StatusInfo &_status) = 0;

			// set next halfedge handle
			virtual void set_next(HalfedgeHandle _heh, HalfedgeHandle _next) = 0;

			// set incident face handle for given halfedge
			virtual void set_face(HalfedgeHandle _heh, FaceHandle _fh) = 0;

			// set vertex texture coordinate
			virtual void set_texcoord(HalfedgeHandle _heh, const Vec2f &_texcoord) = 0;

			// set 3d vertex texture coordinate
			virtual void set_texcoord(VertexHandle _vh, const Vec3f &_texcoord) = 0;

			// set 3d vertex texture coordinate
			virtual void set_texcoord(HalfedgeHandle _heh, const Vec3f &_texcoord) = 0;

			// set halfedge normal
			virtual void set_normal(HalfedgeHandle _heh, const Vec3f &_normal) = 0;

			// set halfedge status
			virtual void set_status(HalfedgeHandle _heh, const OpenMesh::Attributes::StatusInfo &_status) = 0;

			// set edge color
			virtual void set_color(EdgeHandle _eh, const Vec3uc &_color) = 0;

			// set edge color
			virtual void set_color(EdgeHandle _eh, const Vec4uc &_color) = 0;

			// set edge color
			virtual void set_color(EdgeHandle _eh, const Vec3f &_color) = 0;

			// set edge color
			virtual void set_color(EdgeHandle _eh, const Vec4f &_color) = 0;

			// set edge status
			virtual void set_status(EdgeHandle _eh, const OpenMesh::Attributes::StatusInfo &_status) = 0;

			// set face normal
			virtual void set_normal(FaceHandle _fh, const Vec3f &_normal) = 0;

			// set face color
			virtual void set_color(FaceHandle _fh, const Vec3uc &_color) = 0;

			// set face color
			virtual void set_color(FaceHandle _fh, const Vec4uc &_color) = 0;

			// set face color
			virtual void set_color(FaceHandle _fh, const Vec3f &_color) = 0;

			// set face color
			virtual void set_color(FaceHandle _fh, const Vec4f &_color) = 0;

			// set face status
			virtual void set_status(FaceHandle _fh, const OpenMesh::Attributes::StatusInfo &_status) = 0;

			// Store a property in the mesh mapping from an int to a texture file
			// Use set_face_texindex to set the index for each face
			virtual void add_texture_information(int _id, std::string _name) = 0;

			// get reference to base kernel
			virtual BaseKernel *kernel() { return 0; }

			virtual bool is_triangle_mesh()     const { return false; }

			// reserve mem for elements
			virtual void reserve(unsigned int /* nV */,
				unsigned int /* nE */,
				unsigned int /* nF */) {}

			// query number of faces, vertices, normals, texcoords
			virtual size_t n_vertices()   const = 0;
			virtual size_t n_faces()      const = 0;
			virtual size_t n_edges()      const = 0;


			// pre-processing
			virtual void prepare() {}

			// post-processing
			virtual void finish() {}

		private:
			virtual void make_vattribs_impl(AttributeBits attribs) = 0;
			virtual void make_eattribs_impl(AttributeBits attribs) = 0;
			virtual void make_hattribs_impl(AttributeBits attribs) = 0;
			virtual void make_fattribs_impl(AttributeBits attribs) = 0;
		};


		//=============================================================================
	} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================

#undef OMLOG_SOURCE
