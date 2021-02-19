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




 /**
	 \file Attributes.hh
	 This file provides some macros containing attribute usage.
 */


#ifndef OPENMESH_ATTRIBUTES_HH
#define OPENMESH_ATTRIBUTES_HH


 //== INCLUDES =================================================================

#include <string>

#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>


//== NAMESPACES ===============================================================


namespace OpenMesh {

	//== CLASS DEFINITION  ========================================================

	/** Attribute bits
		*
		*  Use the bits to define a standard property at compile time using traits.
		*
		*  \include traits5.cc
		*
		*  \see \ref mesh_type
		*/
	enum class AttributeBits : unsigned
	{
		None = 0,  ///< Clear all attribute bits
		Normal = 1,  ///< Add normals to mesh item (vertices/faces)
		Color = 2,  ///< Add colors to mesh item (vertices/faces/edges)
		PrevHalfedge = 4,  ///< Add storage for previous halfedge (halfedges). The bit is set by default in the DefaultTraits.
		Status = 8,  ///< Add status to mesh item (all items)
		TexCoord1D = 16, ///< Add 1D texture coordinates (vertices, halfedges)
		TexCoord2D = 32, ///< Add 2D texture coordinates (vertices, halfedges)
		TexCoord3D = 64, ///< Add 3D texture coordinates (vertices, halfedges)
		TexCoordAll = TexCoord1D | TexCoord2D | TexCoord3D,
		TextureIndex = 128 ///< Add texture index (faces)
	};

	inline bool operator!(AttributeBits x) {
		return x == AttributeBits::None;
	}

	inline AttributeBits operator~(AttributeBits x) {
		return AttributeBits(~unsigned(x));
	}

	inline AttributeBits operator|(AttributeBits l, AttributeBits r) {
		return AttributeBits(unsigned(l) | unsigned(r));
	}

	inline AttributeBits operator&(AttributeBits l, AttributeBits r) {
		return AttributeBits(unsigned(l) & unsigned(r));
	}

	inline AttributeBits &operator|=(AttributeBits &l, AttributeBits r) {
		l = l | r;
		return l;
	}

	inline AttributeBits &operator&=(AttributeBits &l, AttributeBits r) {
		l = l & r;
		return l;
	}

	inline std::string to_string(const AttributeBits &attrs)
	{
		std::string r;
		if (!!(attrs & AttributeBits::Normal)) r += "Normal|";
		if (!!(attrs & AttributeBits::Color)) r += "Color|";
		if (!!(attrs & AttributeBits::PrevHalfedge)) r += "PrevHalfedge|";
		if (!!(attrs & AttributeBits::Status)) r += "Status|";
		if (!!(attrs & AttributeBits::TexCoord1D)) r += "TexCoord1D|";
		if (!!(attrs & AttributeBits::TexCoord2D)) r += "TexCoord2D|";
		if (!!(attrs & AttributeBits::TexCoord3D)) r += "TexCoord3D|";
		if (!!(attrs & AttributeBits::TextureIndex)) r += "TextureIndex|";
		if (r.length()) r.pop_back();
		return r;
	}

	//=============================================================================

} // namespace OpenMesh
//=============================================================================
#endif // OPENMESH_ATTRIBUTES_HH defined
//=============================================================================
