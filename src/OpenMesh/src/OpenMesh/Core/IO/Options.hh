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




#ifndef OPENMESH_IO_OPTIONS_HH
#define OPENMESH_IO_OPTIONS_HH


 //=== INCLUDES ================================================================


 // OpenMesh
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Attributes.hh>

#include <cstdlib>

//== NAMESPACES ==============================================================


namespace OpenMesh {
	namespace IO {

		//=== IMPLEMENTATION ==========================================================

		/// Definitions of %Options for reading and writing. The options can be
		/// or'ed.
		enum class OptionBits : unsigned {
			None = 0,
			Default = 0x00, ///< No options
			Binary = 0x01, ///< Set binary mode for r/w
			MSB = 0x02, ///< Assume big endian byte ordering
			LSB = 0x04, ///< Assume little endian byte ordering
			Swap = 0x08, ///< Swap byte order in binary mode
			// TODO per-element color settings? (-> attribs)
			// - are these related to the Color typedef in the mesh traits?
			ColorAlpha = 0x10, ///< Has (r) / store (w) alpha values for colors
			ColorFloat = 0x20, ///< Has (r) / store (w) float values for colors (currently only implemented for PLY and OFF files)
			Custom = 0x40, ///< Has (r)             custom properties (currently only implemented in PLY Reader ASCII version)
		};

		inline bool operator!(OptionBits x) {
			return x == OptionBits::None;
		}

		inline OptionBits operator~(OptionBits x) {
			return OptionBits(~unsigned(x));
		}

		inline OptionBits operator|(OptionBits l, OptionBits r) {
			return OptionBits(unsigned(l) | unsigned(r));
		}

		inline OptionBits operator&(OptionBits l, OptionBits r) {
			return OptionBits(unsigned(l) & unsigned(r));
		}

		inline OptionBits & operator|=(OptionBits &l, OptionBits r) {
			l = l | r;
			return l;
		}

		inline OptionBits & operator&=(OptionBits &l, OptionBits r) {
			l = l & r;
			return l;
		}

		/** \name Mesh Reading / Writing
			Option for reader and writer modules.
		*/

		//-----------------------------------------------------------------------------

		/** \brief Set options for reader/writer modules.
		 *
		 *  The class is used in a twofold way.
		 *  -# In combination with reader modules the class is used
		 *     - to pass hints to the reading module, whether the input is
		 *     binary and what byte ordering the binary data has
		 *     - to retrieve information about the file contents after
		 *     succesful reading.
		 *  -# In combination with write modules the class gives directions to
		 *     the writer module, whether to
		 *     - use binary mode or not and what byte order to use
		 *     - store one of the standard properties.
		 *
		 *  The option are defined in \c Options::Flag as bit values and stored in
		 *  an \c int value as a bitset.
		 */
		struct Options
		{
			OptionBits flags = OptionBits::Default;
			AttributeBits vattribs = AttributeBits::None;
			AttributeBits eattribs = AttributeBits::None;
			AttributeBits hattribs = AttributeBits::None;
			AttributeBits fattribs = AttributeBits::None;

			/// Restore state after default constructor.
			void cleanup(void)
			{
				*this = Options();
			}

			/// Clear all bits.
			void clear(void)
			{
				*this = Options();
			}

			/// Returns true if all bits are zero.
			bool is_empty(void) const
			{
				return *this == Options();
			}

			/// Unset options defined in _rhs.
			Options &operator-=(const OptionBits _rhs)
			{
				flags = OptionBits(flags & ~_rhs);
				return *this;
			}

			Options &unset(const OptionBits _rhs)
			{
				return (*this -= _rhs);
			}

			/// Set options defined in _rhs
			Options &operator+=(const OptionBits _rhs)
			{
				flags = OptionBits(flags | _rhs);
				return *this;
			}

			Options &set(const OptionBits _rhs)
			{
				return (*this += _rhs);
			}

			/// Returns true if _rhs has the same options enabled.
			bool operator==(const Options &rhs) const
			{
				// assuming bitwise equality (beware padding!)
				return std::memcmp(this, &rhs, sizeof(Options)) == 0;
			}

			/// Returns true if _rhs does not have the same options enabled.
			bool operator!=(const Options &rhs) const
			{
				return !(*this == rhs);
			}

			// Check if an option or several options are set.
			bool check(const OptionBits _rhs) const
			{
				return (flags & _rhs) == _rhs;
			}

			bool is_binary()               const { return check(OptionBits::Binary); }
			bool color_has_alpha()         const { return check(OptionBits::ColorAlpha); }
			bool color_is_float()          const { return check(OptionBits::ColorFloat); }
			bool vertex_has_normal()       const { return !!(vattribs & AttributeBits::Normal); }
			bool vertex_has_color()        const { return !!(vattribs & AttributeBits::Color); }
			bool vertex_has_texcoord1D()   const { return !!(vattribs & AttributeBits::TexCoord1D); }
			bool vertex_has_texcoord2D()   const { return !!(vattribs & AttributeBits::TexCoord2D); }
			bool vertex_has_texcoord3D()   const { return !!(vattribs & AttributeBits::TexCoord3D); }
			bool vertex_has_texcoord()     const { return !!(vattribs & (AttributeBits::TexCoord1D | AttributeBits::TexCoord2D | AttributeBits::TexCoord3D)); }
			bool vertex_has_status()       const { return !!(vattribs & AttributeBits::Status); }
			bool halfedge_has_normal()     const { return !!(hattribs & AttributeBits::Normal); }
			bool halfedge_has_color()      const { return !!(hattribs & AttributeBits::Color); }
			bool halfedge_has_texcoord1D() const { return !!(hattribs & AttributeBits::TexCoord1D); }
			bool halfedge_has_texcoord2D() const { return !!(hattribs & AttributeBits::TexCoord2D); }
			bool halfedge_has_texcoord3D() const { return !!(hattribs & AttributeBits::TexCoord3D); }
			bool halfedge_has_texcoord()   const { return !!(hattribs & (AttributeBits::TexCoord1D | AttributeBits::TexCoord2D | AttributeBits::TexCoord3D)); }
			bool halfedge_has_status()     const { return !!(hattribs & AttributeBits::Status); }
			bool edge_has_color()          const { return !!(eattribs & AttributeBits::Color); }
			bool edge_has_status()         const { return !!(eattribs & AttributeBits::Status); }
			bool face_has_normal()         const { return !!(fattribs & AttributeBits::Normal); }
			bool face_has_color()          const { return !!(fattribs & AttributeBits::Color); }
			bool face_has_texindex()       const { return !!(fattribs & AttributeBits::TextureIndex); }
			bool face_has_status()         const { return !!(fattribs & AttributeBits::Status); }
		};

		//=============================================================================

	} // namespace IO
} // namespace OpenMesh
//=============================================================================
#endif
//=============================================================================
