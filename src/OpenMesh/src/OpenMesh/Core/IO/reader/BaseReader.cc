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




//=== INCLUDES ================================================================

#include <OpenMesh/Core/IO/reader/BaseReader.hh>

#if defined(OM_CC_MIPS)
#  include <ctype.h>
#else
#endif


//== NAMESPACES ===============================================================


namespace OpenMesh {
namespace IO {


//=== IMPLEMENTATION ==========================================================


static inline char tolower(char c) 
{
  using namespace std;
  return ::tolower(c); 
}


//-----------------------------------------------------------------------------


bool 
BaseReader::
can_u_read(const std::filesystem::path& _filename) const 
{
  // get file extension
  auto extension = _filename.extension().u8string();
  
  if (!extension.empty())
    extension = extension.substr(1);
  else
    extension = _filename.u8string(); //check, if the whole filename defines the extension

  // NOT UNICODE SAFE
  std::transform( extension.begin(), extension.end(),
	  extension.begin(), tolower );

  // locate extension in extension string
  return (get_extensions().find(extension) != std::string::npos);
}


//-----------------------------------------------------------------------------


bool 
BaseReader::
check_extension(const std::filesystem::path& _fname, const std::filesystem::path& _ext) const
{
  auto cmpExt = _ext.u8string();

  // NOT UNICODE SAFE
  std::transform( cmpExt.begin(), cmpExt.end(),  cmpExt.begin(), tolower );

  auto ext = _fname.extension().u8string();

  if (!ext.empty())
  { 
    // extension without dot!
    ext = ext.substr(1);

    // NOT UNICODE SAFE
    std::transform( ext.begin(), ext.end(), ext.begin(), tolower );
    
    return ext == cmpExt;
  }
  return false;  
}


//=============================================================================
} // namespace IO
} // namespace OpenMesh
//=============================================================================
