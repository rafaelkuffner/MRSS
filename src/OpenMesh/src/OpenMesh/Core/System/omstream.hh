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
 //  OpenMesh streams: omlog, omout, omerr
 //
 //=============================================================================

#ifndef OPENMESH_OMSTREAMS_HH
#define OPENMESH_OMSTREAMS_HH


//== INCLUDES =================================================================

//#include <OpenMesh/Core/System/mostream.hh>
#include <OpenMesh/Core/System/config.h>

#include <string>
#include <string_view>
#include <sstream>
#include <type_traits>
#include <utility>

//== CLASS DEFINITION =========================================================

namespace OpenMesh
{
	enum class log_level {
		debug, info, warning, error
	};

	struct log_message {
		std::string_view source;
		std::string_view file;
		int line;
		log_level level;
		std::string message;
	};

	using log_handler_t = void(*)(log_message &&) noexcept;

	OPENMESHDLLEXPORT void set_log_handler(log_handler_t) noexcept;

	OPENMESHDLLEXPORT log_handler_t get_log_handler() noexcept;

	OPENMESHDLLEXPORT void default_log_handler(log_message &&) noexcept;

	class log_submission {
	private:
		log_message m_msg;
		std::ostringstream m_oss;

	public:
		~log_submission()
		{
			m_msg.message = m_oss.str();
			get_log_handler()(std::move(m_msg));
		}

		log_submission(std::string_view source_, std::string_view file_, int line_, log_level level_) :
			m_msg{source_, file_, line_, level_}
		{}

		log_submission(const log_submission &) = delete;
		log_submission & operator=(const log_submission &) = delete;

		template <typename T, typename = std::void_t<decltype(std::declval<std::ostream &>() << std::declval<const T &>())>>
		log_submission & operator<<(const T &t)
		{
			m_oss << t;
			return *this;
		}

	};

}

#define OM_STRINGIZE_IMPL(x) #x
#define OM_STRINGIZE(x) OM_STRINGIZE_IMPL(x)

// #define OMLOG_SOURCE to something before using these
#define OMLOG_DEBUG ::OpenMesh::log_submission(OM_STRINGIZE(OMLOG_SOURCE), __FILE__, __LINE__, ::OpenMesh::log_level::debug)
#define OMLOG_INFO ::OpenMesh::log_submission(OM_STRINGIZE(OMLOG_SOURCE), __FILE__, __LINE__, ::OpenMesh::log_level::info)
#define OMLOG_WARNING ::OpenMesh::log_submission(OM_STRINGIZE(OMLOG_SOURCE), __FILE__, __LINE__, ::OpenMesh::log_level::warning)
#define OMLOG_ERROR ::OpenMesh::log_submission(OM_STRINGIZE(OMLOG_SOURCE), __FILE__, __LINE__, ::OpenMesh::log_level::error)


//=============================================================================
#endif // OPENMESH_OMSTREAMS_HH defined
//=============================================================================
