
#
# BSD 3-Clause License
# 
# Copyright (c) 2013-2019, Benjamin Allen
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

# require new behaviour of:
# CMP0054 (don't dereference quoted variables in if() args)
cmake_minimum_required(VERSION 3.1)

option(CGU_PDB_FAST "CGU build options: use MSVC /debug:fastlink" ON)
option(CGU_LTCG "CGU build options: use link-time codegen" ON)

# Set CMake output dirs to "${CMAKE_BINARY_DIR}/bin" if the current project is the root project
function(cgu_maybe_set_output_dirs)
	# necessary for building shared libs so they all go in the same place and can then be loaded
	# however, the current value takes effect when a target is created, so overriding as a subproject would be bad
	# therefore, we only set the output directories if we are the top-level project
	if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
		message(STATUS "CGU: setting CMake output directories to '${CMAKE_BINARY_DIR}/bin'")
		set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
		set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
		set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	endif()
endfunction()

# Set CGU default build options (warnings etc.)
function(cgu_target_default_build_options target)
	# TODO gcc/clang optimizaton flags?
	# TODO floating-point behaviour (msvc /fp:fast can break things)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
		# note: VS2017 required for C++17
		# Multiprocess build
		target_compile_options(${target} PRIVATE /MP)
		# Disable non-standard behaviour
		target_compile_options(${target} PRIVATE /permissive-)
		# UTF-8 source and execution charsets
		target_compile_options(${target} PRIVATE /utf-8)
		# Full normal warnings
		target_compile_options(${target} PRIVATE /W4)
		# Disable C4800: forcing X to bool (performance warning)
		target_compile_options(${target} PRIVATE /wd4800)
		# Not debug: enable intrinsics
		target_compile_options(${target} PRIVATE "$<$<NOT:$<CONFIG:Debug>>:/Oi>")
		# Function-level linking
		target_compile_options(${target} PRIVATE /Gy)
		# PDB files
		target_compile_options(${target} PRIVATE /Zi)
		target_link_libraries(${target} PRIVATE "$<$<BOOL:${CGU_PDB_FAST}>:-debug:fastlink>" "$<$<NOT:$<BOOL:${CGU_PDB_FAST}>>:-debug:full>")
		if(CGU_LTCG)
			# Link-time codegen
			target_link_libraries(${target} PRIVATE -ltcg:incremental)
			# Whole program optimization
			target_compile_options(${target} PRIVATE /GL)
			# LTCG also has to be specified separately for static libs
			set_property(TARGET ${target} APPEND PROPERTY STATIC_LIBRARY_FLAGS /ltcg)
		endif()
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		# Full normal warnings
		target_compile_options(${target} PRIVATE -Wall -Wextra -pedantic)
		# Threading support
		target_compile_options(${target} PRIVATE -pthread)
		# Promote missing return to error
		target_compile_options(${target} PRIVATE -Werror=return-type)
		# enable coloured output if gcc >= 4.9
		# TODO remove this, use the GCC_COLORS env instead?
		execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
		if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
			target_compile_options(${target} PRIVATE -fdiagnostics-color)
		endif()
	elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
		# Full normal warnings
		target_compile_options(${target} PRIVATE -Wall -Wextra -pedantic)
		# Threading support
		target_compile_options(${target} PRIVATE -pthread)
		# Promote missing return to error
		target_compile_options(${target} PRIVATE -Werror=return-type)
	endif()
endfunction()

# Add MSVC-style compile options
function(cgu_target_msvc_compile_options)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
		target_compile_options(${ARGV})
	endif()
endfunction()

# Add GCC-style compile options
function(cgu_target_gcc_compile_options)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(${ARGV})
	elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
		target_compile_options(${ARGV})
	endif()
endfunction()

# Set source groups to mirror real directory tree
function(cgu_target_source_group_tree target)
	get_target_property(sources ${target} SOURCES)
	foreach(source ${sources})
		get_filename_component(source ${source} ABSOLUTE)
		string(REPLACE "${PROJECT_SOURCE_DIR}/" "" rel ${source})
		if(rel)
			string(REGEX REPLACE "/([^/]*)$" "" rel ${rel})
			if(NOT rel STREQUAL source)
				string(REPLACE "/" "\\\\" rel ${rel})
				source_group(${rel} FILES ${source})
			endif()
		endif()
	endforeach()
endfunction()

# Add sources to a target relative to current CMakeLists.txt directory (CMAKE_CURRENT_SOURCE_DIR)
function(cgu_target_relative_sources target)
	foreach(source ${ARGN})
		target_sources(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
	endforeach()
endfunction()

# Create missing source files relative to current CMakeLists.txt directory (CMAKE_CURRENT_SOURCE_DIR)
function(cgu_create_missing)
	foreach(source ${ARGN})
		if(NOT IS_ABSOLUTE ${source})
			set(path "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
			if (NOT EXISTS ${path})
				message(STATUS "cgu_create_missing: ${path}")
				file(APPEND ${path} "")
			endif()
		endif()
	endforeach()
endfunction()
