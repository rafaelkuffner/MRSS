
# require new behaviour of:
# CMP0054 (don't dereference quoted variables in if() args)
cmake_minimum_required(VERSION 3.1)

project("GREEN" CXX C)

# rename the folder for predefined targets to always be at the top of the solution
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".CMakePredefinedTargets")

# enable IDE project folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# utilities
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(cgu)

# output directories
cgu_maybe_set_output_dirs()

# AssImp defaults this to on; we want it off
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# GLEW
add_subdirectory("ext/glew-2.1.0")
set_property(TARGET glew PROPERTY FOLDER "ext")

# GLFW
option(GLFW_BUILD_EXAMPLES "Build the GLFW example programs" OFF)
option(GLFW_BUILD_TESTS "Build the GLFW test programs" OFF)
option(GLFW_BUILD_DOCS "Build the GLFW documentation" OFF)
option(GLFW_INSTALL "Generate installation target" OFF)
add_subdirectory("ext/glfw-3.3")
set_property(TARGET glfw PROPERTY FOLDER "ext")

# GLM
add_subdirectory("ext/glm-0.9.9.5/glm")
set_property(TARGET glm PROPERTY FOLDER "ext")

# ImGui
add_subdirectory("ext/imgui-1.74")
set_property(TARGET imgui PROPERTY FOLDER "ext")

# fmtlib
add_subdirectory("ext/fmt-7.1.3")
set_property(TARGET fmt PROPERTY FOLDER "ext")

#add_subdirectory("ext/stb")
#set_property(TARGET stb PROPERTY FOLDER "ext")

# AssImp
option(ASSIMP_BUILD_ASSIMP_TOOLS "If the supplementary tools for Assimp are built in addition to the library." OFF)
option(ASSIMP_BUILD_TESTS "If the test suite for Assimp is built in addition to the library." OFF)
option(ASSIMP_BUILD_ZLIB "Build your own zlib" ON)
add_subdirectory("ext/assimp-5.0.1")
set_property(TARGET assimp PROPERTY FOLDER "ext")
set_property(TARGET IrrXML PROPERTY FOLDER "ext")
set_property(TARGET zlib PROPERTY FOLDER "ext")
set_property(TARGET zlibstatic PROPERTY FOLDER "ext")
if(TARGET UpdateAssimpLibsDebugSymbolsAndDLLs)
	set_property(TARGET UpdateAssimpLibsDebugSymbolsAndDLLs PROPERTY FOLDER "ext")
endif()

# allow including headers from src tree with <>
include_directories("src")

# io lib
add_subdirectory("src/iobuffer")

# OpenMesh (modified)
add_subdirectory("src/OpenMesh")
include_directories("src/OpenMesh/src")

# cgu lib
add_subdirectory("src/cgu")

# main exe
add_subdirectory("src/green")

# resources
add_subdirectory("res")
