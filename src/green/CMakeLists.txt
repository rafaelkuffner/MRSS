
# for openmp
cmake_minimum_required(VERSION 3.9)

set(sources
	CMakeLists.txt
	"${PROJECT_SOURCE_DIR}/cmake/gitver.cmake"
	gitverdetail.h.in
	gitver.hpp
	gitver.cpp
	main.hpp
	main.cpp
	config.hpp
	config.cpp
	model.hpp
	model.cpp
	model_entity.hpp
	model_entity.cpp
	assimp_openmesh.hpp
	assimp_openmesh.cpp
	entity.hpp
	entity.cpp
	simd.h
	decimate.hpp
	decimate.cpp
	saliency.hpp
	saliency.cpp
	neighborhood.hpp
	neighborhood.cpp
	curvature.hpp
	curvature.cpp
	meshutils.hpp
	meshutils.cpp
	dialog.hpp
	dialog_win32.cpp
	dialog_simple.cpp
	imguiex.hpp
	imguiex.cpp
	clipp.h
	uilocale.hpp
	uilocale.cpp
	uilocale_en.hpp
	uilocale_jp.hpp
	model.glsl
	deferred.glsl
	about_licences.txt
)

cgu_create_missing(${sources})

if(WIN32)
	list(APPEND sources
		green.rc
		green.ico
	)
endif()

add_executable(green ${sources})
target_compile_features(green PUBLIC cxx_std_17)
target_link_libraries(green PRIVATE cgu fmt assimp OpenMeshCore OpenMeshTools)
cgu_target_default_build_options(green)
cgu_target_source_group_tree(green)

if(WIN32)
	target_compile_definitions(green PRIVATE GREEN_DIALOG_WIN32)
else()
	target_compile_definitions(green PRIVATE GREEN_DIALOG_SIMPLE)
endif()

# shaders
cgu_target_embed_string(green model.glsl glsl_green_model)
cgu_target_embed_string(green deferred.glsl glsl_green_deferred)

# other strings
cgu_target_embed_string(green about_licences.txt about_licences)

# for openmesh
target_compile_definitions(green PUBLIC _USE_MATH_DEFINES)

cgu_target_msvc_compile_options(green PRIVATE "/arch:AVX2")
cgu_target_gcc_compile_options(green PRIVATE "-mavx2")

# msvc openmp conflicts with two-phase name lookup support
cgu_target_msvc_compile_options(green PRIVATE "/Zc:twoPhase-")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	target_link_libraries(green PRIVATE OpenMP::OpenMP_CXX)
endif()

set(gitversrc "${CMAKE_CURRENT_SOURCE_DIR}/gitverdetail.h.in")
set(gitverdst "${CMAKE_CURRENT_BINARY_DIR}/gitverdetail.h")

add_custom_target(
	green_gitver
	"${CMAKE_COMMAND}" "-Dsrc=${gitversrc}" "-Ddst=${gitverdst}" -P "${PROJECT_SOURCE_DIR}/cmake/gitver.cmake"
	DEPENDS "${gitversrc}"
	BYPRODUCTS "${gitverdst}"
	WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
	VERBATIM
)

add_dependencies(green green_gitver)
