
# for openmp
cmake_minimum_required(VERSION 3.9)

set(sources
	CMakeLists.txt
	main.cpp
	model.hpp
	model.cpp
	assimp_openmesh.hpp
	assimp_openmesh.cpp
	entity.hpp
	entity.cpp
	simd.h
	plf_colony.h
	saliency.hpp
	saliency.cpp
	neighborhood.hpp
	neighborhood.cpp
	curvature.hpp
	curvature.cpp
	meshutils.hpp
	meshutils.cpp
	Histogram.h
	Histogram.cpp
	FloatHistogram.h
	FloatHistogram.cpp
)

cgu_create_missing(${sources})

add_executable(green ${sources})
target_compile_features(green PUBLIC cxx_std_17)
target_link_libraries(green PRIVATE cgu assimp OpenMeshCore OpenMeshTools)
cgu_target_default_build_options(green)
cgu_target_source_group_tree(green)

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
