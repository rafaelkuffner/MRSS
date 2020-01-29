
set(sources
	CMakeLists.txt
	"${PROJECT_SOURCE_DIR}/cmake/cgu.cmake"
	opengl.hpp
	opengl.cpp
	util.hpp
	mesh.hpp
	mesh.cpp
	shader.hpp
	shader.cpp
)

cgu_create_missing(${sources})

add_library(cgu ${sources})
target_compile_features(cgu PUBLIC cxx_std_17)
target_link_libraries(cgu PUBLIC glew glfw glm imgui) # stb
cgu_target_default_build_options(cgu)
cgu_target_source_group_tree(cgu)

# glm uses unnamed structs without turning off the warning about non-standardness
cgu_target_msvc_compile_options(cgu PUBLIC /wd4201)
