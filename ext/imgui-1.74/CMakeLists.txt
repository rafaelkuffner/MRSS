
cmake_minimum_required(VERSION 3.1)

project("imgui" CXX C)

add_library(imgui
	"include/imstb_truetype.h"
	"include/imconfig.h"
	"include/imgui.h"
	"include/imgui_internal.h"
	"include/imstb_rectpack.h"
	"include/imstb_textedit.h"
	"include/imgui_impl_glfw.h"
	"include/imgui_impl_opengl3.h"
	"source/imgui_widgets.cpp"
	"source/imgui.cpp"
	"source/imgui_demo.cpp"
	"source/imgui_draw.cpp"
	"source/imgui_impl_glfw.cpp"
	"source/imgui_impl_opengl3.cpp"
)

target_include_directories(imgui PUBLIC "include")

target_link_libraries(imgui PUBLIC glfw glew)
target_compile_definitions(imgui PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLEW)
