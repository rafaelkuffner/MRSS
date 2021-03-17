
set(sources
	CMakeLists.txt
	iobuffer.hpp
	iobuffer.cpp
)

add_library(iobuffer ${sources})
target_compile_features(iobuffer PUBLIC cxx_std_17)
cgu_target_default_build_options(iobuffer)
cgu_target_source_group_tree(iobuffer)
