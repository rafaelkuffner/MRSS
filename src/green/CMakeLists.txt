
set(sources
	CMakeLists.txt
	main.cpp
)

cgu_create_missing(${sources})

add_executable(green ${sources})
target_compile_features(green PUBLIC cxx_std_17)
target_link_libraries(green PRIVATE cgu)
cgu_target_default_build_options(green)
cgu_target_source_group_tree(green)