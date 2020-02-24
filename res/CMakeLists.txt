# add resources as targets so they show up in IDEs

set(sources
	CMakeLists.txt
	model.glsl
	deferred.glsl
)

cgu_create_missing(${sources})

add_custom_target(
	green_res
	SOURCES
	${sources}
)

cgu_target_source_group_tree(green_res)
