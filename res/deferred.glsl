
//! #version 330 core

#ifdef _FRAGMENT_

uniform sampler2D u_sampler_color;
uniform sampler2D u_sampler_depth;
uniform isampler2D u_sampler_id;

// hover entity, hover vertex, select entity, select vertex
uniform ivec4 u_selection;

// these colors are linear, not srgb
uniform vec4 u_hover_entity_color = vec4(0.8, 0.8, 0, 1);
uniform vec4 u_hover_vertex_color = vec4(1, 0.3, 0.3, 1);
uniform vec4 u_select_entity_color = vec4(1, 0.3, 0, 1);
uniform vec4 u_select_vertex_color = vec4(0, 1, 0, 1);

void main() {
	gl_FragDepth = texture(u_sampler_depth, v_tex_coord).r;
	vec2 dtx = vec2(dFdx(v_tex_coord.x), dFdy(v_tex_coord.y));
	const int searchpx = 4;
	vec4 dist = vec4(9001);
	bvec4 selmask = greaterThanEqual(u_selection, ivec4(0));
	if (any(selmask)) {
		for (float x = -searchpx * dtx.x; x <= searchpx * dtx.x; x += dtx.x) {
			for (float y = -searchpx * dtx.y; y <= searchpx * dtx.y; y += dtx.y) {
				ivec2 id = texture(u_sampler_id, v_tex_coord + vec2(x, y)).xy;
				dist = mix(dist, min(dist, vec4(length(vec2(x, y)))), equal(id.xyxy, u_selection) && selmask);
			}
		}
	}
	vec4 color0 = texture(u_sampler_color, v_tex_coord);
	f_color = color0;
	f_color = mix(f_color, color0 * 0.5 + u_select_entity_color, bvec4(dist.z > 0 && dist.z < 1));
	f_color = mix(f_color, color0 * 0.5 + u_hover_entity_color, bvec4(dist.x > 0 && dist.x < 1));
	f_color = mix(f_color, color0 * 0.5 + u_select_vertex_color, bvec4(dist.w > 0 && dist.w < 1));
	f_color = mix(f_color, color0 * 0.5 + u_hover_vertex_color, bvec4(dist.y > 0 && dist.y < 1));
}

#endif
