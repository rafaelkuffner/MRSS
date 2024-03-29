
//! #version 330 core

#ifdef _FRAGMENT_

uniform sampler2D u_sampler_color;
uniform sampler2D u_sampler_depth;
uniform isampler2D u_sampler_id;

// hover entity, hover vertex, select entity, select vertex
uniform ivec4 u_selection;
uniform bvec4 u_show_selection;

// these colors are linear, not srgb
uniform vec4 u_hover_entity_color = vec4(0.8, 0.8, 0, 1);
uniform vec4 u_hover_vertex_color = vec4(0.8, 0.8, 0, 1);
//uniform vec4 u_hover_vertex_color = vec4(1, 0.3, 0.3, 1);
uniform vec4 u_select_entity_color = vec4(1, 0.4, 0, 1);
uniform vec4 u_select_vertex_color = vec4(0.2, 0.8, 0, 1);
//uniform vec4 u_select_vertex_color = vec4(0, 1, 0, 1);

// && doesn't operator on vectors, & only operates on integers (+vectors)
// AMD is strict on this, NVidia is more relaxed. the GLSL spec is silly.
bvec4 bvec_and(bvec4 a, bvec4 b) {
	return bvec4(a.x && b.x, a.y && b.y, a.z && b.z, a.w && b.w);
}

bvec4 bvec_and(bvec4 a, bvec4 b, bvec4 c) {
	return bvec_and(a, bvec_and(b, c));
}

void main() {
	gl_FragDepth = texture(u_sampler_depth, v_tex_coord).r;
	vec2 dtx = vec2(dFdx(v_tex_coord.x), dFdy(v_tex_coord.y));
	const int searchpx = 2;
	vec4 dist = vec4(9001);
	bvec4 selmask = greaterThanEqual(u_selection, ivec4(0));
	if (any(selmask)) {
		// TODO this is kinda slow
		for (float x = -searchpx * dtx.x; x <= searchpx * dtx.x; x += dtx.x) {
			for (float y = -searchpx * dtx.y; y <= searchpx * dtx.y; y += dtx.y) {
				ivec2 id = texture(u_sampler_id, v_tex_coord + vec2(x, y)).xy;
				dist = mix(dist, min(dist, vec4(length(vec2(x, y)))), bvec_and(equal(id.xyxy, u_selection), equal(id.xxxx, u_selection.xxzz), selmask));
			}
		}
	}
	vec4 color0 = texture(u_sampler_color, v_tex_coord);
	f_color = color0;
	f_color = mix(f_color, color0 * 0.5 + u_select_entity_color, bvec4(u_show_selection.z && dist.z > 0 && dist.z < 1));
	f_color = mix(f_color, color0 * 0.5 + u_hover_entity_color, bvec4(u_show_selection.x && dist.x > 0 && dist.x < 1));
	f_color = mix(f_color, color0 * 0.5 + u_select_vertex_color, bvec4(u_show_selection.w && dist.w >= 0 && dist.w < 1));
	f_color = mix(f_color, color0 * 0.5 + u_hover_vertex_color, bvec4(u_show_selection.y && dist.y >= 0 && dist.y < 1));
}

#endif
