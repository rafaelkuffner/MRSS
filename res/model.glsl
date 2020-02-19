
//! #version 330 core

uniform mat4 u_modelview;
uniform mat4 u_projection;
uniform vec4 u_color;
uniform vec4 u_color_hover;
uniform vec4 u_color_select;
uniform vec3 u_pos_bias = vec3(0);
uniform float u_shading = 0;
uniform float u_depth_bias = 0;
uniform int u_entity_id;
uniform ivec4 u_selection;
uniform bool u_use_vert_color = false;
uniform int u_vert_color_map = 0;

bool vertex_hovered(int vid) {
	return any(equal(u_selection.xy, ivec2(u_entity_id, vid)));
}

bool vertex_selected(int vid) {
	return any(equal(u_selection.zw, ivec2(u_entity_id, vid)));
}

#ifdef _VERTEX_

layout(location=0) in vec3 a_pos_m;
layout(location=1) in vec3 a_norm_m;
layout(location=2) in vec4 a_color;

out VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
} v_out;

vec4 gamma_decode(vec4 c) {
	return vec4(pow(c.rgb, vec3(2.2)), c.a);
}

vec4 map_color_zbrush(float k, bool h, bool s) {
	// high color: red or hover/select
	vec4 c0 = (h ? u_color_hover * 1.5 : s ? u_color_select * 1.5 : vec4(1, 0, 0, 1));
	// cyan to white
	vec4 cl = vec4(k * 2, 1, 1, 1);
	// white to <high>
	vec4 ch = mix(vec4(1), c0, 2 * k - 1);
	// we have to interpolate in gamma-encoded space for it to look 'correct'
	return gamma_decode(mix(cl, ch, bvec4(k > 0.5)));
}

vec4 map_color(vec4 c, bool h, bool s) {
	switch (u_vert_color_map) {
	case 0:
		// direct color
		return h ? u_color_hover : s ? u_color_select : c;
	case 1:
		// scalar attribute * uniform color
		return c.r * (h ? u_color_hover : s ? u_color_select : u_color);
	case 2:
		// zbrush style (with different colors for hover/select)
		return map_color_zbrush(c.r, h, s);
	default:
		return c;
	}
}

void main() {
	vec4 pos_v = u_modelview * vec4(a_pos_m + u_pos_bias, 1);
	gl_Position = u_projection * pos_v;
	v_out.pos_v = pos_v.xyz;
	v_out.norm_v = normalize(vec3(u_modelview * vec4(a_norm_m + vec3(0,0.0001,0), 0)));
	v_out.color = map_color(u_use_vert_color ? a_color : u_color, vertex_hovered(gl_VertexID), vertex_selected(gl_VertexID));
	v_out.id = gl_VertexID;
}

#endif

#ifdef _FRAGMENT_

in VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
} v_in;

layout(location=0) out vec4 f_color;
layout(location=1) out ivec2 f_id;

void main() {
	vec3 norm_v = normalize(v_in.norm_v);
	f_color = vec4(v_in.color.rgb * (abs(norm_v.z) * u_shading + 1.0 - u_shading), 1);
	gl_FragDepth = frag_depth() + u_depth_bias;
	f_id.x = u_entity_id;
	f_id.y = v_in.id;
}

#endif
