
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

#ifdef _VERTEX_

layout(location=0) in vec3 a_pos_m;
layout(location=1) in vec3 a_norm_m;

out VertexData {
	vec3 pos_v;
	vec3 norm_v;
	flat int id;
} v_out;

void main() {
	vec4 pos_v = u_modelview * vec4(a_pos_m + u_pos_bias, 1);
	gl_Position = u_projection * pos_v;
	v_out.pos_v = pos_v.xyz;
	v_out.norm_v = normalize(vec3(u_modelview * vec4(a_norm_m + vec3(0,0.0001,0), 0)));
	v_out.id = gl_VertexID;
}

#endif

#ifdef _FRAGMENT_

in VertexData {
	vec3 pos_v;
	vec3 norm_v;
	flat int id;
} v_in;

layout(location=0) out vec4 f_color;
layout(location=1) out ivec2 f_id;

bool hovered() {
	return any(equal(u_selection.xy, ivec2(u_entity_id, v_in.id)));
}

bool selected() {
	return any(equal(u_selection.zw, ivec2(u_entity_id, v_in.id)));
}

void main() {
	vec3 norm_v = normalize(v_in.norm_v);
	vec4 color = hovered() ? u_color_hover : selected() ? u_color_select : u_color;
	f_color = vec4(color.xyz * (abs(norm_v.z) * u_shading + 1.0 - u_shading), 1);
	gl_FragDepth = frag_depth() + u_depth_bias;
	f_id.x = u_entity_id;
	f_id.y = v_in.id;
}

#endif
