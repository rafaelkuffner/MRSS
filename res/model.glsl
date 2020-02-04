
//! #version 330 core

uniform mat4 u_modelview;
uniform mat4 u_projection;
uniform vec4 u_color;
uniform vec3 u_pos_bias = vec3(0);
uniform float u_shading = 0;
uniform float u_depth_bias = 0;

#ifdef _VERTEX_

layout(location=0) in vec3 a_pos_m;
layout(location=1) in vec3 a_norm_m;

out VertexData {
	vec3 pos_v;
	vec3 norm_v;
} v_out;

void main() {
	vec4 pos_v = u_modelview * vec4(a_pos_m + u_pos_bias, 1);
	gl_Position = u_projection * pos_v;
	v_out.pos_v = pos_v.xyz;
	v_out.norm_v = normalize(vec3(u_modelview * vec4(a_norm_m + vec3(0,0.0001,0), 0)));
}

#endif

#ifdef _FRAGMENT_

in VertexData {
	vec3 pos_v;
	vec3 norm_v;
} v_in;

layout(location=0) out vec4 f_color;

void main() {
	vec3 norm_v = normalize(v_in.norm_v);
	f_color = vec4(u_color.xyz * (abs(norm_v.z) * u_shading + 1.0 - u_shading), 1);
	gl_FragDepth = frag_depth() + u_depth_bias;
}

#endif
