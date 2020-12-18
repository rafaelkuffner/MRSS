
//! #version 330 core

uniform mat4 u_modelview;
uniform mat4 u_projection;
uniform vec4 u_color;
uniform vec3 u_pos_bias = vec3(0);
uniform float u_shading = 0;
uniform float u_depth_bias = 0;
uniform int u_entity_id;
uniform int u_vert_color_map = 0;
uniform bool u_show_samples = false;

#ifdef _VERTEX_

layout(location=0) in vec3 a_pos_m;
layout(location=1) in vec3 a_norm_m;

// saliency in R or color in RGB; sampled flag in A
layout(location=2) in vec4 a_color;

out VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
	// this causes problems with intel igpus, and we're not using it atm anyway
	//flat bool should_discard;
} v_out;

vec4 gamma_decode(vec4 c) {
	return vec4(pow(c.rgb, vec3(2.2)), c.a);
}

vec4 map_color_zbrush(float k) {
	k = clamp(k, 0.0, 1.0);
	// high color: red
	vec4 c0 = vec4(1, 0, 0, 1);
	// cyan to white
	vec4 cl = vec4(k * 2, 1, 1, 1);
	// white to <high>
	vec4 ch = mix(vec4(1), c0, 2 * k - 1);
	// we have to interpolate in gamma-encoded space for it to look 'correct'
	return gamma_decode(mix(cl, ch, bvec4(k > 0.5)));
}

vec4 map_color(vec4 c) {
	switch (u_vert_color_map) {
	case 0:
		// uniform color
		return u_color;
	case 1:
		// attribute color
		return c;
	case 2:
		// scalar attribute * uniform color
		return c.r * u_color;
	case 3:
		// zbrush style
		return map_color_zbrush(c.r);
	case 4:
		// difference, zbrush style
		return map_color_zbrush(c.r * 0.5 + 0.5);
	default:
		return c;
	}
}

void main() {
	vec4 pos_v = u_modelview * vec4(a_pos_m + u_pos_bias, 1);
	gl_Position = u_projection * pos_v;
	v_out.pos_v = pos_v.xyz;
#ifdef SHADE_SMOOTH
	v_out.norm_v = normalize(vec3(u_modelview * vec4(a_norm_m + vec3(0,0.0001,0), 0)));
#endif
	v_out.color = map_color(a_color);
	v_out.id = gl_VertexID;
	//v_out.should_discard = u_show_samples && a_color.a < 0.5;
}

#endif

#ifdef _GEOMETRY_

layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;

in VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
	//flat bool should_discard;
} v_in[];

out VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
	//flat bool should_discard;
} v_out;

void main() {
#ifdef SHADE_FLAT
	// calc triangle face normal
	vec3 fnorm_v = normalize(cross(v_in[1].pos_v - v_in[0].pos_v, v_in[2].pos_v - v_in[1].pos_v));
#endif
	for (int i = 0; i < 3; i++) {
		gl_Position = gl_in[i].gl_Position;
		v_out.pos_v = v_in[i].pos_v;
#ifdef SHADE_FLAT
		v_out.norm_v = fnorm_v;
#else
		v_out.norm_v = v_in[i].norm_v;
#endif
		v_out.color = v_in[i].color;
		v_out.id = v_in[i].id;
		EmitVertex();
	}
	EndPrimitive();
}

#endif

#ifdef _FRAGMENT_

in VertexData {
	vec3 pos_v;
	vec3 norm_v;
	vec4 color;
	flat int id;
	//flat bool should_discard;
} v_in;

layout(location=0) out vec4 f_color;
layout(location=1) out ivec2 f_id;

void main() {
	//if (v_in.should_discard) discard;
	vec3 norm_v = normalize(v_in.norm_v);
	f_color = vec4(v_in.color.rgb * (abs(norm_v.z) * u_shading + 1.0 - u_shading), 1);
	gl_FragDepth = frag_depth() + u_depth_bias;
	f_id.x = u_entity_id;
	f_id.y = v_in.id;
}

#endif
