
/*
BSD 3-Clause License

Copyright (c) 2013-2019, Benjamin Allen and Joshua Scott
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cstddef>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

#include "mesh.hpp"

namespace cgu {
	
	void gl_mesh_slice::draw() const {
		assert(indexcount >= 0);
		GLint old_poly_mode;
		glBindVertexArray(vao);
		glDrawElements(mode, indexcount, GL_UNSIGNED_INT, (void *)(indexstart * sizeof(int)));
		glBindVertexArray(0);
	}

	gl_mesh mesh_builder_base::build_impl(const void *data, size_t size, size_t stride, gl_mesh m) const {

		// TODO do i need to MAP_INVALIDATE_BUFFER_BIT when streaming? (bufhint)

		bool makevao = false;

		if (!m.vao) makevao = true, m.vao = gl_object::gen_vertex_array();
		if (!m.vbo) m.vbo = gl_object::gen_buffer();
		if (!m.ibo) m.ibo = gl_object::gen_buffer();

		glBindVertexArray(m.vao);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), bufhint);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it

		glBindBuffer(GL_ARRAY_BUFFER, m.vbo);
		glBufferData(GL_ARRAY_BUFFER, size, data, bufhint);

		if (makevao) {
			for (auto &attrib : attribs) {
				int i = &attrib - attribs.data();
				switch (attrib.atttype) {
				case vertex_attrib_type::float_unnormalized:
					glEnableVertexAttribArray(i);
					glVertexAttribPointer(i, attrib.components, attrib.srctype, false, stride, (void *) attrib.offset);
					break;
				case vertex_attrib_type::float_normalized:
					glEnableVertexAttribArray(i);
					glVertexAttribPointer(i, attrib.components, attrib.srctype, true, stride, (void *) attrib.offset);
					break;
				case vertex_attrib_type::integer:
					glEnableVertexAttribArray(i);
					glVertexAttribIPointer(i, attrib.components, attrib.srctype, stride, (void *) attrib.offset);
					break;
				default:
					break;
				}
			}
		}

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		m.mode = mode;
		m.indexcount = int(indices.size());

		return m;
	}

	mesh_builder<> mesh_load_obj(const std::filesystem::path &filepath) {

		using namespace std;
		using namespace glm;

		// struct for storing wavefront index data
		struct obj_vertex {
			int p = 0, n = 0, t = 0;
		};

		// create reading buffers
		vector<vec3> positions;
		vector<vec3> normals;
		vector<vec2> uvs;
		vector<obj_vertex> wv_vertices;

		// open file
		ifstream obj_file(filepath);
		if (!obj_file.is_open()) {
			cerr << "Error: could not open " << filepath << endl;
			throw runtime_error("could not open " + filepath.string());
		}

		// good() means that failbit, badbit and eofbit are all not set
		while (obj_file.good()) {

			// pull out line from file
			string line;
			getline(obj_file, line);
			istringstream obj_line(line);

			// pull out mode from line
			string mode;
			obj_line >> mode;

			// reading like this means whitespace at the start of the line is fine
			// attempting to read from an empty string/line will set the failbit
			if (obj_line.good()) {
				if (mode == "v") {
					vec3 v;
					obj_line >> v.x >> v.y >> v.z;
					positions.push_back(v);
				} else if (mode == "vn") {
					vec3 vn;
					obj_line >> vn.x >> vn.y >> vn.z;
					normals.push_back(vn);
				} else if (mode == "vt") {
					vec2 vt;
					obj_line >> vt.x >> vt.y;
					uvs.push_back(vt);
				} else if (mode == "f") {
					std::vector<obj_vertex> face;
					while (obj_line.good()) {
						obj_vertex v;
						// scan in position index
						obj_line >> v.p;
						if (obj_line.fail()) break;
						// look ahead for a match
						if (obj_line.peek() == '/') {	
							// ignore the '/' character
							obj_line.ignore(1);			
							// scan in uv (texture coord) index (if it's there)
							if (obj_line.peek() != '/') {
								obj_line >> v.t;
							}
							// scan in normal index (if it's there)
							if (obj_line.peek() == '/') {
								obj_line.ignore(1);
								obj_line >> v.n;
							}
						}
						// subtract one because of wavefront indexing
						v.p -= 1;
						v.n -= 1;
						v.t -= 1;
						face.push_back(v);
					}
					// IFF we have 3 vertices, construct a triangle
					if (face.size() == 3) {
						for (int i = 0; i < 3; ++i) {
							wv_vertices.push_back(face[i]);
						}
					} else {
						cerr << "OBj file " << filepath << " contains non-triangle face" << endl;
					}
				}
			}
		}

		// if we don't have any normals, create them naively
		if (normals.empty()) {
			// create the normals as 3d vectors of 0
			normals.resize(positions.size(), vec3(0));
			// add the normal for every face to each vertex-normal
			for (size_t i = 0; i < wv_vertices.size() / 3; i++) {
				obj_vertex &a = wv_vertices[i * 3];
				obj_vertex &b = wv_vertices[i * 3 + 1];
				obj_vertex &c = wv_vertices[i * 3 + 2];
				// set the normal index to be the same as position index
				a.n = a.p;
				b.n = b.p;
				c.n = c.p;
				// calculate the face normal
				vec3 ab = positions[b.p] - positions[a.p];
				vec3 ac = positions[c.p] - positions[a.p];
				vec3 face_norm = cross(ab, ac);
				// contribute the face norm to each vertex
				float l = length(face_norm);
				if (l > 0) {
					face_norm / l;
					normals[a.n] += face_norm;
					normals[b.n] += face_norm;
					normals[c.n] += face_norm;
				}
			}
			// normalize the normals
			for (auto &n : normals) {
				n = normalize(n);
			}
		}

		// TODO create spherical UV's if they don't exist

		// create mesh data
		mesh_builder<> mb;
		for (unsigned int i = 0; i < wv_vertices.size(); ++i) {
			mb.push_index(i);
			mb.push_vertex({
				positions[wv_vertices[i].p],
				normals[wv_vertices[i].n],
				uvs[wv_vertices[i].t]
			});
		}
		return mb;

	}

}
