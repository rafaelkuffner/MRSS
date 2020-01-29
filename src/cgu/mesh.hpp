
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

#pragma once

#include <cassert>
#include <cstdint>
#include <vector>
#include <utility>
#include <initializer_list>
#include <filesystem>

#include <glm/glm.hpp>

#include "opengl.hpp"

namespace cgu {

	// non-owning slice of a gl_mesh.
	// draws some subrange of the indices in the mesh.
	struct gl_mesh_slice {
		GLuint vao = 0;
		int indexstart = 0;
		int indexcount = 0;
		GLenum mode = 0;

		// use glPolygonMode if wireframe is desired
		void draw() const;
	};

	struct [[nodiscard]] gl_mesh {
		gl_object vao;
		gl_object vbo;
		gl_object ibo;
		int indexcount = 0;
		GLenum mode = 0;

		// use glPolygonMode if wireframe is desired
		void draw() const {
			slice(0, indexcount).draw();
		}

		explicit operator bool() const {
			return bool(vao);
		}

		bool operator!() const {
			return !bool(*this);
		}

		gl_mesh_slice slice(int start, int count) const {
			assert(0 <= start && start + count <= indexcount);
			gl_mesh_slice ms;
			ms.mode = mode;
			ms.indexstart = start;
			ms.indexcount = count;
			ms.vao = vao;
			return ms;
		}
	};

	enum class vertex_attrib_type : uint8_t {
		float_unnormalized, float_normalized, integer
	};

	struct vertex_attrib {
		vertex_attrib_type atttype = vertex_attrib_type::float_unnormalized;
		int8_t components;
		GLenum srctype = GL_FLOAT;
		size_t offset = 0;
	};

	struct mesh_vertex {
		glm::vec3 pos{0};
		glm::vec3 norm{0};
		glm::vec2 uv{0};
	};

	template <typename T>
	struct mesh_vertex_traits {};

	template <>
	struct mesh_vertex_traits<mesh_vertex> {
		static std::vector<vertex_attrib> attribs() {
			return {
				{vertex_attrib_type::float_unnormalized, 3, GL_FLOAT, offsetof(mesh_vertex, pos)},
				{vertex_attrib_type::float_unnormalized, 3, GL_FLOAT, offsetof(mesh_vertex, norm)},
				{vertex_attrib_type::float_unnormalized, 2, GL_FLOAT, offsetof(mesh_vertex, uv)}
			};
		}
	};

	struct mesh_builder_base {
		GLenum mode = GL_TRIANGLES;
		GLenum bufhint = GL_STATIC_DRAW;
		std::vector<vertex_attrib> attribs;
		std::vector<GLuint> indices;

		void push_index(GLuint i) {
			indices.push_back(i);
		}

		void push_indices(std::initializer_list<GLuint> inds) {
			indices.insert(indices.end(), inds);
		}

		gl_mesh build_impl(const void *data, size_t size, size_t stride, gl_mesh m) const;
	};

	template <typename T = mesh_vertex>
	struct mesh_builder : mesh_builder_base {
		std::vector<T> vertices;
		
		mesh_builder() {
			attribs = mesh_vertex_traits<T>::attribs();
		}

		explicit mesh_builder(GLenum mode_) {
			mode = mode_;
			attribs = mesh_vertex_traits<T>::attribs();
		}

		void clear() {
			vertices.clear();
			indices.clear();
		}

		GLuint push_vertex(const T &v) {
			auto size = vertices.size();
			assert(size == decltype(size)(GLuint(size)));
			vertices.push_back(v);
			return GLuint(size);
		}

		gl_mesh build(gl_mesh m = {}) {
			return build_impl(vertices.data(), vertices.size() * sizeof(T), sizeof(T), std::move(m));
		}
	};

	mesh_builder<> mesh_load_obj(const std::filesystem::path &);

}
