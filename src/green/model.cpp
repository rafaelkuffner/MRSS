
#include "model.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <glm/gtc/type_ptr.hpp>

#include <cgu/shader.hpp>

namespace green {

	Model::Model(const std::filesystem::path &fpath) {
		m_trimesh.request_face_normals();
		m_trimesh.request_vertex_normals();
		// load custom "quality" property, if it exists
		OpenMesh::IO::Options readOptions = OpenMesh::IO::Options::Custom;
		// TODO unicode filenames...
		std::cerr << "Loading model " << fpath << std::endl;
		if (OpenMesh::IO::read_mesh(m_trimesh, fpath.string(), readOptions)) {
			std::cerr << "Loaded" << std::endl;
		} else {
			std::cerr << "Failed" << std::endl;
			throw std::runtime_error("failed to load model");
		}
		m_trimesh.triangulate();
		// calculate normals if missing (note: need face normals to calc vertex normals)
		if (!readOptions.face_has_normal()) m_trimesh.update_face_normals();
		if (!readOptions.vertex_has_normal()) m_trimesh.update_vertex_normals();
		// bounding box
		for (auto vit = m_trimesh.vertices_begin(); vit != m_trimesh.vertices_end(); ++vit) {
			m_bound_min = min(m_bound_min, om2glm(m_trimesh.point(*vit)));
			m_bound_max = max(m_bound_max, om2glm(m_trimesh.point(*vit)));
		}

		// TODO loading on background thread
		// gl mesh upload has to happen on main thread

	}

	void Model::update_vao() {
		bool makevao = false;
		if (!m_vao) makevao = true, m_vao = cgu::gl_object::gen_vertex_array();
		if (!m_ibo) m_ibo = cgu::gl_object::gen_buffer();
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		const size_t ntris = m_trimesh.n_faces();
		assert(ntris <= size_t(INT_MAX));
		m_vao_ntris = ntris;
		const size_t size_ibo = ntris * 3 * sizeof(GLuint);
		glBindVertexArray(m_vao);
		// the GL_ELEMENT_ARRAY_BUFFER binding sticks to the VAO so we shouldn't unbind it
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_ibo, nullptr, GL_STATIC_DRAW);
		if (ntris) {
			auto pibo = reinterpret_cast<GLuint *>(
				glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, size_ibo, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
			);
			for (auto fit = m_trimesh.faces_begin(); fit != m_trimesh.faces_end(); ++fit) {
				auto hit = m_trimesh.cfv_iter(*fit);
				// assume triangles
				*pibo++ = GLuint(hit++->idx());
				*pibo++ = GLuint(hit++->idx());
				*pibo++ = GLuint(hit++->idx());
			}
			glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		}
		if (makevao) {
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, false, 12, 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, false, 12, 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vbos() {
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		m_trimesh.points();
		m_trimesh.vertex_normals();
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * 12, m_trimesh.points(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * 12, m_trimesh.vertex_normals(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void Model::draw() const {
		if (!m_vao) return;
		glBindVertexArray(m_vao);
		glDrawElements(GL_TRIANGLES, m_vao_ntris * 3, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	void Model::draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar) const {

		static GLuint prog = 0;
		if (!prog) {
			std::ostringstream hdr;
			hdr << "#ifdef _FRAGMENT_\n#define CGU_NEED_FRAG_DEPTH\n" << cgu::glsl_frag_depth_source << "\n#endif\n";
			prog = cgu::make_shader_program_from_files(
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{hdr.str()},
				{"./res/model.glsl"}
			).release();
		}

		glUseProgram(prog);
		glUniform1f(glGetUniformLocation(prog, cgu::glsl_uniform_zfar_name), zfar);
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(modelview));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(projection));

		draw();
	}

}
