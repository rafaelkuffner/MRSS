
#include "model.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <chrono>
#include <mutex>

#include <glm/gtc/type_ptr.hpp>

#include <cgu/shader.hpp>

#include <imgui.h>

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
		const size_t nverts = m_trimesh.n_vertices();
		assert(ntris <= size_t(INT_MAX));
		assert(nverts <= size_t(INT_MAX));
		m_vao_ntris = ntris;
		m_vao_nverts = nverts;
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

	void Model::draw(GLenum polymode) const {
		if (!m_vao) return;
		glBindVertexArray(m_vao);
		if (polymode == GL_POINT) {
			glDrawArrays(GL_POINTS, 0, m_vao_nverts);
		} else {
			glPolygonMode(GL_FRONT_AND_BACK, polymode);
			glDrawElements(GL_TRIANGLES, m_vao_ntris * 3, GL_UNSIGNED_INT, 0);
		}
		glBindVertexArray(0);
	}

	void Model::draw(const glm::mat4 &modelview, const glm::mat4 &projection, float zfar, const model_draw_params &params, GLenum polymode) const {

		static GLuint prog = 0;
		if (!prog) {
			// note: needs custom depth env
			prog = cgu::make_shader_program_from_files(
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{cgu::glsl_frag_depth_source},
				{"./res/model.glsl"}
			).release();
		}

		glUseProgram(prog);
		glUniform1f(glGetUniformLocation(prog, cgu::glsl_uniform_zfar_name), zfar);
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(modelview));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(projection));
		glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, value_ptr(params.color));
		glUniform4fv(glGetUniformLocation(prog, "u_color_hover"), 1, value_ptr(params.color_hover));
		glUniform4fv(glGetUniformLocation(prog, "u_color_select"), 1, value_ptr(params.color_select));
		glm::vec3 pos_bias = -bound_center();
		glUniform3fv(glGetUniformLocation(prog, "u_pos_bias"), 1, value_ptr(pos_bias));
		glUniform1f(glGetUniformLocation(prog, "u_shading"), params.shading);
		float bias = 0;
		// biases adjusted for log-depth
		if (polymode == GL_LINE) bias = -0.00001;
		if (polymode == GL_POINT) bias = -0.00002;
		glUniform1f(glGetUniformLocation(prog, "u_depth_bias"), bias);
		glUniform1i(glGetUniformLocation(prog, "u_entity_id"), params.entity_id);
		glUniform4iv(glGetUniformLocation(prog, "u_selection"), 1, (GLint *) &params.sel);

		draw(polymode);
	}

	ModelEntity::ModelEntity() {
		
	}

	void ModelEntity::load(const std::filesystem::path &fpath) {
		static std::mutex load_mtx;
		m_model.reset();
		m_fpath = fpath;
		m_pending = std::async([=]() {
			// apparently openmesh load is not threadsafe?
			// TODO check assimp too
			std::lock_guard lock(load_mtx);
			return std::make_unique<Model>(fpath);
		});
	}

	glm::mat4 ModelEntity::transform() const {
		if (!m_model) return glm::mat4(1);
		glm::mat4 transform(1);
		transform = glm::translate(transform, m_translation);
		transform = glm::scale(transform, glm::vec3(m_scale));
		// TODO rotation

		glm::mat3 basis(m_basis_vectors[basis_right].v, m_basis_vectors[basis_up].v, m_basis_vectors[basis_back].v);
		return transform * glm::mat4(transpose(basis));
	}

	void ModelEntity::draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, const selection &sel) {
		if (m_pending.valid() && m_pending.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
			try {
				// is this the best place to do this?
				m_model = m_pending.get();
				// gl stuff has to happen on main thread
				m_model->update_vao();
				m_model->update_vbos();
				m_scale = m_model->unit_bound_scale() * 4;
				m_translation = glm::vec3(0);
			} catch (...) {
				// load failed
			}
		}
		if (ImGui::Begin("Models")) {
			// TODO unicode...
			if (ImGui::CollapsingHeader(m_fpath.filename().string().c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
				if (ImGui::IsItemHovered()) ImGui::SetTooltip(m_fpath.string().c_str());
				ImGui::PushID(m_model.get());
				if (m_pending.valid()) {
					ImGui::Text("Loading...");
				} else if (m_model) {
					ImGui::Text("%zd vertices, %zd triangles", m_model->trimesh().n_vertices(), m_model->trimesh().n_faces());
				} else {
					ImGui::Text("Failed to load model");
				}
				ImGui::Checkbox("Faces", &m_show_faces);
				ImGui::SameLine();
				ImGui::Checkbox("Edges", &m_show_edges);
				ImGui::SameLine();
				ImGui::Checkbox("Verts", &m_show_verts);
				auto pick_basis = [&](const char *label, int *basis) {
					ImGui::SetNextItemWidth(ImGui::GetTextLineHeight() * 3.5f);
					ImGui::Combo(
						label, basis,
						[](void *data, int item, const char **out_text) {
							auto vs = reinterpret_cast<basis_vector *>(data);
							*out_text = vs[item].name;
							return true;
						}, m_basis_vectors, 6
					);
				};
				pick_basis("Right", &basis_right);
				ImGui::SameLine();
				pick_basis("Up", &basis_up);
				ImGui::SameLine();
				pick_basis("Back", &basis_back);
				ImGui::SliderFloat("Scale", &m_scale, 0, 1000, "%.4f", 8);
				ImGui::SliderFloat3("Translation", value_ptr(m_translation), -10, 10);
				if (ImGui::Button("Reset Scale")) m_scale = m_model->unit_bound_scale() * 4;
				ImGui::SameLine();
				if (ImGui::Button("Remove")) ImGui::OpenPopup("##remove");
				if (ImGui::BeginPopup("##remove", ImGuiWindowFlags_Modal)) {
					ImGui::Text("Remove model \"%s\" ?", m_fpath.filename().string().c_str());
					if (ImGui::IsItemHovered()) ImGui::SetTooltip(m_fpath.string().c_str());
					if (ImGui::Button("Remove")) {
						m_dead = true;
						ImGui::CloseCurrentPopup();
					}
					ImGui::SameLine();
					if (ImGui::Button("Cancel")) ImGui::CloseCurrentPopup();
					ImGui::EndPopup();
				}
				ImGui::PopID();
			}
		}
		ImGui::End();
		if (m_model) {
			model_draw_params params;
			params.sel = sel;
			params.entity_id = id();
			params.color = {0.6f, 0.6f, 0.5f, 1};
			params.color_hover = {0.7f, 0.7f, 0.3f, 1};
			params.color_select = {0.7f, 0.4f, 0.1f, 1};
			glEnable(GL_CULL_FACE);
			glColorMaski(1, GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
			if (m_show_faces) m_model->draw(view * transform(), proj, zfar, params, GL_FILL);
			params.shading = 0;
			params.color = {0.03f, 0.03f, 0.03f, 1};
			params.color_hover = params.color;
			params.color_select = params.color;
			glDisable(GL_CULL_FACE);
			glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			if (m_show_edges) m_model->draw(view * transform(), proj, zfar, params, GL_LINE);
			params.color = {0.5f, 0, 0, 1};
			params.color_hover = {1, 1, 0, 1};
			params.color_select = {0.5f, 1, 0, 1};
			params.sel.hover_entity = -1;
			params.sel.select_entity = -1;
			glColorMaski(1, GL_TRUE, GL_TRUE, GL_FALSE, GL_FALSE);
			if (m_show_verts) m_model->draw(view * transform(), proj, zfar, params, GL_POINT);
			glEnable(GL_CULL_FACE);
		}
	}

	ModelEntity::~ModelEntity() {

	}

}
