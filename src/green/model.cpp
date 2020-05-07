
#include "model.hpp"

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <chrono>
#include <mutex>
#include <algorithm>

#include <glm/gtc/type_ptr.hpp>

#include <cgu/shader.hpp>

#include <imgui.h>

#include "main.hpp"
#include "model.glsl.hpp"

namespace {

	auto & mesh_io_mutex() {
		static std::mutex m;
		return m;
	}

	struct saliency_clipboard_data {
		green::model_saliency_data sd;
		std::vector<float> data;
	};

	auto & saliency_clipboard() {
		static saliency_clipboard_data d;
		return d;
	}
}

namespace green {

	Model::Model(const std::filesystem::path &fpath) {
		m_trimesh.request_face_normals();
		m_trimesh.request_vertex_normals();
		m_trimesh.request_vertex_colors();
		// load custom "quality" property and vertex colors if they exist
		// FIXME what else do we need to ask for and preserve on export? texcoords?
		OpenMesh::IO::Options readOptions = OpenMesh::IO::Options::Custom | OpenMesh::IO::Options::VertexColor;

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

		std::cout << "Computing vertex areas" << std::endl;
		m_prop_vertex_area = computeVertexAreas(m_trimesh);

		std::cout << "Computing edge lengths" << std::endl;
		m_prop_edge_length = computeEdgeLengths(m_trimesh);

		// copy original vertex colors
		// (because we need to be able to overwrite the actual vertex colors during export)
		if (readOptions.check(OpenMesh::IO::Options::VertexColor)) {
			std::cout << "Found vertex colors" << std::endl;
			m_trimesh.add_property(m_prop_vcolor_original);
			for (auto vIt = m_trimesh.vertices_begin(); vIt != m_trimesh.vertices_end(); ++vIt) {
				m_trimesh.property(m_prop_vcolor_original, *vIt) = m_trimesh.color(*vIt);
			}
		}

		using sprop_t = decltype(model_saliency_data::prop_saliency);
		std::vector<std::pair<sprop_t, model_saliency_data>> sprops;

		// check for 'raw' saliency properties
		for (auto it = m_trimesh.vprops_begin(); it != m_trimesh.vprops_end(); ++it) {
			// TODO any way to do this other than dynamic cast?
			auto prop = dynamic_cast<OpenMesh::PropertyT<float> *>(*it);
			if (!prop) continue;
			if (prop->name().substr(0, 7) == "quality") {
				std::cout << "Found quality property " << prop->name() << std::endl;
				auto ph0 = sprop_t{int(it - m_trimesh.vprops_begin())};
				if (&m_trimesh.property(ph0) != prop) abort();
				model_saliency_data sd;
				sd.filename = fpath.filename().u8string();
				sd.propname = prop->name();
				sprops.emplace_back(ph0, std::move(sd));
			}
		}

		// move 'raw' saliency to unnamed properties
		// need 2 passes because we cant modify the properties while iterating them
		for (auto &p : sprops) {
			sprop_t ph;
			m_trimesh.add_property(ph);
			m_trimesh.property(ph).data_vector() = std::move(m_trimesh.property(p.first).data_vector());
			m_trimesh.remove_property(p.first);
			auto sd = std::move(p.second);
			sd.prop_saliency = ph;
			m_original_saliency.push_back(std::move(sd));
		}

	}

	void Model::save(const std::filesystem::path &fpath, OpenMesh::VPropHandleT<float> prop_saliency) {

		if (m_trimesh.has_vertex_colors() && prop_saliency.is_valid()) {
			std::cout << "Colorizing saliency" << std::endl;
			for (auto vIt = m_trimesh.vertices_begin(); vIt != m_trimesh.vertices_end(); ++vIt) {
				TriMesh::Color col;
				const float s = m_trimesh.property(prop_saliency, *vIt);
				OpenMesh::Vec3f v;
				mapScalarToColor(v, s, TransferFunction::ZBRUSH);
				v = v * 255;
				col[0] = v[0];
				col[1] = v[1];
				col[2] = v[2];
				m_trimesh.set_color(*vIt, col);
			}
		}

		// TODO save arbitrary saliency properties with arbitrary names
		auto prop_quality = decltype(prop_saliency){};
		if (prop_saliency.is_valid()) {
			m_trimesh.add_property(prop_quality, "quality");
			m_trimesh.property(prop_quality).set_persistent(true);
			m_trimesh.property(prop_quality).data_vector() = m_trimesh.property(prop_saliency).data_vector();
		}

		OpenMesh::IO::Options opts = OpenMesh::IO::Options::Custom;
		if (m_trimesh.has_vertex_colors()) opts = opts | OpenMesh::IO::Options::VertexColor;

		// TODO unicode filenames...
		auto res = OpenMesh::IO::write_mesh(m_trimesh, fpath.string(), opts);

		m_trimesh.remove_property(prop_quality);

		if (!res) {
			std::cerr << "Could not write mesh file " << fpath << std::endl;
			throw std::runtime_error("failed to save model");
		}
	}

	void Model::update_vao() {
		bool makevao = false;
		if (!m_vao) makevao = true, m_vao = cgu::gl_object::gen_vertex_array();
		if (!m_ibo) m_ibo = cgu::gl_object::gen_buffer();
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
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
			glVertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, false, sizeof(glm::vec3), 0);
			glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 4, GL_FLOAT, false, sizeof(glm::vec4), 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	void Model::update_vbos() {
		if (!m_vbo_pos) m_vbo_pos = cgu::gl_object::gen_buffer();
		if (!m_vbo_norm) m_vbo_norm = cgu::gl_object::gen_buffer();
		if (!m_vbo_col) m_vbo_col = cgu::gl_object::gen_buffer();
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec3), m_trimesh.points(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_norm);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec3), m_trimesh.vertex_normals(), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_col);
		glBufferData(GL_ARRAY_BUFFER, m_trimesh.n_vertices() * sizeof(glm::vec4), nullptr, GL_DYNAMIC_DRAW);
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
			prog = cgu::make_shader_program(
				"green::model",
				"330 core",
				{GL_VERTEX_SHADER, GL_FRAGMENT_SHADER},
				{cgu::glsl_frag_depth_source, cgu::strings::glsl_green_model}
			).release();
		}

		glUseProgram(prog);
		glUniform1f(glGetUniformLocation(prog, cgu::glsl_uniform_zfar_name), zfar);
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_modelview"), 1, false, value_ptr(modelview));
		glUniformMatrix4fv(glGetUniformLocation(prog, "u_projection"), 1, false, value_ptr(projection));
		glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, value_ptr(params.color));
		glm::vec3 pos_bias = -bound_center();
		glUniform3fv(glGetUniformLocation(prog, "u_pos_bias"), 1, value_ptr(pos_bias));
		glUniform1f(glGetUniformLocation(prog, "u_shading"), params.shading);
		float bias = 0;
		// biases adjusted for log-depth
		if (polymode == GL_LINE) bias = -0.00001;
		if (polymode == GL_POINT) bias = -0.00002;
		glUniform1f(glGetUniformLocation(prog, "u_depth_bias"), bias);
		glUniform1i(glGetUniformLocation(prog, "u_entity_id"), params.entity_id);
		glUniform1i(glGetUniformLocation(prog, "u_vert_color_map"), params.vert_color_map);
		glUniform1i(glGetUniformLocation(prog, "u_show_samples"), params.show_samples);

		draw(polymode);
	}

	ModelEntity::ModelEntity() {
		
	}

	void ModelEntity::load(const std::filesystem::path &fpath) {
		m_model.reset();
		m_fpath_load = fpath;
		m_pending_load = std::async([=]() {
			// apparently openmesh load is not threadsafe?
			// TODO check assimp too
			std::lock_guard lock(mesh_io_mutex());
			return std::make_unique<Model>(fpath);
		});
	}

	void ModelEntity::save(const std::filesystem::path &fpath) {
		const auto salout = m_saliency_index < m_saliency_outputs.size() ? m_saliency_outputs[m_saliency_index] : model_saliency_data{};
		if (salout.prop_saliency.is_valid()) m_model->trimesh().request_vertex_colors();
		m_fpath_save = fpath;
		m_pending_save = std::async([=]() {
			// TODO is this also not threadsafe? idk
			// ... can it run parallel with load?
			std::lock_guard lock(mesh_io_mutex());
			m_model->save(fpath, salout.prop_saliency);
			return true;
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

	void ModelEntity::update_vbo() {
		if (!m_model) return;
		if (!m_saliency_vbo_dirty) return;
		if (m_color_mode == color_mode::saliency && m_saliency_index >= m_saliency_outputs.size()) return;
		if (m_color_mode == color_mode::vcolor && !m_model->prop_vcolor_original().is_valid()) return;
		if (m_color_mode == color_mode::saliency_comparison && m_saliency_index >= m_saliency_outputs.size()) return;
		if (m_color_mode == color_mode::saliency_comparison && m_saliency_baseline_index >= m_saliency_outputs.size()) return;
		GLuint vbo_col = m_model->vbo_color();
		if (!vbo_col) return;
		const auto &mesh = m_model->trimesh();
		const auto nverts = m_model->vao_nverts();
		glBindBuffer(GL_ARRAY_BUFFER, vbo_col);
		auto *data = reinterpret_cast<glm::vec4 *>(
			glMapBufferRange(GL_ARRAY_BUFFER, 0, nverts * sizeof(glm::vec4), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
		);
		if (m_color_mode == color_mode::saliency) {
			auto salprop = m_saliency_outputs[m_saliency_index].prop_saliency;
			auto sampledprop = m_saliency_outputs[m_saliency_index].prop_sampled;
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float s = mesh.property(salprop, v);
				const bool sampled = sampledprop.is_valid() && mesh.property(sampledprop, v);
				data[i] = glm::vec4(s, 0, 0, sampled);
			}
		} else if (m_color_mode == color_mode::saliency_comparison) {
			auto baseprop = m_saliency_outputs[m_saliency_baseline_index].prop_saliency;
			auto salprop = m_saliency_outputs[m_saliency_index].prop_saliency;
			auto sampledprop = m_saliency_outputs[m_saliency_index].prop_sampled;
			m_saliency_errors.min = 9001;
			m_saliency_errors.max = -9001;
			float sse = 0;
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const float b = mesh.property(baseprop, v);
				const float s = mesh.property(salprop, v);
				const bool sampled = mesh.property(sampledprop, v);
				const float e = s - b;
				m_saliency_errors.min = std::min(m_saliency_errors.min, e);
				m_saliency_errors.max = std::max(m_saliency_errors.max, e);
				sse += e * e;
				// TODO do scale on gpu
				data[i] = glm::vec4(e * m_saliency_error_scale, 0, 0, sampled);
			}
			m_saliency_errors.rms = sqrt(sse / nverts);
		} else if (m_color_mode == color_mode::vcolor) {
			for (size_t i = 0; i < nverts; i++) {
				const auto v = OpenMesh::VertexHandle(i);
				const auto col = mesh.property(m_model->prop_vcolor_original(), v);
				// need to gamma decode the color because the shader expects linear
				data[i] = glm::vec4(pow(glm::vec3(col[0], col[1], col[2]) / 255.f, glm::vec3(2.2f)), 1);
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		m_saliency_vbo_dirty = false;
	}

	void ModelEntity::draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar) {
		if (m_pending_load.valid() && m_pending_load.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
			try {
				// is this the best place to do this?
				m_model = m_pending_load.get();
				// gl stuff has to happen on main thread
				m_model->update_vao();
				m_model->update_vbos();
				m_scale = m_model->unit_bound_scale() * 4;
				m_translation = glm::vec3(0);
				// add saliency results for loaded saliency
				for (auto &sd : m_model->original_saliency()) {
					m_saliency_outputs.push_back(sd);
					m_saliency_vbo_dirty = true;
				}
			} catch (...) {
				// load failed
			}
		}
		auto &sel = ui_selection();
		const bool selected = sel.select_entity == id();
		if (ImGui::Begin("Models")) {
			// TODO unicode...
			if (selected) ImGui::PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
			if (ImGui::CollapsingHeader(m_fpath_load.filename().string().c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
				if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", m_fpath_load.string().c_str());
				ImGui::PushID(m_model.get());
				if (m_pending_load.valid()) {
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
				ImGui::SameLine();
				if (ImGui::Button("Remove")) ImGui::OpenPopup("##remove");
				ImGui::SameLine();
				if (ImGui::Button("Select")) {
					sel.select_entity = id();
					sel.select_vertex = -1;
				}
				if (ImGui::BeginPopup("##remove", ImGuiWindowFlags_Modal)) {
					ImGui::Text("Remove model \"%s\" ?", m_fpath_load.filename().string().c_str());
					if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", m_fpath_load.string().c_str());
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
			if (selected) ImGui::PopStyleColor();
		}
		ImGui::End();
		if (ImGui::Begin("Selection")) {
			if (selected) {
				
				ImGui::PushStyleColor(ImGuiCol_Header, {0.7f, 0.4f, 0.1f, 1});
				// TODO unicode?
				ImGui::Selectable(m_fpath_load.filename().string().c_str(), true);
				ImGui::PopStyleColor();
				if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", m_fpath_load.string().c_str());
				ImGui::PushID(m_model.get());

				if (ImGui::Button("Export...")) ImGui::OpenPopup("Export##export");
				ImGui::SameLine();
				if (m_pending_save.valid()) {
					if (m_pending_save.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
						try {
							m_save_ok = m_pending_save.get();
						} catch (...) {
							// save failed
							m_save_ok = false;
						}
					} else {
						// TODO unicode?
						ImGui::Text("Exporting %s...", m_fpath_save.filename().u8string().c_str());
					}
				} else if (!m_fpath_save.empty()) {
					ImGui::Text("Export %s %s", m_fpath_save.filename().u8string().c_str(), m_save_ok ? "succeeded" : "failed");
				} else {
					ImGui::Text("");
				}
				// TODO unicode?
				if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", m_fpath_save.u8string().c_str());

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
				if (ImGui::Button("Reset##scale")) m_scale = m_model->unit_bound_scale() * 4;
				ImGui::SameLine();
				ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x * 0.6f);
				ImGui::SliderFloat("Scale", &m_scale, 0, 1000, "%.4f", 8);
				if (ImGui::Button("Reset##translation")) m_translation = glm::vec3{0};
				ImGui::SameLine();
				ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x * 0.6f);
				ImGui::InputFloat3("Translation", value_ptr(m_translation));

				ImGui::SetNextWindowSize({400, 150}, ImGuiCond_Appearing);
				if (ImGui::BeginPopupModal("Export##export")) {
					ImGui::Selectable(m_fpath_load.filename().string().c_str());
					if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", m_fpath_load.string().c_str());
					if (m_color_mode == color_mode::saliency && m_saliency_index < m_saliency_outputs.size()) {
						auto &salout = m_saliency_outputs[m_saliency_index];
						ImGui::Text("Saliency: %s", std::string(salout).c_str());
					} else {
						ImGui::Text("No saliency will be exported");
					}
					auto cwd = std::filesystem::current_path();
					ImGui::Text("CWD: %s", cwd.u8string().c_str());
					// TODO better?
					static char buf[1024];
					if (ImGui::IsWindowAppearing()) {
						snprintf(buf, sizeof(buf), "%s", m_fpath_save.u8string().c_str());
					}
					ImGui::InputText("Path", buf, sizeof(buf));
					// TODO unicode?
					auto fpath = std::filesystem::u8path(buf);
					auto stat = std::filesystem::status(fpath);
					bool cansave = !fpath.empty();
					if (std::filesystem::is_directory(stat)) {
						ImGui::Text("Path is a directory");
						cansave = false;
					} else if (std::filesystem::exists(stat)) {
						ImGui::Text("Path exists! Save will overwrite.");
					}
					if (ImGui::Button("Save") && cansave && !m_pending_save.valid()) {
						save(std::filesystem::absolute(fpath));
						ImGui::CloseCurrentPopup();
					}
					if (ImGui::IsItemHovered()) {
						if (m_pending_save.valid()) {
							ImGui::SetTooltip("Export already in progress");
						} else if (!cansave) {
							ImGui::SetTooltip("Can't save to this path");
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Cancel")) ImGui::CloseCurrentPopup();
					ImGui::EndPopup();
				}
				
				ImGui::Separator();

				//ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
				int cur_color_mode = int(m_color_mode);
				if (ImGui::Combo("Color Mode", &cur_color_mode, "None\0Vertex Color\0Saliency\0Saliency Comparison\0")) {
					m_color_mode = color_mode(cur_color_mode);
					if (m_color_mode != color_mode::none) m_saliency_vbo_dirty = true;
				}
				
				if (ImGui::Combo(
					"Saliency", &m_saliency_index,
					[](void *data, int item, const char **out_text) {
						auto ds = reinterpret_cast<model_saliency_data *>(data);
						const auto &so = ds[item];
						// TODO better?
						static std::string str;
						str = std::string(so);
						*out_text = str.c_str();
						return true;
					}, m_saliency_outputs.data(), m_saliency_outputs.size()
				)) {
					m_saliency_vbo_dirty = true;
				}

				if (m_color_mode == color_mode::saliency_comparison && m_saliency_baseline_index < m_saliency_outputs.size()) {
					ImGui::Text("Baseline: %s", m_saliency_outputs[m_saliency_baseline_index].str().c_str());
					if (m_saliency_index < m_saliency_outputs.size()) {
						auto &err = m_saliency_errors;
						ImGui::Text("Errors: min=%.3f, max=%.3f, rmse=%.3f", err.min, err.max, err.rms);
						if (ImGui::SliderFloat("Error Scale", &m_saliency_error_scale, 1, 100, "%.3f", 2)) {
							m_saliency_vbo_dirty = true;
						}
					}
				}

				if (ImGui::Button("Paste")) {
					auto &clip = saliency_clipboard();
					if (clip.data.size() == m_model->trimesh().n_vertices()) {
						model_saliency_data sd = std::move(clip.sd);
						m_model->trimesh().add_property(sd.prop_saliency);
						m_model->trimesh().property(sd.prop_saliency).data_vector() = std::move(clip.data);
						clip.sd = {};
						clip.data.clear();
						m_saliency_outputs.push_back(std::move(sd));
						m_saliency_vbo_dirty = true;
						// references/iterators into saliency outputs are invalidated
					} else {
						ImGui::OpenPopup("Paste Error##pasteerror");
					}
				}

				if (m_saliency_index < m_saliency_outputs.size()) {
					auto &salout = m_saliency_outputs[m_saliency_index];
					ImGui::SameLine();
					if (ImGui::Button("Copy")) {
						auto &clip = saliency_clipboard();
						clip.sd = salout;
						clip.sd.prop_saliency.reset();
						clip.data = m_model->trimesh().property(salout.prop_saliency).data_vector();
					}
					ImGui::SameLine();
					if (ImGui::Button("Remove")) ImGui::OpenPopup("##remove");
					ImGui::SameLine();
					if (ImGui::Button("Baseline")) m_saliency_baseline_index = m_saliency_index;
					if (salout.filename.empty()) {
						ImGui::draw_saliency_params(salout.uparams);
						ImGui::draw_saliency_progress(salout.progress);
						if (ImGui::Button("Reload Parameters")) ui_saliency_user_params() = salout.uparams;
					}
					if (ImGui::BeginPopupModal("Paste Error##pasteerror")) {
						auto &clip = saliency_clipboard();
						if (clip.data.size() == m_model->trimesh().n_vertices()) {
							ImGui::CloseCurrentPopup();
						} else {
							ImGui::Text("Can't paste saliency for %d vertices into model with %d vertices", clip.data.size(), m_model->trimesh().n_vertices());
							if (ImGui::Button("OK")) ImGui::CloseCurrentPopup();
						}
					}
					if (ImGui::BeginPopup("##remove", ImGuiWindowFlags_Modal)) {
						ImGui::Text("Remove saliency result?");
						if (ImGui::Button("Remove")) {
							auto it = m_saliency_outputs.begin() + m_saliency_index;
							m_model->trimesh().remove_property(it->prop_saliency);
							m_saliency_outputs.erase(it);
							// references/iterators into saliency outputs are invalidated
							if (m_saliency_index >= m_saliency_outputs.size()) {
								m_saliency_index = std::max(0, m_saliency_index - 1);
								m_saliency_vbo_dirty = true;
							}
							ImGui::CloseCurrentPopup();
						}
						ImGui::SameLine();
						if (ImGui::Button("Cancel")) ImGui::CloseCurrentPopup();
						ImGui::EndPopup();
					}
				}
				ImGui::PopID();
			}
		}
		ImGui::End();
		if (m_model) {
			update_vbo();
			model_draw_params params;
			params.sel = sel;
			params.entity_id = id();
			// faces
			params.color = {0.6f, 0.6f, 0.5f, 1};
			const bool sal_valid = m_saliency_index < m_saliency_outputs.size();
			if (m_color_mode == color_mode::saliency && sal_valid) params.vert_color_map = 3;
			if (m_color_mode == color_mode::saliency_comparison && sal_valid && m_saliency_baseline_index < m_saliency_outputs.size()) params.vert_color_map = 4;
			if (m_color_mode == color_mode::vcolor && m_model->prop_vcolor_original().is_valid()) params.vert_color_map = 1;
			glEnable(GL_CULL_FACE);
			glColorMaski(1, GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
			if (m_show_faces) m_model->draw(view * transform(), proj, zfar, params, GL_FILL);
			// edges
			params.vert_color_map = 0;
			params.shading = 0;
			params.color = {0.03f, 0.03f, 0.03f, 1};
			glDisable(GL_CULL_FACE);
			glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			if (m_show_edges) m_model->draw(view * transform(), proj, zfar, params, GL_LINE);
			// verts
			params.color = {0.5f, 0, 0, 1};
			params.sel.hover_entity = -1;
			params.sel.select_entity = -1;
			params.show_samples = sal_valid && (m_color_mode == color_mode::saliency || m_color_mode == color_mode::saliency_comparison);
			glColorMaski(1, GL_TRUE, GL_TRUE, GL_FALSE, GL_FALSE);
			if (m_show_verts) m_model->draw(view * transform(), proj, zfar, params, GL_POINT);
			glEnable(GL_CULL_FACE);
		}
	}

	std::future<saliency_result> ModelEntity::compute_saliency_async(const saliency_user_params &uparams, saliency_progress &progress) {
		if (!m_model) return {};
		saliency_mesh_params mparams;
		mparams.mesh = &m_model->trimesh();
		mparams.prop_vertex_area = m_model->prop_vertex_area();
		mparams.prop_edge_length = m_model->prop_edge_length();
		// create properties
		mparams.prop_saliency_levels.resize(uparams.levels);
		mparams.mesh->add_property(mparams.prop_curvature);
		mparams.mesh->add_property(mparams.prop_saliency);
		mparams.mesh->add_property(mparams.prop_sampled);
		for (int i = 0; i < uparams.levels; i++) {
			mparams.mesh->add_property(mparams.prop_saliency_levels[i]);
		}
		// cleanup will be run when the result is received from the future
		mparams.cleanup = [=, pprogress=&progress](bool r) mutable noexcept {
			// destroy temp properties
			mparams.mesh->remove_property(mparams.prop_curvature);
			for (int i = 0; i < uparams.levels; i++) {
				mparams.mesh->remove_property(mparams.prop_saliency_levels[i]);
			}
			if (r) {
				// if preview, remove any previous preview results
				m_saliency_outputs.erase(std::remove_if(m_saliency_outputs.begin(), m_saliency_outputs.end(), [](const auto &sd) { return sd.uparams.preview; }), m_saliency_outputs.end());
				// save user params, progress output and actual saliency mesh property
				m_saliency_outputs.push_back({"", "", uparams, *pprogress, mparams.prop_saliency, mparams.prop_sampled});
				// give focus to this result
				m_saliency_index = m_saliency_outputs.size() - 1;
				m_saliency_vbo_dirty = true;
			} else {
				// cancelled, destroy saliency property too
				mparams.mesh->remove_property(mparams.prop_saliency);
				mparams.mesh->remove_property(mparams.prop_sampled);
			}
		};
		return green::compute_saliency_async(mparams, uparams, progress);
	}

	ModelEntity::~ModelEntity() {

	}

}
