
#include "model_entity.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/euler_angles.hpp>

#include <imgui.h>

#include "imguiex.hpp"
#include "main.hpp"

namespace {
	using namespace green;

	auto & mesh_io_mutex() {
		static std::mutex m;
		return m;
	}

	struct saliency_clipboard_data {
		model_saliency_data sd;
		std::vector<float> data;
	};

	auto & saliency_clipboard() {
		static saliency_clipboard_data d;
		return d;
	}

	bool combo_color_mode(model_color_mode &mode) {
		using namespace ImGui;
		const char *optstr =
			"None\0"
			"Vertex Color\0"
			"Saliency\0"
			"Saliency Comparison\0"
			"Curvature (DoN)\0"
			"Surface Normal\0"
			"UVs\0"
			"Checkerboard\0";
		int x = int(mode);
		if (Combo("Color Mode", &x, optstr)) {
			mode = model_color_mode(x);
			return true;
		}
		return false;
	}

	GLuint checkerboard_texture() {
		static GLuint tex = 0;
		if (!tex) {
			glGenTextures(1, &tex);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, tex);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 1);
			constexpr size_t shift = 6;
			constexpr size_t sx = 1 << shift;
			std::vector<GLubyte> data(sx * sx * 4);
			for (size_t i = 0; i < sx * sx; i++) {
				const auto a = GLubyte((i & 1u) - 1u);
				const auto b = GLubyte(((i >> shift) & 1u) - 1u);
				const auto k = (a ^ b) | GLubyte(15);
				data[i * 4 + 0] = k;
				data[i * 4 + 1] = k;
				data[i * 4 + 2] = k;
				data[i * 4 + 3] = 0xFF;
			}
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, sx, sx, 0, GL_RGBA, GL_UNSIGNED_BYTE, data.data());
			glBindTexture(GL_TEXTURE_2D, 0);
		}
		return tex;
	}
}

namespace green {

	ModelEntity::ModelEntity() {

	}

	void ModelEntity::load(const std::filesystem::path &fpath) {
		// lock in case of premature closure
		// mainly because the pending_load future's dtor will block otherwise
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async([=, lock=std::move(lock)]() {
			// apparently openmesh load is not threadsafe?
			// TODO check assimp too
			std::lock_guard iolock(mesh_io_mutex());
			return std::make_unique<Model>(fpath);
		});
	}

	void ModelEntity::save(const std::filesystem::path &fpath) {
		if (!m_model) return;
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		const auto &saliency = m_model->saliency();
		const auto sdcolor = m_saliency_export_index < saliency.size() ? saliency[m_saliency_export_index] : model_saliency_data{};
		const auto sdbase = m_saliency_baseline_index < saliency.size() ? saliency[m_saliency_baseline_index] : model_saliency_data{};
		m_fpath_save = fpath;
		m_pending_save = std::async([=, lock=std::move(lock)]() {
			// TODO is this also not threadsafe? idk
			// ... can it run parallel with load?
			std::lock_guard iolock(mesh_io_mutex());
			model_save_params sparams;
			sparams.prop_saliency = sdcolor.prop_saliency;
			sparams.prop_saliency_baseline = sdbase.prop_saliency;
			sparams.color_mode = m_exp_color_mode;
			sparams.original_vids = m_save_original_vids;
			sparams.binary = m_save_binary;
			m_model->save(fpath, sparams);
			return true;
		});
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m) {
		m_model.reset();
		m_translation = glm::vec3(0);
		m_fpath_load = fpath;
		m_pending_load = std::async(std::launch::deferred, [m=std::move(m)]() mutable {
			return std::make_unique<Model>(std::move(m));
		});
		// non-blocking wait will never succeed otherwise
		m_pending_load.wait();
	}

	void ModelEntity::load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress) {
		load(fpath, std::move(m));
		m_decimated = true;
		m_src_dec_uparams = dec_uparams;
		m_src_dec_progress = dec_progress;
	}

	void ModelEntity::move_by(const glm::vec3 &d) {
		m_translation += d;
		if (d != glm::vec3(0)) invalidate_scene();
	}

	glm::mat4 ModelEntity::transform() const {
		if (!m_model) return glm::mat4(1);
		glm::mat4 transform(1);
		transform = glm::translate(transform, m_translation);
		transform *= glm::eulerAngleYXZ(m_rotation_euler_yxz.y, m_rotation_euler_yxz.x, m_rotation_euler_yxz.z);
		transform = glm::scale(transform, glm::vec3(m_scale));
		glm::mat3 basis(m_basis_vectors[m_basis_right].v, m_basis_vectors[m_basis_up].v, m_basis_vectors[m_basis_back].v);
		return transform * glm::mat4(transpose(basis));
	}

	void ModelEntity::update_vbo() {
		if (!m_model) return;
		if (!m_saliency_vbo_dirty) return;
		std::shared_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) return;
		auto &saliency_outputs = m_model->saliency();
		model_color_params cparams;
		if (m_saliency_index < saliency_outputs.size()) {
			auto &sd = saliency_outputs[m_saliency_index];
			cparams.prop_saliency = sd.prop_saliency;
			cparams.sal_gamma = sd.gamma;
		}
		if (m_saliency_baseline_index < saliency_outputs.size()) {
			auto &sd = saliency_outputs[m_saliency_baseline_index];
			cparams.prop_saliency_baseline = sd.prop_saliency;
		}
		cparams.color_mode = m_disp_color_mode;
		cparams.error_scale = m_saliency_error_scale;
		m_model->update_color(cparams, &m_saliency_errors);
		m_saliency_vbo_dirty = false;
		invalidate_scene();
	}

	void ModelEntity::draw_window_models(bool selected) {
		using namespace ImGui;
		if (Begin("Models")) {
			PushID(this);
			bool dirty = false;
			draw_select_header(selected);
			if (m_pending_load.valid()) {
				Text("Loading...");
			} else if (m_model) {
				// ok
			} else {
				Text("Failed to load model");
			}
			dirty |= Checkbox("Faces", &m_show_faces);
			SameLine();
			dirty |= Checkbox("Edges", &m_show_edges);
			SameLine();
			dirty |= Checkbox("Verts", &m_show_verts);
			SameLine();
			SetCursorPosX(GetCursorPosX() + std::max(0.f, GetContentRegionAvail().x - 20));
			PushStyleColor(ImGuiCol_Button, ImVec4{0.6f, 0.3f, 0.3f, 1});
			if (Button("X", {-1, 0}) || m_try_kill) {
				OpenPopup("##close");
				m_try_kill = false;
			}
			if (IsItemHovered()) SetTooltip("Close model");
			PopStyleColor();
			if (BeginPopupModal("##close", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoDecoration)) {
				Text("Close model \"%s\" ?", name().c_str());
				SetHoveredTooltip("%s", make_name_tooltip().c_str());
				if (Button("Close")) {
					std::unique_lock lock(m_modelmtx, std::defer_lock);
					if (lock.try_lock()) {
						m_dead = true;
						invalidate_scene();
					} else {
						spawn_locked_notification();
					}
					CloseCurrentPopup();
				}
				SameLine();
				if (Button("Cancel")) {
					CloseCurrentPopup();
				}
				EndPopup();
			}
			Separator();
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
			if (dirty) invalidate_scene();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_selection() {
		using namespace ImGui;
		if (!m_model) return;
		if (Begin("Selection")) {
			PushID(this);
			draw_select_header(true);
			SetHoveredTooltip("%s", make_name_tooltip().c_str());
			SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);

			Text("%zd vertices, %zd triangles", m_model->mesh().n_vertices(), m_model->mesh().n_faces());
			if (m_decimated) Text(
				"Decimated %d vertices (%.1f%%) in %.3fs",
				m_src_dec_progress.completed_collapses,
				100.f * m_src_dec_progress.completed_collapses / float(m_src_dec_uparams.targetverts + m_src_dec_progress.completed_collapses),
				m_src_dec_progress.elapsed_time / std::chrono::duration<double>(1.0)
			);

			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			if (m_pending_save.valid()) {
				if (m_pending_save.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
					bool save_ok = false;
					std::string msg;
					try {
						save_ok = m_pending_save.get();
						msg = fmt::format("Export \"{}\" succeeded", m_fpath_save.u8string());
					} catch (std::exception &e) {
						std::cerr << "failed to save model: " << e.what() << std::endl;
						save_ok = false;
						msg = fmt::format("Export \"{}\" failed : {}", m_fpath_save.u8string(), e.what());
					} catch (...) {
						std::cerr << "failed to save model" << std::endl;
						save_ok = false;
						msg = fmt::format("Export \"{}\" failed", m_fpath_save.u8string());
					}
					ui_spawn([=, name=name(), nametooltip=make_name_tooltip(), msg=std::move(msg)](bool *p_open) {
						const char *title = save_ok ? "Export succeeded" : "Export failed";
						if (*p_open) OpenPopup(title);
						if (BeginPopupModal(title, p_open, ImGuiWindowFlags_AlwaysAutoResize)) {
							SelectHeader(true, name.c_str(), nametooltip.c_str());
							TextUnformatted(msg);
							if (Button("  OK  ")) {
								*p_open = false;
								CloseCurrentPopup();
							}
							EndPopup();
						}
					});
				} else {
					TextDisabled("Exporting %s...", m_fpath_save.filename().u8string().c_str());
				}
			}

			if (CollapsingHeader("Transform")) {
				auto pick_basis = [&](const char *label, int *basis) {
					SetNextItemWidth(GetTextLineHeight() * 3.5f);
					Combo(
						label, basis,
						[](void *data, int item, const char **out_text) {
							auto vs = reinterpret_cast<basis_vector *>(data);
							*out_text = vs[item].name;
							return true;
						}, m_basis_vectors, 6
					);
				};
				pick_basis("Right", &m_basis_right);
				SameLine();
				pick_basis("Up", &m_basis_up);
				SameLine();
				pick_basis("Back", &m_basis_back);
				if (Button("Reset##scale")) m_scale = m_model->unit_bound_scale() * 4;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderFloat("Scale", &m_scale, 0, 1000, "%.4f", 8);
				if (Button("Reset##translation")) m_translation = glm::vec3{0};
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				InputFloat3("Translation", value_ptr(m_translation));
				if (Button("Reset##roty")) m_rotation_euler_yxz.y = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Yaw", &m_rotation_euler_yxz.y, -180, 180);
				if (Button("Reset##rotx")) m_rotation_euler_yxz.x = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Pitch", &m_rotation_euler_yxz.x, -180, 180);
				if (Button("Reset##rotz")) m_rotation_euler_yxz.z = 0;
				SameLine();
				SetNextItemWidth(GetContentRegionAvail().x * 0.6f);
				SliderAngle("Roll", &m_rotation_euler_yxz.z, -180, 180);
			}

			if (CollapsingHeader("Rendering")) {
				bool r = false;
				if (RadioButton("Smooth Shading", !m_shade_flat)) {
					m_shade_flat = false;
					r = true;
				}
				SameLine();
				if (RadioButton("Flat Shading", m_shade_flat)) {
					m_shade_flat = true;
					r = true;
				}
				r |= Checkbox("Color Faces", &m_color_faces);
				SetHoveredTooltip("Apply color mode to faces");
				SameLine();
				r |= Checkbox("Color Verts", &m_color_verts);
				SetHoveredTooltip("Apply color mode to vertices");
				r |= Checkbox("Cull Faces", &m_cull_faces);
				// TODO culling 'back' edges not supported atm
				//SameLine();
				//r |= Checkbox("Cull Edges", &m_cull_edges);
				r |= SliderInt("Point Size", &m_vert_point_size, 1, 5);
				m_vert_point_size = std::max(1, m_vert_point_size);
				if (r) invalidate_scene();
			}

			Separator();

			if (combo_color_mode(m_disp_color_mode)) invalidate_saliency_vbo();

			Separator();

			std::vector<model_saliency_data> empty_saliency_outputs;
			auto &saliency_outputs = m_model ? m_model->saliency() : empty_saliency_outputs;

			if (Combo(
				"Saliency", &m_saliency_index,
				[](void *data, int item, const char **out_text) {
					auto ds = reinterpret_cast<model_saliency_data *>(data);
					const auto &so = ds[item];
					// TODO better?
					static std::string str;
					str = std::string(so);
					*out_text = str.c_str();
					return true;
				}, 
				saliency_outputs.data(), saliency_outputs.size()
					)) {
				invalidate_saliency_vbo();
			}

			if (m_disp_color_mode == model_color_mode::saliency_comparison && m_saliency_baseline_index < saliency_outputs.size()) {
				Text("Baseline: %s", saliency_outputs[m_saliency_baseline_index].str().c_str());
				if (m_saliency_index < saliency_outputs.size()) {
					auto &err = m_saliency_errors;
					Text("Errors: min=%.3f, max=%.3f, rmse=%.3f", err.min, err.max, err.rms);
					if (SliderFloat("Error Scale", &m_saliency_error_scale, 1, 100, "%.3f", 2)) {
						invalidate_saliency_vbo();
					}
				}
			}

			if (Button("Paste")) {
				auto &clip = saliency_clipboard();
				std::unique_lock lock(m_modelmtx, std::defer_lock);
				if (!lock.try_lock()) {
					spawn_locked_notification();
				} else if (clip.data.size() == m_model->mesh().n_vertices()) {
					model_saliency_data sd = std::move(clip.sd);
					m_model->mesh().add_property(sd.prop_saliency);
					m_model->mesh().property(sd.prop_saliency).data_vector() = std::move(clip.data);
					clip.sd = {};
					clip.data.clear();
					saliency_outputs.push_back(std::move(sd));
					invalidate_saliency_vbo();
					// references/iterators into saliency outputs are invalidated
				} else {
					OpenPopup("Paste Error##pasteerror");
				}
			}

			if (m_saliency_index < saliency_outputs.size()) {
				auto &salout = saliency_outputs[m_saliency_index];
				SameLine();
				if (Button("Copy")) {
					std::shared_lock lock(m_modelmtx, std::defer_lock);
					if (!lock.try_lock()) {
						spawn_locked_notification();
					} else {
						auto &clip = saliency_clipboard();
						clip.sd = salout;
						clip.sd.prop_saliency.reset();
						clip.data = m_model->mesh().property(salout.prop_saliency).data_vector();
					}
				}
				SameLine();
				if (Button("Remove")) OpenPopup("##remove");
				SameLine();
				if (Button("Rename")) OpenPopup("Rename Saliency##rename");
				SameLine();
				if (Button("Baseline")) m_saliency_baseline_index = m_saliency_index;
				if (SliderFloat("Gamma", &salout.gamma, 0, 2)) {
					invalidate_saliency_vbo();
				}
				Checkbox("Persistent", &salout.persistent);
				SetHoveredTooltip("Persistent properties will be preserved when decimating\nand will be exported by default");
				Separator();
				if (CollapsingHeader("Saliency Parameters")) {
					if (salout.uparams_known) {
						draw_saliency_params(salout.uparams);
						if (Button("Reload")) m_sal_uparams = salout.uparams;
					} else {
						TextDisabled("Parameters unknown");
					}
				}
				if (salout.filename.empty()) {
					if (CollapsingHeader("Saliency Progress")) {
						draw_saliency_progress(salout.progress);
					}
				}
				if (BeginPopupModal("Paste Error##pasteerror")) {
					auto &clip = saliency_clipboard();
					if (clip.data.size() == m_model->mesh().n_vertices()) {
						CloseCurrentPopup();
					} else {
						Text("Can't paste saliency for %d vertices into model with %d vertices", clip.data.size(), m_model->mesh().n_vertices());
						if (Button("OK")) CloseCurrentPopup();
					}
				}
				bool rename_open = true;
				SetNextWindowSize({300, 100}, ImGuiCond_Appearing);
				if (BeginPopupModal("Rename Saliency##rename", &rename_open)) {
					Text("New saliency property name");
					// TODO better?
					static char buf[1024]{};
					if (IsWindowAppearing()) {
						strncpy(buf, salout.dispname.c_str(), sizeof(buf));
						SetKeyboardFocusHere();
					}
					SetNextItemWidth(GetContentRegionAvail().x);
					const char *badchars = " \t@[]";
					if (InputText("", buf, sizeof(buf), ImGuiInputTextFlags_EnterReturnsTrue)) {
						if (std::string_view(buf).find_first_of(badchars) == std::string::npos) {
							salout.dispname = buf;
							CloseCurrentPopup();
						}
					}
					if (std::string_view(buf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of: %s", badchars);
					}
					EndPopup();
				}
				if (BeginPopupModal("##remove")) {
					Text("Remove saliency result?");
					if (Button("Remove")) {
						std::unique_lock lock(m_modelmtx, std::defer_lock);
						if (!lock.try_lock()) {
							spawn_locked_notification();
						} else {
							auto it = saliency_outputs.begin() + m_saliency_index;
							m_model->mesh().remove_property(it->prop_saliency);
							saliency_outputs.erase(it);
							// references/iterators into saliency outputs are invalidated
							// (so this section should come last)
							if (m_saliency_index >= saliency_outputs.size()) {
								m_saliency_index = std::max(0, m_saliency_index - 1);
							}
							invalidate_saliency_vbo();
						}
						CloseCurrentPopup();
					}
					SameLine();
					if (Button("Cancel")) CloseCurrentPopup();
					EndPopup();
				}
			} else {
				SameLine();
				ButtonDisabled("Copy");
				SameLine();
				ButtonDisabled("Remove");
				SameLine();
				ButtonDisabled("Rename");
				SameLine();
				ButtonDisabled("Baseline");
			}
			Separator();
			PopID();
		}
		End();
	}

	void ModelEntity::draw_window_export() {
		using namespace ImGui;
		if (!m_model) return;
		auto &saliency_outputs = m_model->saliency();
		if (m_try_export) {
			m_try_export = false;
			m_saliency_export_index = m_saliency_index;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				sd.should_export = sd.persistent;
			}
			OpenPopup("Export##export");
		}
		SetNextWindowSize({500, 400}, ImGuiCond_Appearing);
		bool export_window_open = true;
		if (BeginPopupModal("Export##export", &export_window_open)) {
			PushID(this);
			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			draw_select_header(true);
			// export path
			const auto &pathhint = m_fpath_save.empty() ? m_fpath_load : m_fpath_save;
			// TODO better?
			static char pathbuf[1024]{};
			if (IsWindowAppearing()) snprintf(pathbuf, sizeof(pathbuf), "%s", pathhint.u8string().c_str());
			auto &fpath_save_fut = ui_save_path(pathhint, Button("Browse"));
			if (fpath_save_fut.valid() && fpath_save_fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				auto s = fpath_save_fut.get().u8string();
				if (!s.empty()) snprintf(pathbuf, sizeof(pathbuf), "%s", s.c_str());
			}
			SameLine();
			InputText("Path", pathbuf, sizeof(pathbuf));
			// export color mode
			combo_color_mode(m_exp_color_mode);
			if (m_exp_color_mode == model_color_mode::saliency_comparison) {
				if (m_saliency_baseline_index < m_model->saliency().size()) {
					Text("Baseline: %s", m_model->saliency()[m_saliency_baseline_index].str().c_str());
				} else {
					TextColored(badcol, "Invalid baseline");
				}
			}
			Separator();
			Text("Saliency Properties");
			if (BeginChild("saliency", {0, -115}, false, ImGuiWindowFlags_AlwaysVerticalScrollbar)) {
				for (int i = 0; i < saliency_outputs.size(); i++) {
					PushID(i);
					auto &sd = saliency_outputs[i];
					Checkbox("##shouldexport", &sd.should_export);
					SetHoveredTooltip("Export this saliency property?");
					SameLine();
					if (RadioButton("##colorize", i == m_saliency_export_index)) m_saliency_export_index = i;
					SetHoveredTooltip("Colorize this saliency property?");
					const char *badchars = " \t@[]";
					char propbuf[256]{};
					strncpy(propbuf, sd.expname.c_str(), sizeof(propbuf));
					SameLine();
					SetNextItemWidth(150);
					if (InputText("##expname", propbuf, sizeof(propbuf))) {
						if (std::string_view(propbuf).find_first_of(badchars) == std::string::npos) {
							sd.expname = propbuf;
						}
					}
					SetHoveredTooltip("Export property name.\nDefault when empty is the current display name.");
					SameLine();
					if (sd.expname.empty()) {
						auto p0 = GetCursorPos();
						SetCursorPosX(p0.x - 150);
						TextDisabled("%s", sd.dispname.c_str());
						SameLine();
						SetCursorPosX(p0.x);
					}
					if (Button("<")) sd.expname = sd.dispname;
					SetHoveredTooltip("Copy display name to export name");
					SameLine();
					Text("%s", sd.str().c_str());
					SetHoveredTooltip(sd.uparams_known ? sd.uparams.str().c_str() : "Parameters unknown");
					if (std::string_view(propbuf).find_first_of(badchars) != std::string::npos) {
						TextColored(badcol, "Property name must not contain any of \"%s\"", badchars);
					}
					PopID();
				}
			}
			EndChild();
			Separator();
			// other options
			Checkbox("Binary", &m_save_binary);
			SetHoveredTooltip("Save file as binary if supported");
			SameLine();
			Checkbox("Original Vertex IDs", &m_save_original_vids);
			SetHoveredTooltip("Export extra property with original vertex ids from before decimation");
			TextDisabled("Supported formats are currently PLY and OBJ.\nColorization and property export are only supported for PLY.");
			// validate path and maybe export
			auto fpath = std::filesystem::u8path(pathbuf);
			auto stat = std::filesystem::status(fpath);
			bool cansave = !fpath.empty();
			const char *badchars = "\\/:*?\"<>|";
			if (fpath.is_relative()) {
				TextColored(badcol, "Path must be absolute");
				cansave = false;
			} else if (std::filesystem::is_directory(stat)) {
				TextColored(badcol, "Path is a directory");
				cansave = false;
			} else if (!std::filesystem::is_directory(fpath.parent_path())) {
				TextColored(badcol, "Directory does not exist");
				cansave = false;
			} else if (fpath.filename().u8string().find_first_of(badchars) != std::string::npos) {
				TextColored(badcol, "File name must not contain any of \"%s\"", badchars);
				cansave = false;
			} else if (std::filesystem::exists(stat)) {
				TextColored(badcol, "Path exists! Save will overwrite");
				//TextColored(badcol, u8"保存先が存在します！保存したら上書きします");
			}
			if (Button("Save") && cansave && !m_pending_save.valid()) {
				save(std::filesystem::absolute(fpath));
				CloseCurrentPopup();
			}
			if (IsItemHovered()) {
				if (m_pending_save.valid()) {
					SetTooltip("Export already in progress");
				} else if (!cansave) {
					SetTooltip("Can't save to this path");
				}
			}
			SameLine();
			if (Button("Cancel")) CloseCurrentPopup();
			PopID();
			EndPopup();
		}
	}

	void ModelEntity::draw_window_saliency() {
		using namespace ImGui;
		if (ui_saliency_window_open()) {
			if (Begin("Saliency")) {
				PushID(this);
				draw_select_header(true);
				m_sal_need_preview |= Checkbox("Preview", &m_sal_uparams.preview);
				SetHoveredTooltip("Enable interactive preview\nUse on small models only! (< ~500k vertices)\nCoarse subsampling will be activated.");
				SameLine();
				bool go = false;
				if (m_sal_uparams.preview || m_sal_future.valid()) {
					ButtonDisabled("GO", {-1, 0});
				} else {
					if (Button("GO", {-1, 0})) go = true;
				}
				Separator();
				// TODO easily obtain settings from other models
				m_sal_need_preview |= edit_saliency_params(m_sal_uparams);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				if (m_sal_need_preview && m_sal_uparams.preview) {
					// user params edited and preview enabled, try launch preview
					if (m_sal_future.valid()) {
						m_sal_progress.should_cancel = true;
					} else {
						go = true;
						m_sal_need_preview = false;
					}
				}
				if (go) {
					// launch calculation
					auto real_uparams = m_sal_uparams;
					if (real_uparams.preview) {
						real_uparams.subsample_auto = true;
						real_uparams.subsample_manual = false;
						real_uparams.samples_per_neighborhood = 10;
					}
					m_sal_progress = {};
					m_sal_progress.levels.resize(real_uparams.levels);
					saliency_async(real_uparams);
				}
				PopID();
			}
			End();
		}
	}

	void ModelEntity::draw_window_saliency_progress(bool selected) {
		using namespace ImGui;
		// check for saliency completion
		if (m_sal_future.valid()) {
			if (m_sal_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
				// getting the result will trigger cleanup (from its dtor)
				if (!m_sal_future.get()) {
					// cancelled
				}
			}
		}
		if (ui_saliency_window_open() && (saliency_state::idle < m_sal_progress.state && (m_sal_progress.state < saliency_state::done || selected))) {
			if (Begin("Saliency")) {
				PushID(this);
				draw_select_header(selected);
				draw_saliency_progress(m_sal_progress);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				PopID();
			}
			End();
		}
	}

	void ModelEntity::draw_window_decimation() {
		using namespace ImGui;
		if (ui_decimation_window_open()) {
			if (Begin("Decimation")) {
				FmtTextColored(extra_colors().bad, "-- Work in Progress --");
				if (m_decimated) FmtTextColored(extra_colors().bad, "Warning: model has already been decimated");
				auto *sd = selected_saliency();
				if (sd) FmtText("Saliency: {}", sd->str());
				bool go = false;
				if (sd && m_dec_uparams.use_saliency) {
					if (Button("GO with saliency", {-1, 0})) go = true;
				} else if (!m_dec_uparams.use_saliency) {
					if (Button("GO without saliency", {-1, 0})) go = true;
				} else {
					TextDisabled("Select a saliency result or uncheck 'use saliency'");
				}
				Separator();
				edit_decimate_params(m_dec_uparams);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				if (go) {
					// launch
					decimate_async(m_dec_uparams);
				}
			}
			End();
		}
	}

	void ModelEntity::draw_window_decimation_progress(bool selected) {
		using namespace ImGui;
		if (m_decimated && ui_decimation_window_open() && (selected || m_src_dec_progress.state < decimation_state::done)) {
			if (Begin("Decimation")) {
				PushID(this);
				draw_select_header(selected);
				draw_decimate_progress(m_src_dec_progress);
				Separator();
				SetCursorPosY(GetCursorPosY() + GetStyle().ItemSpacing.y);
				PopID();
			}
			End();
		}
	}

	bool ModelEntity::draw_select_header(bool selected) {
		using namespace ImGui;
		bool r = false;
		if (SelectHeader(selected, name().c_str(), make_name_tooltip().c_str())) {
			auto &sel = ui_selection();
			sel.select_entity = id();
			sel.select_vertex = -1;
			r = true;
		}
		return r;
	}

	void ModelEntity::spawn_locked_notification() const {
		using namespace ImGui;
		ui_spawn([nametooltip=make_name_tooltip(), name=name()](bool *p_open) {
			if (*p_open) OpenPopup("Model in use");
			if (BeginPopupModal("Model in use", p_open, ImGuiWindowFlags_AlwaysAutoResize)) {
				SelectHeader(true, name.c_str(), nametooltip.c_str());
				Text("The operation cannot be performed because\nthe model is currently in use.");
				if (Button("  OK  ")) {
					*p_open = false;
					CloseCurrentPopup();
				}
				EndPopup();
			}
		});
	}

	std::string ModelEntity::make_name_tooltip() const {
		auto s = m_fpath_load.u8string();
		if (m_decimated) {
			s += "\n";
			s += m_src_dec_uparams.str();
		}
		return s;
	}

	void ModelEntity::invalidate_saliency_vbo() {
		m_saliency_vbo_dirty = true;
		invalidate_scene();
	}

	void ModelEntity::pre_draw() {
		auto &sel = ui_selection();
		const bool selected = sel.select_entity == id();
		if (!dead() && selected) ui_select_model(this);
		if (!dead() && sel.hover_entity == id()) ui_hover_model(this);
		if (!selected) return;
		draw_window_saliency();
		draw_window_decimation();
	}

	void ModelEntity::draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) {
		const auto xform0 = transform();
		if (m_pending_load.valid() && m_pending_load.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
			try {
				// is this the best place to do this?
				m_model = m_pending_load.get();
			} catch (std::exception &e) {
				std::cerr << "failed to load model: " << e.what() << std::endl;
			} catch (...) {
				std::cerr << "failed to load model" << std::endl;
			}
			// TODO improve load failure handling (destroy entity, show popup?)
			if (m_model) {
				// gl stuff has to happen on main thread (and dont catch exceptions from it)
				m_model->update_vaos();
				m_model->update_vbos();
				m_scale = m_model->unit_bound_scale() * 4;
				invalidate_saliency_vbo();
			}
		}
		auto &sel = ui_selection();
		const bool selected = sel.select_entity == id();
		draw_window_models(selected);
		draw_window_saliency_progress(selected);
		draw_window_decimation_progress(selected);
		if (m_model) {
			if (selected) draw_window_selection();
			if (selected) draw_window_export();
			if (draw_scene) {
				update_vbo();
				GLuint tex = 0;
				// determine color map to apply in shader
				auto &saliency_outputs = m_model->saliency();
				int color_map = 0;
				const bool sal_valid = m_saliency_index < saliency_outputs.size();
				if (m_disp_color_mode == model_color_mode::saliency && sal_valid) color_map = 3;
				if (m_disp_color_mode == model_color_mode::saliency_comparison && sal_valid && m_saliency_baseline_index < saliency_outputs.size()) color_map = 4;
				if (m_disp_color_mode == model_color_mode::vcolor && m_model->prop_vcolor_original().is_valid()) color_map = 1;
				if (m_disp_color_mode == model_color_mode::doncurv) color_map = 3;
				if (m_disp_color_mode == model_color_mode::normal) color_map = 7;
				if (m_disp_color_mode == model_color_mode::uv && m_model->has_texcoords2d()) color_map = 5;
				if (m_disp_color_mode == model_color_mode::checkerboard && m_model->has_texcoords2d()) {
					color_map = 6;
					tex = checkerboard_texture();
				}
				// prepare to draw
				model_draw_params params;
				params.sel = sel;
				params.entity_id = id();
				// faces
				params.polymode = GL_FILL;
				params.shade_flat = m_shade_flat;
				params.color = {0.6f, 0.6f, 0.5f, 1};
				params.vert_color_map = m_color_faces ? color_map : 0;
				glCullFace(GL_BACK);
				if (m_cull_faces) glEnable(GL_CULL_FACE);
				glColorMaski(1, GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_faces) {
					glActiveTexture(GL_TEXTURE0);
					if (tex) glBindTexture(GL_TEXTURE_2D, tex);
					m_model->draw(view * transform(), proj, zfar, params);
					glBindTexture(GL_TEXTURE_2D, 0);
				}
				// edges
				params.polymode = GL_LINE;
				params.shade_flat = false;
				params.shading = 0;
				params.color = {0.03f, 0.03f, 0.03f, 1};
				params.vert_color_map = 0;
				// culling not supported atm
				//set_cull_faces(m_cull_edges);
				glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_edges) m_model->draw(view * transform(), proj, zfar, params);
				// boundary edges
				params.boundaries = true;
				params.shade_flat = false;
				params.shading = 0;
				params.color = {1, 0, 0, 1};
				params.vert_color_map = 0;
				// culling not supported atm
				//set_cull_faces(m_cull_edges);
				glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
				if (m_show_faces) m_model->draw(view * transform(), proj, zfar, params);
				// verts
				params.boundaries = false;
				params.polymode = GL_POINT;
				params.color = {0.5f, 0, 0, 1};
				params.vert_color_map = m_color_verts ? color_map : 0;
				params.sel.hover_entity = -1;
				params.sel.select_entity = -1;
				//params.show_samples = sal_valid && (m_color_mode == color_mode::saliency || m_color_mode == color_mode::saliency_comparison);
				glDisable(GL_CULL_FACE);
				glPointSize(m_vert_point_size);
				glColorMaski(1, GL_TRUE, GL_TRUE, GL_FALSE, GL_FALSE);
				if (m_show_verts) m_model->draw(view * transform(), proj, zfar, params);

			}
		}
		const auto xform1 = transform();
		if (xform0 != xform1) invalidate_scene();
	}

	void ModelEntity::saliency_async(const saliency_user_params &uparams0) {
		if (!m_model) return;
		// can't run multiple saliency calculations for one model at once (because we only track one future).
		// otherwise, multiple saliency calculations can run simultaneously (although suboptimally).
		if (m_sal_future.valid()) return;
		std::unique_lock lock(m_modelmtx, std::defer_lock);
		if (!lock.try_lock()) {
			spawn_locked_notification();
			return;
		}
		saliency_user_params uparams = uparams0;
		saliency_mesh_params mparams;
		m_model->init_saliency_params(mparams, uparams);
		// cleanup will be run when the result is received from the future
		// can use shared lock for actual computation because only property contents are being used
		// HACK make_shared is an ugly workaround for std::function needing to be copyable
		lock.unlock();
		mparams.cleanup = [=, pprogress=&m_sal_progress, slock=std::make_shared<std::shared_lock<std::shared_mutex>>(m_modelmtx)](bool r) mutable noexcept {
			// now need to upgrade to unique lock
			slock->unlock();
			std::unique_lock lock(m_modelmtx);
			if (!m_model) return;
			auto &saliency_outputs = m_model->saliency();
			m_model->cleanup_saliency_params(mparams, r);
			if (r) {
				// if preview, remove any previous preview results
				saliency_outputs.erase(std::remove_if(saliency_outputs.begin(), saliency_outputs.end(), [](const auto &sd) { return sd.uparams.preview; }), saliency_outputs.end());
				// save user params, progress output and actual saliency mesh property
				model_saliency_data sd;
				sd.uparams = uparams;
				sd.progress = *pprogress;
				sd.prop_saliency = mparams.prop_saliency;
				sd.prop_sampled = mparams.prop_sampled;
				sd.uparams_known = true;
				saliency_outputs.push_back(sd);
				// give focus to this result
				m_saliency_index = saliency_outputs.size() - 1;
				invalidate_saliency_vbo();
			}
		};
		m_sal_future = green::compute_saliency_async(mparams, uparams, m_sal_progress);
	}

	void ModelEntity::decimate_async(const decimate_user_params &uparams) {
		// this function isn't const because it needs to lock the mutex - mutable?
		std::shared_lock lock1(m_modelmtx, std::defer_lock);
		if (!lock1.try_lock()) {
			spawn_locked_notification();
			return;
		}
		auto e = std::make_unique<ModelEntity>();
		std::unique_lock lock2(e->m_modelmtx);
		e->m_fpath_load = m_fpath_load;
		e->m_decimated = true;
		e->m_src_dec_uparams = uparams;
		e->m_basis_right = m_basis_right;
		e->m_basis_up = m_basis_up;
		e->m_basis_back = m_basis_back;
		e->m_rotation_euler_yxz = m_rotation_euler_yxz;
		e->m_cull_faces = m_cull_faces;
		e->m_cull_edges = m_cull_edges;
		e->m_pending_load = std::async([=, e=e.get(), lock1=std::move(lock1), lock2=std::move(lock2)]() mutable {
			// prepare mesh copy to decimate
			// this is slow enough on big models to be async
			auto &saliency_outputs = m_model->saliency();
			std::vector<model_saliency_data> sdv;
			saliency_prop_t srcprop;
			for (int i = 0; i < saliency_outputs.size(); i++) {
				auto &sd = saliency_outputs[i];
				if (i == m_saliency_index || sd.persistent) {
					srcprop = sd.prop_saliency;
					sdv.push_back(sd);
				} else if (sd.persistent) {
					sdv.push_back(sd);
				}
			}
			auto m = m_model->prepare_decimate(srcprop, sdv);
			// unlock source model
			lock1.unlock();
			// now decimate
			if (!m.decimate(e->m_src_dec_uparams, e->m_src_dec_progress)) throw std::runtime_error("decimation was cancelled");
			return std::make_unique<Model>(std::move(m));
		});
		spawn_entity(std::move(e));
	}

	ModelEntity::~ModelEntity() {
		m_sal_progress.should_cancel = true;
		if (m_sal_future.valid()) m_sal_future.wait();
	}

}
