#pragma once

#ifndef GREEN_MODEL_ENTITY_HPP
#define GREEN_MODEL_ENTITY_HPP

#include "model.hpp"

namespace green {

	class ModelEntity : public Entity {
	private:
		// unique: adding/removing properties, vertices etc
		// shared: reading/writing contents of existing properties
		std::shared_mutex m_modelmtx;
		std::unique_ptr<Model> m_model;

		// loading
		std::filesystem::path m_fpath_load;
		std::future<std::unique_ptr<Model>> m_pending_load;
		decimate_user_params m_src_dec_uparams;
		decimate_progress m_src_dec_progress;
		bool m_decimated = false;

		// computation
		int m_sal_preset = 0;
		saliency_user_params m_sal_uparams_custom;
		saliency_user_params m_sal_uparams_globality;
		saliency_progress m_sal_progress;
		std::future<saliency_result> m_sal_future;
		decimate_user_params m_dec_uparams;
		bool m_sal_want_preview = false;
		bool m_sal_need_preview = false;

		// saving
		std::filesystem::path m_fpath_save;
		model_color_mode m_exp_color_mode = model_color_mode::saliency;
		bool m_save_binary = true;
		bool m_save_original_vids = false;
		std::future<bool> m_pending_save;

		// selected saliency results
		int m_saliency_index = 0;
		int m_saliency_export_index = 0;
		int m_saliency_baseline_index = 0;
		model_saliency_errors m_saliency_errors;

		// transform
		int m_basis_right = 0, m_basis_up = 2, m_basis_back = 4;
		float m_scale = 1;
		glm::vec3 m_translation{0};
		// TODO store as quat, only use euler for editing?
		glm::vec3 m_rotation_euler_yxz{0};

		// rendering options
		model_color_mode m_disp_color_mode = model_color_mode::saliency;
		model_color_map m_disp_color_map = model_color_map::uniform;
		float m_saliency_error_scale = 1;
		float m_dec_err_gamma = 0.35f;
		int m_vert_point_size = 3;
		bool m_color_faces = true;
		bool m_color_verts = false;
		bool m_show_faces = true;
		bool m_show_edges = false;
		bool m_show_verts = false;
		bool m_cull_faces = true;
		bool m_cull_edges = true;
		bool m_shade_flat = false;

		// action/status
		bool m_try_export = false;
		bool m_try_kill = false;
		bool m_saliency_vbo_dirty = false;
		bool m_dead = false;

		struct basis_vector {
			glm::vec3 v;
			const char *name;
		};

		// TODO make this static
		basis_vector m_basis_vectors[6]{
			{{+1, 0, 0}, "+X"},
			{{-1, 0, 0}, "-X"},
			{{0, +1, 0}, "+Y"},
			{{0, -1, 0}, "-Y"},
			{{0, 0, +1}, "+Z"},
			{{0, 0, -1}, "-Z"}
		};

		void update_vbo();
		void draw_window_models(bool selected);
		void draw_window_selection();
		void draw_window_export();
		void draw_window_saliency();
		void draw_window_saliency_progress(bool selected);
		void draw_window_decimation();
		void draw_window_decimation_progress(bool selected);
		bool draw_select_header(bool selected);
		void spawn_locked_notification() const;
		std::string make_name_tooltip() const;
		void invalidate_saliency_vbo();
		void saliency_async(const saliency_user_params &uparams);
		void decimate_async(const decimate_user_params &uparams);

	public:
		ModelEntity();

		const Model * model() const {
			return m_model.get();
		}

		Model * model() {
			return m_model.get();
		}

		void load(const std::filesystem::path &fpath);

		void save(const std::filesystem::path &fpath);

		void load(const std::filesystem::path &fpath, Model m);

		void load(const std::filesystem::path &fpath, Model m, const decimate_user_params &dec_uparams, const decimate_progress &dec_progress);

		void try_export() {
			m_try_export = true;
		}

		bool show_entity_outline() const {
			return m_show_faces;
		}

		bool decimated() const {
			return m_decimated;
		}

		const model_saliency_data * selected_saliency() const {
			if (!m_model) return nullptr;
			if (m_saliency_index >= m_model->saliency().size()) return nullptr;
			return &m_model->saliency()[m_saliency_index];
		}

		virtual std::string name() const override {
			std::string s = m_fpath_load.filename().u8string();
			if (m_decimated) {
				s += " [dec ";
				s += m_src_dec_uparams.str();
				s += "]";
			}
			return s;
		}

		virtual void move_by(const glm::vec3 &d) override;

		virtual void try_kill() override {
			m_try_kill = true;
		}

		virtual bool dead() const override {
			return m_dead;
		}

		// basis transform from world space to model space
		glm::mat3 basis() const {
			return glm::mat3(m_basis_vectors[m_basis_right].v, m_basis_vectors[m_basis_up].v, m_basis_vectors[m_basis_back].v);
		}

		virtual glm::mat4 transform() const override;

		virtual void pre_draw() override;

		virtual void draw(const glm::mat4 &view, const glm::mat4 &proj, float zfar, bool draw_scene) override;

		virtual ~ModelEntity();
	};

}

#endif
