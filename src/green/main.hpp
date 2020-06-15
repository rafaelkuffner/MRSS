
#pragma once

#ifndef GREEN_MAIN_HPP
#define GREEN_MAIN_HPP

#include <future>
#include <filesystem>
#include <functional>

#include "entity.hpp"
#include "model.hpp"
#include "saliency.hpp"

namespace green {

	entity_selection & ui_selection();

	saliency_user_params & ui_saliency_user_params();

	const saliency_progress & ui_saliency_progress();

	void ui_select_model(ModelEntity *);

	ModelEntity * ui_select_model();

	void ui_hover_model(ModelEntity *);

	ModelEntity * ui_hover_model();

	std::future<std::filesystem::path> & ui_save_path(const std::filesystem::path &hint, bool prompt);

	// spawn persistent ui component. persists until *p_open is set to false.
	void ui_spawn(std::function<void(bool *p_open)>);

}

namespace ImGui {

	void draw_saliency_progress(const green::saliency_progress &progress);

	bool edit_saliency_params(green::saliency_user_params &uparams);

	void draw_saliency_params(const green::saliency_user_params &uparams);

}

#endif
