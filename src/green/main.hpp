
#pragma once

#ifndef GREEN_MAIN_HPP
#define GREEN_MAIN_HPP

#include <future>
#include <filesystem>

#include "entity.hpp"
#include "model.hpp"
#include "saliency.hpp"

namespace green {

	entity_selection & ui_selection();

	saliency_user_params & ui_saliency_user_params();

	const saliency_progress & ui_saliency_progress();

	void ui_current_model(ModelEntity *);

	ModelEntity * ui_current_model();

	std::future<std::filesystem::path> & ui_save_path(const std::filesystem::path &hint, bool prompt);

}

namespace ImGui {

	void draw_saliency_progress(const green::saliency_progress &progress);

	bool edit_saliency_params(green::saliency_user_params &uparams);

	void draw_saliency_params(const green::saliency_user_params &uparams);

}

#endif
