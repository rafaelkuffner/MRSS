
#pragma once

#ifndef GREEN_CONFIG_HPP
#define GREEN_CONFIG_HPP

#include "saliency.hpp"

#include <filesystem>

namespace green {

	bool load_config(const std::filesystem::path &);

	bool save_presets(const std::filesystem::path &);

	bool save_saliency_presets(const std::filesystem::path &, const std::vector<saliency_preset> &);

}

#endif
