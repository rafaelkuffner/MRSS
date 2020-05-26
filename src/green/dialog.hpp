
#pragma once

#ifndef GREEN_DIALOG_HPP
#define GREEN_DIALOG_HPP

#include <vector>
#include <filesystem>
#include <future>

#include <GLFW/glfw3.h>

namespace green {

	std::future<std::vector<std::filesystem::path>> open_file_dialog(GLFWwindow *win, const std::filesystem::path &hint, bool multi);

	std::future<std::filesystem::path> save_file_dialog(GLFWwindow *win, const std::filesystem::path &hint);

}

#endif
