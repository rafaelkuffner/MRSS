
// TODO implement simple dialogs with imgui
#ifdef GREEN_DIALOG_SIMPLE

#include "dialog.hpp"

namespace stdfs = std::filesystem;

namespace green {

	std::future<std::vector<std::filesystem::path>> open_file_dialog(GLFWwindow *win, const std::filesystem::path &hint, bool multi) {
		return {};
	}

	std::future<std::filesystem::path> save_file_dialog(GLFWwindow *win, const std::filesystem::path &hint) {
		return {};
	}

}

#endif
