
#include "imguiex.hpp"

#include <imgui_internal.h>

namespace ImGui {

	bool ButtonDisabled(const char *label, const ImVec2 &size_arg) {
		PushStyleColor(ImGuiCol_Text, GetStyle().Colors[ImGuiCol_TextDisabled]);
		bool r = ButtonEx(label, size_arg, ImGuiButtonFlags_Disabled);
		PopStyleColor();
		return r;
	}

}
