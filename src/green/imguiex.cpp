
#include "imguiex.hpp"

#include <imgui_internal.h>

namespace {

	ImGui::extra_colors_t init_extra_colors() {
		using namespace ImGui;
		extra_colors_t c;
		c.bad = {0.9f, 0.4f, 0.4f, 1};
		c.bad_hov = {c.bad.x, c.bad.y, c.bad.z, GetStyle().Colors[ImGuiCol_FrameBg].w};
		c.bad_bg = {c.bad.x * 0.7f, c.bad.y * 0.7f, c.bad.z * 0.7f, c.bad_hov.w};
		return c;
	}

}

namespace ImGui {

	extra_colors_t & extra_colors() {
		static extra_colors_t c = init_extra_colors();
		return c;
	}

	bool ButtonDisabled(const char *label, const ImVec2 &size_arg) {
		PushStyleColor(ImGuiCol_Text, GetStyle().Colors[ImGuiCol_TextDisabled]);
		bool r = ButtonEx(label, size_arg, ImGuiButtonFlags_Disabled);
		PopStyleColor();
		return r;
	}

}
