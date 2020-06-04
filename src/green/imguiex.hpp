
#pragma once

#ifndef GREEN_IMGUIEX_HPP
#define GREEN_IMGUIEX_HPP

#include <imgui.h>

namespace ImGui {

	bool ButtonDisabled(const char* label, const ImVec2& size_arg = ImVec2{0, 0});

	template <typename ...Ts>
	void SetHoveredTooltip(const char *fmt, const Ts &...args) {
		if (IsItemHovered()) SetTooltip(fmt, args...);
	}

}

#endif
