
#pragma once

#ifndef GREEN_IMGUIEX_HPP
#define GREEN_IMGUIEX_HPP

#include <imgui.h>

namespace ImGui {

	bool ButtonDisabled(const char* label, const ImVec2& size_arg = ImVec2{0, 0});

	template <typename ...Ts>
	inline void SetHoveredTooltip(const char *fmt, const Ts &...args) {
		if (IsItemHovered()) SetTooltip(fmt, args...);
	}

	inline bool SliderAny(const char *label, int *v, int v_min, int v_max, const char *format = "%d") {
		return SliderInt(label, v, v_min, v_max, format);
	}

	inline bool SliderAny(const char *label, float *v, float v_min, float v_max, const char *format = "%.3f", float power = 1.f) {
		return SliderFloat(label, v, v_min, v_max, format, power);
	}

	template <typename T>
	struct param_widgets {
		T *defparams = nullptr;
		T *realparams = nullptr;

		template <typename U, typename ...Args>
		bool slider(const char *label, U T::*param, U vmin, U vmax, const Args &...args) {
			const ImVec4 badcol{0.9f, 0.4f, 0.4f, 1};
			const ImVec4 badcolhov{badcol.x, badcol.y, badcol.z, GetStyle().Colors[ImGuiCol_FrameBg].w};
			const ImVec4 badcolbg{badcol.x * 0.7f, badcol.y * 0.7f, badcol.z * 0.7f, badcolhov.w};
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			auto &val = realparams->*param;
			const bool oor = !(val >= vmin && val <= vmax);
			if (oor) {
				PushStyleColor(ImGuiCol_FrameBg, badcolbg);
				PushStyleColor(ImGuiCol_FrameBgHovered, badcolhov);
				PushStyleColor(ImGuiCol_SliderGrab, badcol);
				// can't be out-of-range and active at the same time
			}
			bool r = SliderAny(label, &val, vmin, vmax, args...);
			if (oor) PopStyleColor(3);
			PopID();
			return r;
		}

		bool checkbox(const char *label, bool T::*param) {
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			bool r = Checkbox(label, &(realparams->*param));
			PopID();
			return r;
		}

		bool inputint(const char *label, int T::*param, int step = 1, int step_fast = 100) {
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			bool r = InputInt(label, &(realparams->*param), step, step_fast);
			PopID();
			return r;
		}
	};

	template <typename T>
	param_widgets(T *, T *) -> param_widgets<T>;

}

#endif
