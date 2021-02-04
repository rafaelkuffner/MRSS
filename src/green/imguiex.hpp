
#pragma once

#ifndef GREEN_IMGUIEX_HPP
#define GREEN_IMGUIEX_HPP

#include <fmt/format.h>

#include <imgui.h>

#include "uilocale.hpp"

namespace ImGui {

	using green::uilocale;
	using green::uistring;

	struct extra_colors_t {
		ImVec4 bad;
		ImVec4 bad_hov;
		ImVec4 bad_bg;
	};

	extra_colors_t & extra_colors();

	bool ButtonDisabled(const char* label, const ImVec2& size_arg = ImVec2{0, 0});

	template <typename S, typename ...Ts>
	inline void FmtText(const S &s, const Ts &...args) {
		::fmt::memory_buffer buf;
		::fmt::vformat_to(buf, s, ::fmt::make_format_args(args...));
		TextUnformatted(buf.begin(), buf.end());
	}

	template <typename ...Ts>
	inline void FmtText(const green::uilocale &loc, green::uistring s, const Ts &...args) {
		::fmt::memory_buffer buf;
		loc[s].format_to(buf, args...);
		TextUnformatted(buf.begin(), buf.end());
	}

	template <typename S, typename ...Ts>
	inline void FmtTextColored(ImVec4 col, const S &s, const Ts &...args) {
		::fmt::memory_buffer buf;
		::fmt::vformat_to(buf, s, ::fmt::make_format_args(args...));
		PushStyleColor(ImGuiCol_Text, col);
		TextUnformatted(buf.begin(), buf.end());
		PopStyleColor();
	}

	template <typename ...Ts>
	inline void FmtTextColored(ImVec4 col, const green::uilocale &loc, green::uistring s, const Ts &...args) {
		::fmt::memory_buffer buf;
		loc[s].format_to(buf, args...);
		PushStyleColor(ImGuiCol_Text, col);
		TextUnformatted(buf.begin(), buf.end());
		PopStyleColor();
	}

	template <typename S, typename ...Ts>
	inline void FmtTooltip(const S &s, const Ts &...args) {
		::fmt::memory_buffer buf;
		::fmt::vformat_to(buf, s, ::fmt::make_format_args(args...));
		BeginTooltip();
		TextUnformatted(buf.begin(), buf.end());
		EndTooltip();
	}

	template <typename ...Ts>
	inline void FmtTooltip(const green::uilocale &loc, green::uistring s, const Ts &...args) {
		::fmt::memory_buffer buf;
		loc[s].format_to(buf, args...);
		BeginTooltip();
		TextUnformatted(buf.begin(), buf.end());
		EndTooltip();
	}

	template <typename S, typename ...Ts>
	inline void FmtHoveredTooltip(const S &s, const Ts &...args) {
		if (IsItemHovered()) FmtTooltip(s, args...);
	}

	template <typename ...Ts>
	inline void SetHoveredTooltip(const char *fmt, const Ts &...args) {
		if (IsItemHovered()) SetTooltip(fmt, args...);
	}

	template <typename ...Ts>
	inline void SetHoveredTooltip(const uilocale &loc, uistring fmt, const Ts &...args) {
		if (IsItemHovered()) SetTooltip(loc[fmt].c_str(), args...);
	}

	template <typename T>
	struct param_widgets {
		const uilocale *loc = nullptr;
		T *defparams = nullptr;
		T *realparams = nullptr;
		float label_width = 80;

		void draw_label(const char *label) {
			const auto p0 = GetCursorPos();
			const auto ap0 = GetCursorScreenPos();
			PushClipRect(ap0, {ap0.x + label_width, ap0.y + GetTextLineHeightWithSpacing()}, true);
			TextUnformatted(label);
			PopClipRect();
			SameLine(p0.x + label_width);
		}

		void hovered_tooltip(uistring tooltip) {
			if (!tooltip || !IsItemHovered()) return;
			FmtTooltip(*loc, tooltip);
		}

		bool slider_any(const char *label, int *v, int v_min, int v_max, const char *format = "%d") {
			draw_label(label);
			SetNextItemWidth(GetContentRegionAvail().x);
			return SliderInt(label, v, v_min, v_max, format);
		}

		bool slider_any(const char *label, float *v, float v_min, float v_max, const char *format = "%.3f", float power = 1.f) {
			draw_label(label);
			SetNextItemWidth(GetContentRegionAvail().x);
			return SliderFloat("", v, v_min, v_max, format, power);
		}

		template <typename U, typename ...Args>
		bool slider(uistring label, uistring tooltip, U T::*param, U vmin, U vmax, const Args &...args) {
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			auto &val = realparams->*param;
			const bool oor = !(val >= vmin && val <= vmax);
			if (oor) {
				PushStyleColor(ImGuiCol_FrameBg, extra_colors().bad_bg);
				PushStyleColor(ImGuiCol_FrameBgHovered, extra_colors().bad_hov);
				PushStyleColor(ImGuiCol_SliderGrab, extra_colors().bad);
				// can't be out-of-range and active at the same time
			}
			bool r = slider_any((*loc)[label].c_str(), &val, vmin, vmax, args...);
			hovered_tooltip(tooltip);
			if (oor) PopStyleColor(3);
			PopID();
			return r;
		}

		bool checkbox(uistring label, uistring tooltip, bool T::*param) {
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			bool r = Checkbox((*loc)[label].c_str(), &(realparams->*param));
			hovered_tooltip(tooltip);
			PopID();
			return r;
		}

		bool inputint(uistring label, uistring tooltip, int T::*param, int step = 1, int step_fast = 100) {
			PushID(label);
			if (Button("Reset")) realparams->*param = defparams->*param;
			SameLine();
			bool r = InputInt((*loc)[label].c_str(), &(realparams->*param), step, step_fast);
			hovered_tooltip(tooltip);
			PopID();
			return r;
		}
	};

	template <typename T>
	param_widgets(const uilocale *, T *, T *) -> param_widgets<T>;

}

#endif
