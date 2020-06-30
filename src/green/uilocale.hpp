
#pragma once

#ifndef GREEN_UILOCALE_HPP
#define GREEN_UILOCALE_HPP

#include <string>
#include <string_view>

namespace green {

	namespace uistrings {
		enum uistring {
			// common
			help_input,
			help_output,
			help_threads,
			help_cli_opts_common,
			help_cli_version,
			help_cli_help,
			help_cli_ascii,
			help_cli_gui,
			param_threads,
			// saliency
			help_cli_opts_saliency,
			help_sal_go,
			help_sal_area,
			help_sal_levels,
			help_sal_normpower,
			help_sal_curvweight,
			help_sal_noisefilter,
			help_sal_noiseheight,
			help_sal_samplespern,
			help_sal_subsample,
			help_sal_full,
			param_sal_area,
			param_sal_levels,
			param_sal_normpower,
			param_sal_curvweight,
			param_sal_noisefilter,
			param_sal_noiseheight,
			param_sal_samplespern,
			param_sal_subsample,
			// decimation
			help_cli_opts_decimate,
			help_dec_go,
			help_dec_targetverts,
			help_dec_weight,
			help_dec_power,
			help_dec_bins,
			help_dec_propname,
			help_dec_usesaliency,
			help_dec_nosaliency,
			param_dec_usesaliency,
			param_dec_targetverts,
			param_dec_weight,
			param_dec_power,
			param_dec_bins,

			// ---
			_count
		};
	}

	using uistrings::uistring;

	class uistring_impl : public std::string_view {
	public:
		// init to empty string; never null
		constexpr uistring_impl() : std::string_view("") {}

		// string literals only!
		// null termination is required, and no copy is made.
		template <size_t N>
		constexpr explicit uistring_impl(const char(&str)[N]) noexcept : std::string_view(str) {}

		// string literals only!
		// null termination is required, and no copy is made.
		template <size_t N>
		constexpr uistring_impl & operator=(const char(&str)[N]) noexcept {
			std::string_view::operator=(std::string_view(str));
			return *this;
		}

		constexpr const char * c_str() const noexcept {
			// this is safe because we require the string view to be constructed from a literal
			return data();
		}

		std::string clone() const {
			return std::string{*this};
		}
	};

	class uilocale {
	private:
		uistring_impl m_strings[uistrings::_count];

	public:
		constexpr uilocale() {}

		constexpr const uistring_impl & operator[](uistring s) const noexcept {
			return m_strings[s];
		}

		constexpr uistring_impl & operator[](uistring s) noexcept {
			return m_strings[s];
		}

		constexpr bool any_empty() const noexcept {
			for (auto &s : m_strings) {
				if (s.empty()) return true;
			}
			return false;
		}

	};

	const uilocale & uilocale_en();
	const uilocale & uilocale_jp();

}

#endif