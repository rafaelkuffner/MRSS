
#pragma once

#include "uilocale.hpp"
#include "uilocale_en.hpp"

namespace green {

	constexpr uilocale uilocale_init_jp() {
		using namespace uistrings;
		// default to english if not specified
		uilocale loc = uilocale_init_en();
		loc[help_cli_version] = "バージョンを見せます";
		// TODO obviously
		return loc;
	}

	const uilocale & uilocale_jp() {
		static constexpr uilocale loc = uilocale_init_jp();
		static_assert(!loc.any_empty(), "empty string in localization");
		return loc;
	}

}
