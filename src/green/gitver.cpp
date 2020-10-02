
#include "gitver.hpp"

#include <string_view>

#include "gitverdetail.h"

namespace green {

	const char * git_describe() {
		return GREEN_GIT_DESCRIBE;
	}

	const char * git_revision() {
		return GREEN_GIT_REV;
	}

	const char * git_timestamp() {
		return GREEN_GIT_TIMESTAMP;
	}

	bool git_has_changes() {
		return GREEN_GIT_HAS_CHANGES;
	}

	bool git_is_release() {
		// proper releases just look like v0.2
		// anything else looks like v0.2-N-gHASH
		return std::string_view(git_describe()).find("-g") == std::string::npos;
	}

}
