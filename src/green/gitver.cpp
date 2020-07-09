
#include "gitver.hpp"

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

}
