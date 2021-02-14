
// use this header instead of getverdetail.h to improve build times

#pragma once

#ifndef GREEN_GITVER_HPP
#define GREEN_GITVER_HPP

namespace green {

	const char * git_describe();

	const char * git_branch();

	const char * git_revision();

	const char * git_timestamp();

	bool git_has_changes();

	bool git_is_release();

}

#endif
