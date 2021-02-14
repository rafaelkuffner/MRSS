
cmake_minimum_required(VERSION 3.1)

message(STATUS "Retrieving git revision info")

set(git_rev unknown)
set(git_branch unknown)
set(git_describe unknown)
set(git_timestamp unknown)
set(git_status unknown)
set(git_has_changes 0)

find_package(Git)

set(descfile "${dst}.rev")
set(prevdesc "")
if(EXISTS "${descfile}")
	file(READ "${descfile}" prevdesc)
endif()

if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
		OUTPUT_VARIABLE git_rev
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
		OUTPUT_VARIABLE git_branch
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} describe HEAD
		OUTPUT_VARIABLE git_describe
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd HEAD
		OUTPUT_VARIABLE git_timestamp
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} --no-optional-locks status --porcelain
		OUTPUT_VARIABLE git_status
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	if(NOT "${git_status}" STREQUAL "")
		set(git_has_changes 1)
	endif()
	message(STATUS "Git rev: " ${git_rev})
	message(STATUS "Git branch: " ${git_branch})
	message(STATUS "Git describe: " ${git_describe})
	message(STATUS "Git timestamp: " ${git_timestamp})
	message(STATUS "Git has changes: " ${git_has_changes})
	message(STATUS "Git status:\n" ${git_status})
else()
	message(STATUS "Git not found")
endif()

set(desc "${git_describe} @${git_branch} (${git_rev}) ${git_has_changes}")
if("${desc}" STREQUAL "${prevdesc}")
	message(STATUS "Git metadata unchanged, skipping configure")
else()
	message(STATUS "configuring ${src} -> ${dst}")
	configure_file("${src}" "${dst}" @ONLY)
	file(WRITE "${descfile}" "${desc}")
endif()
