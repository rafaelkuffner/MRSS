
cmake_minimum_required(VERSION 3.1)

message(STATUS "Retrieving git revision info")

set(git_rev unknown)
set(git_describe unknown)
set(git_timestamp unknown)

find_package(Git)

set(revfile "${dst}.rev")
set(prev "")
if(EXISTS "${revfile}")
	file(READ "${revfile}" prev)
endif()

if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
		OUTPUT_VARIABLE git_rev
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
	message(STATUS "Git rev: " ${git_rev})
	message(STATUS "Git describe: " ${git_describe})
	message(STATUS "Git timestamp: " ${git_timestamp})
else()
	message(STATUS "Git not found")
endif()

if("${git_rev}" STREQUAL "${prev}")
	message(STATUS "Git rev up-to-date, skipping configure")
else()
	message(STATUS "configuring ${src} -> ${dst}")
	configure_file("${src}" "${dst}" @ONLY)
	file(WRITE "${revfile}" ${git_rev})
endif()
