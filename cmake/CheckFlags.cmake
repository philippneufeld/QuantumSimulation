# Philipp Neufeld, 2021-2022

#	This function checks if a program (SOURCE) can be compiled given the flags (FLAGS)
function(try_compile_with_flags RESULT SOURCE FLAGS)
		set(TEST_FILE_CPP "${CMAKE_BINARY_DIR}/CMakeFiles/try_compile_with_flags.cpp")
		set(TEST_FILE_EXE "${CMAKE_BINARY_DIR}/CMakeFiles/try_compile_with_flags${CMAKE_EXECUTABLE_SUFFIX}")
		file(WRITE "${TEST_FILE_CPP}" "${SOURCE}")

		set(CMAKE_TRY_COMPILE_TARGET_TYPE, EXECUTABLE)
		try_compile("${RESULT}"
				"${CMAKE_BINARY_DIR}"
				"${TEST_FILE_CPP}"
				OUTPUT_VARIABLE TRY_COMPILE_OUTPUT
				COMPILE_DEFINITIONS "${FLAGS}"
				COPY_FILE "${TEST_FILE_EXE}"
		)

		#	message("TRY_COMPILE_OUTPUT: ${TRY_COMPILE_OUTPUT}")

		file(REMOVE "${TEST_FILE_CPP}")
		file(REMOVE "${TEST_FILE_EXE}")
endfunction()

#	Returns first flag that is passed that makes SOURCE compileable
function(select_flag RESULT SOURCE OLD_FLAGS)

		set(FLAG_FOUND FALSE)
		set(${RESULT} "" PARENT_SCOPE)

		set(FLAG_OPTIONS "")
		foreach(testflag IN LISTS ARGN)
				if("${FLAG_OPTIONS}" STREQUAL "")
						set(FLAG_OPTIONS "\"${testflag}\"")
				else()
						set(FLAG_OPTIONS "${FLAG_OPTIONS}, \"${testflag}\"")
				endif()

				if (NOT FLAG_FOUND)
						unset(TRY_COMPILE_FOUND)
						try_compile_with_flags(TRY_COMPILE_FOUND "${SOURCE}" "${OLD_FLAGS} ${testflag}")
						
						if (TRY_COMPILE_FOUND)
								set("${RESULT}" "${testflag}" PARENT_SCOPE)
								set(FLAG_FOUND TRUE)
						endif()
				endif()
		endforeach()

		if(NOT FLAG_FOUND)
				message(FATAL_ERROR "No working flag found in ${FLAG_OPTIONS}")
		endif()
endfunction()
