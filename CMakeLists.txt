# Philipp Neufeld, 2021

cmake_minimum_required (VERSION 3.18)

# don't allow build in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory! Please use cmake [<options>] -S <path-to-source> -B <path-to-build>")
endif()

# define project
project ("QuantumSimulation" VERSION 1.0)

# use C++14
set (CMAKE_CXX_STANDARD 14)

# include directories
set(QSIM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}")
set(QSIM_PUBLIC_INCLUDE_DIRS, "${QSIM_INCLUDE_DIR}")

# find required libs
find_package (Threads REQUIRED)

# compiler options
set(QSIM_COMPILER_OPTIONS "")
set(QSIM_COMPILER_DEFINITIONS "")
set(QSIM_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")

# executable suffix
set(QSIM_EXE_SUFFIX ".out")

function(qsim_add_app APP_TARGET)
    add_subdirectory("apps/${APP_TARGET}")
    target_include_directories("${APP_TARGET}" PUBLIC "${QSIM_PUBLIC_INCLUDE_DIRS}")
    target_compile_options("${APP_TARGET}" PUBLIC "${QSIM_COMPILER_OPTIONS}")
    target_compile_definitions("${APP_TARGET}" PUBLIC "${QSIM_COMPILER_DEFINITIONS}")
    target_link_libraries("${APP_TARGET}" "${QSIM_LINK_LIBRARIES}")
    set_target_properties("${APP_TARGET}" PROPERTIES SUFFIX "${QSIM_EXE_SUFFIX}")
endfunction()

# add apps
qsim_add_app("Test")
