# Philipp Neufeld, 2021

cmake_minimum_required (VERSION 3.17)

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

# source files
set(QSIM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/QSim")
set(QSIM_SOURCES
    "${QSIM_SOURCE_DIR}/Python/PythonInterpreter.cpp"
    "${QSIM_SOURCE_DIR}/Python/Plotting.cpp"
    "${QSIM_SOURCE_DIR}/Util/Argparse.cpp"
    "${QSIM_SOURCE_DIR}/Util/CalcApp.cpp"
    "${QSIM_SOURCE_DIR}/Util/Memory.cpp"
    "${QSIM_SOURCE_DIR}/Util/DataFile.cpp"
    "${QSIM_SOURCE_DIR}/Util/CLIProgressBar.cpp")

# find required libs
find_package(Threads REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# compiler options
set(QSIM_PUBLIC_COMPILER_DEFINITIONS "")
set(QSIM_PUBLIC_INCLUDE_DIRS ${QSIM_INCLUDE_DIR} ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
set(QSIM_PUBLIC_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${Python3_LIBRARIES})

# GNUC on an ARM system (disable ABI warning)
if (CMAKE_COMPILER_IS_GNUCC AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    message("ARM system and GNUC compiler detected. Disabling ABI warnings (-Wno-psabi)")
    set(QSIM_PUBLIC_COMPILER_OPTIONS "-Wno-psabi")
else()
    set(QSIM_PUBLIC_COMPILER_OPTIONS "")
endif()

# executable suffix
set(QSIM_EXE_SUFFIX ".out")

function(qsim_add_app APP_TARGET)
    add_subdirectory("apps/${APP_TARGET}")
    target_include_directories("${APP_TARGET}" PUBLIC ${QSIM_PUBLIC_INCLUDE_DIRS})
    target_compile_options("${APP_TARGET}" PUBLIC ${QSIM_PUBLIC_COMPILER_OPTIONS})
    target_compile_definitions("${APP_TARGET}" PUBLIC ${QSIM_PUBLIC_COMPILER_DEFINITIONS})
    target_link_libraries("${APP_TARGET}" ${QSIM_PUBLIC_LINK_LIBRARIES})
    set_target_properties("${APP_TARGET}" PROPERTIES SUFFIX ${QSIM_EXE_SUFFIX})
endfunction()

# add apps
qsim_add_app("Test")
qsim_add_app("Quadrature")
qsim_add_app("Rb87_2lvl")
qsim_add_app("Rb87_EIT")
qsim_add_app("Rb87_SAS")
qsim_add_app("Rb87_MW")
qsim_add_app("NO_RydEx")
