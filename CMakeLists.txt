# Philipp Neufeld, 2021-2022

cmake_minimum_required (VERSION 3.17)

# don't allow build in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory! Please use cmake [<options>] -S <path-to-source> -B <path-to-build>")
endif()

# define project
project ("QuantumSimulation" VERSION 1.0)

# include directories
set(QSIM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}")

# add static library
set(QSIM_LIB_TARGET "${CMAKE_PROJECT_NAME}")
set(QSIM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/QSim")
set(QSIM_SOURCES
    
    "${QSIM_SOURCE_DIR}/Execution/ThreadPool.cpp"
    "${QSIM_SOURCE_DIR}/Python/PythonInterpreter.cpp"
    "${QSIM_SOURCE_DIR}/Python/Plotting.cpp"
    "${QSIM_SOURCE_DIR}/Util/Argparse.cpp"
    "${QSIM_SOURCE_DIR}/Util/SimulationApp.cpp"
    "${QSIM_SOURCE_DIR}/Util/Memory.cpp"
    "${QSIM_SOURCE_DIR}/Util/DataFile.cpp"
    "${QSIM_SOURCE_DIR}/Util/Progress.cpp"
    "${QSIM_SOURCE_DIR}/Util/ProgressBar.cpp"
    )
add_library("${QSIM_LIB_TARGET}" STATIC "${QSIM_SOURCES}")

# require C++17 for library and all children using it
target_compile_features("${QSIM_LIB_TARGET}" PUBLIC cxx_std_17)

# set include directory as interface
target_include_directories("${QSIM_LIB_TARGET}" INTERFACE "${QSIM_INCLUDE_DIR}")

# threading
find_package(Threads REQUIRED)
if (Threads_FOUND)
    target_link_libraries("${QSIM_LIB_TARGET}" PUBLIC "${CMAKE_THREAD_LIBS_INIT}")
endif()

# hdf5 (standardized file format for storing data)
find_package(HDF5 REQUIRED)
if (HDF5_FOUND)
    target_include_directories("${QSIM_LIB_TARGET}" PUBLIC "${HDF5_INCLUDE_DIRS}")
    target_link_libraries("${QSIM_LIB_TARGET}" PUBLIC "${HDF5_LIBRARIES}")
    # target_compile_definitions("${QSIM_LIB_TARGET}" PUBLIC "${HDF5_DEFINITIONS}")
endif()

# eigen3 (linear algebra library)
find_package(Eigen3 REQUIRED NO_MODULE)
if (EIGEN3_FOUND)
    # fixes weird bug in FindEigen3 on the PI5 simulation machines (OpenSUSE)
    string(REGEX REPLACE "^/usr//usr" "/usr" EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR}")
    message("Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
    target_include_directories("${QSIM_LIB_TARGET}" PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()

# python (optional)
find_package(Python3 COMPONENTS Interpreter Development NumPy)
if (Python3_FOUND)
    target_include_directories("${QSIM_LIB_TARGET}" PUBLIC "${Python3_INCLUDE_DIRS}" "${Python3_NumPy_INCLUDE_DIRS}")
    target_link_libraries("${QSIM_LIB_TARGET}" PUBLIC "${Python3_LIBRARIES}")
    target_compile_definitions("${QSIM_LIB_TARGET}" PUBLIC "-DQSIM_PYTHON3")
endif()

# GNUC on an ARM system (disable ABI warning)
if (CMAKE_COMPILER_IS_GNUCC AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    message("ARM system and GNUC compiler detected. Disabling ABI warnings (-Wno-psabi)")
    target_compile_options("${QSIM_LIB_TARGET}" PUBLIC "-Wno-psabi")
endif()

# put the built library file into the .../lib/ folder
set(BUILD_LIB_DIR "${CMAKE_BINARY_DIR}/lib")
set_target_properties("${QSIM_LIB_TARGET}" 
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${BUILD_LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${BUILD_LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${BUILD_LIB_DIR}"
)

set(QSIM_LIB_NAME "${CMAKE_STATIC_LIBRARY_PREFIX}${QSIM_LIB_TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(QSIM_LIB_PATH "${BUILD_LIB_DIR}/${QSIM_LIB_NAME}")

# function that adds a target for a qsim app
function(qsim_add_app APP_TARGET)
    add_subdirectory("apps/${APP_TARGET}")
    target_link_libraries("${APP_TARGET}" PUBLIC "${QSIM_LIB_TARGET}")
endfunction()

# add apps
qsim_add_app("Test")
qsim_add_app("Quadrature")
qsim_add_app("CurveFit")
qsim_add_app("Rb87_2lvl")
qsim_add_app("Rb87_EIT")
qsim_add_app("Rb87_SAS")
qsim_add_app("NO_RydEx")
