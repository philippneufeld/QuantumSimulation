#!/usr/bin/env bash

# initialize cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_COMPILER="$(which clang)" -DCMAKE_CXX_COMPILER="$(which clang++)"

# compile for the first time
cd build && make && cd ..
