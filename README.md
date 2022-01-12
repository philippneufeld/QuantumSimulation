# QuantumSimulation

This project includes methods of calculating the dynamics as well as the steady-state of simple N-level systems.

Various apps are included that showcase the capabilities of the library.

## Compilation

First navigate into the project directory and initialize the build directory by running

    cmake -S . -B build -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_COMPILER="$(which clang)" -DCMAKE_CXX_COMPILER="$(which clang++)"

Here, clang/clang++ can be replaced by another compiler (e.g. gcc/g++). <br/>
Next, navigate into the build directory and compile all apps by using

    make

Alternatively the apps can also be compiled separately by running the above `make` command in the build folder of a specific app.