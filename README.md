# QuantumSimulation

This project includes methods of calculating the dynamics as well as the steady-state of simple N-level systems.

Various apps are included that showcase the capabilities of the library.

## Dependencies

This project depends on a bunch of third-party libraries, that should be installed on the system:

- [Eigen3](http://eigen.tuxfamily.org)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [Python 3](https://www.python.org/) (optional)


## Compilation

First navigate into the project directory and initialize the build directory by running

    cmake -S . -B build -DCMAKE_BUILD_TYPE="Release" -DCMAKE_C_COMPILER="$(which clang)" -DCMAKE_CXX_COMPILER="$(which clang++)"

Here, clang/clang++ can be replaced by another compiler (e.g. gcc/g++).
Next, navigate into the build directory and compile all apps by using

    cd build && make

Alternatively the apps can also be compiled separately by running the above `make` command in the build folder of a specific app.
