# Nitric oxide simulation

Code to run density matrix simulations and calculation of Stark maps for nitric oxide, developed during my Master's thesis time.
This project includes methods of calculating the dynamics as well as the steady-state of simple N-level systems.
Various apps are included that showcase the capabilities of the library.

## Dependencies

This project depends on a bunch of third-party libraries, that should be installed on the system:

- [Eigen3](http://eigen.tuxfamily.org)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)


## Compilation

If you have just cloned the repository and cmake did not run yet execute `init.sh`.
Otherwise, navigate into the build directory and compile all apps by using

    cd build && make

Alternatively the apps can also be compiled separately by running the above `make` command in the build folder of a specific app.
