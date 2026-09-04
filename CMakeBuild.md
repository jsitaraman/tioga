
# Building and installing TIOGA using CMake 

TIOGA has been configured to use CMake to configure, build, and install the
library. The user can choose standard CMake options as well as additional
TIOGA-specific options to customize the build and installation process. A brief
description of the CMake-based build process and configuration options are
described in this document.

The minimal dependencies to build TIOGA on your system are CMake, a working C,
C++, and Fortran compilers as well as an MPI library along with its headers. If
the dependencies are satisfied, then execute the following commands to clone and
build TIOGA:

```
git clone <TIOGA_GITHUB_URL>
cd tioga 

# Create a build directory
mkdir build

# Configure build using auto-discovered parameters
cmake ../

# Build the library
make
```

When the steps are successfully executed, the compiled static library is located
in `tioga/build/src/libtioga.a`. 

## Building `driver` and `gridGen` executables

By default, CMake does not build the `tioga.exe` driver code or the `buildGrid`
executable. To enable these at configure phase:

```
cmake -DBUILD_TIOGA_EXE:BOOL=ON -DBUILD_GRIDGEN_EXE:BOOL=ON ../
```

followed by `make`. The executables will be located in `build/driver/tioga.exe`
and `build/gridGen/buildGrid` respectively.

## Customizing compilers 

To use different compilers other than what is detected by CMake use the
following configure command:

```
CC=mpicc CXX=mpicxx FC=mpif90 cmake ../ 
```

## Release, Debug, and other compilation options

Use `-DCMAKE_BUILD_TYPE` with `Release`, `Debug` or `RelWithDebInfo` to build
with different optimization or debugging flags. For example,

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

You can also use `CMAKE_CXX_FLAGS`, `CMAKE_C_FLAGS`, and `CMAKE_Fortran_FLAGS`
to specify additional compile time flags of your choosing. For example,

```
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_Fortran_FLAGS="-fbounds-check -fbacktrace" \
  ../
```

## Custom install location

Finally, it is usually desirable to specify the install location when using
`make install` when using TIOGA with other codes.

```
# Configure TIOGA several options
CC=mpicc CXX=mpicxx FC=mpif90 cmake \
  -DCMAKE_INSTALL_PREFIX=${HOME}/software/ \
  -DBUILD_TIOGA_EXE=ON \
  -DBUILD_GRIDGEN_EXE=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_Fortran_FLAGS="-fbounds-check -fbacktrace" \
  ../
  
# Compile library and install at user-defined location
make && make install
```

## GPU donor search

The donor search in `MeshBlock::search()` has four interchangeable
implementations, selected at configure time with `TIOGA_SEARCH_BACKEND`:

| Backend       | What runs on the device                                     |
|---------------|-------------------------------------------------------------|
| `cpu`         | nothing; the host search (default)                          |
| `adt_gpu`     | the ADT query loop, over the host-built tree                |
| `cubql`       | a cuBQL BVH over all cells of a block, built and queried    |
| `cubql_batch` | as `cubql`, but one BVH per rank across all of its blocks   |

All three GPU backends need `-DTIOGA_ENABLE_CUDA=ON` and a
[cuBQL](https://github.com/NVIDIA/cuBQL) checkout, which is header only here:

```
cmake \
  -DTIOGA_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90 \
  -DTIOGA_CUBQL_DIR=/path/to/cuBQL \
  -DTIOGA_SEARCH_BACKEND=cubql \
  -DTIOGA_ENABLE_UNIQUEID=off \
  ../
```

A backend may decline at run time, for instance for high order elements whose
containment test goes through host callbacks. The host search runs instead, so
the result is always correct whatever the configuration.

### De-duplication of query points

`TIOGA_ENABLE_UNIQUEID` (on by default) controls whether repeated query points
are detected on the host so that only the first occurrence is searched. The GPU
backends search every point regardless, so they only pay off once that pass is
compiled out, and they stand down entirely while it is on. Configuring a GPU
backend with `TIOGA_ENABLE_UNIQUEID=on` warns and leaves the host search in
place.

### Comparing the backends

Because the backend is a configure time choice, backends cannot be compared
within a single run. `scripts/compare_search_backends.sh` builds one tree per
backend, runs `case/` through each, checks that they all find the same donors
as the host search, and prints the cost of each search phase:

```
scripts/compare_search_backends.sh -n 8 -c /path/to/cuBQL
```

Use `-u` to repeat the comparison with de-duplication left on, which checks
that the host path is unaffected and that the GPU backends stand down.

The same instrumentation is available directly, in any build:

| Variable               | Effect                                            |
|------------------------|---------------------------------------------------|
| `TIOGA_DONOR_DUMP=dir` | write the donor found for every query point       |
| `TIOGA_SEARCH_TIMERS=1`| print the per-phase search cost, reduced over ranks|
| `TIOGA_SEARCH_REPEAT=n` | repeat the search `n` times and time the last pass |

`TIOGA_SEARCH_REPEAT` exists because the drivers call `performConnectivity()`
once, which measures the cuBQL backends on a cold tree. Passes after the first
mark the coordinates dirty, exactly as a moving mesh does, so the BVH is refit
rather than rebuilt and the steady state cost is what gets reported.
