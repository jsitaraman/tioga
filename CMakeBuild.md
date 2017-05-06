
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
