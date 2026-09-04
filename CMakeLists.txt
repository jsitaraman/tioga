cmake_minimum_required(VERSION 2.8)

project(TIOGA CXX C Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" on)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" on)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" on)
option(TIOGA_HAS_NODEGID "Support node global IDs (default: on)" on)
option(TIOGA_ENABLE_UNIQUEID
  "De-duplicate repeated donor search query points on the host (default: on)" on)
option(TIOGA_ENABLE_TIMERS "Track timing information for TIOGA (default: off)" OFF)
option(TIOGA_OUTPUT_STATS "Output statistics for TIOGA holecutting (default: off)" OFF)
option(TIOGA_ENABLE_CUDA "Build the cuBQL based GPU donor search (default: off)" OFF)
set(TIOGA_SEARCH_BACKEND "cpu" CACHE STRING
  "Donor search backend: cpu, adt_gpu, cubql or cubql_batch")
set_property(CACHE TIOGA_SEARCH_BACKEND PROPERTY STRINGS
  cpu adt_gpu cubql cubql_batch)
set(TIOGA_CUBQL_DIR "/scratch/software/cuBQL" CACHE PATH
  "Root of the cuBQL source tree (the directory containing cuBQL/bvh.h)")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Ensure C++11 standard is enabled
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "-g -rdynamic -std=c++0x ${CMAKE_CXX_FLAGS}")
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
endif()

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fbounds-check -fbacktrace -fdefault-real-8 -ffree-line-length-none")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128")
endif()

if (APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

if (BUILD_SHARED_LIBS)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()

if (TIOGA_HAS_NODEGID)
  add_definitions(-DTIOGA_HAS_NODEGID)
endif()

if (TIOGA_ENABLE_UNIQUEID)
  add_definitions(-DTIOGA_ENABLE_UNIQUEID)
endif()

if (TIOGA_ENABLE_TIMERS)
  add_definitions(-DTIOGA_ENABLE_TIMERS)
endif()

if (TIOGA_OUTPUT_STATS)
  add_definitions(-DTIOGA_OUTPUT_STATS)
endif()

# CUDA / cuBQL for the GPU donor search
if (TIOGA_ENABLE_CUDA)
  if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    if (${CMAKE_VERSION} VERSION_LESS 3.24)
      set(CMAKE_CUDA_ARCHITECTURES "80")
    else()
      set(CMAKE_CUDA_ARCHITECTURES "native")
    endif()
  endif()
  enable_language(CUDA)
  if (NOT EXISTS "${TIOGA_CUBQL_DIR}/cuBQL/bvh.h")
    message(FATAL_ERROR
      "TIOGA_ENABLE_CUDA is ON but cuBQL was not found at TIOGA_CUBQL_DIR="
      "${TIOGA_CUBQL_DIR} (expected ${TIOGA_CUBQL_DIR}/cuBQL/bvh.h). "
      "Point -DTIOGA_CUBQL_DIR=<path> at a cuBQL checkout.")
  endif()
  add_definitions(-DTIOGA_ENABLE_CUDA)
  message(STATUS "#tioga: CUDA search enabled, arch=${CMAKE_CUDA_ARCHITECTURES}, "
                 "cuBQL=${TIOGA_CUBQL_DIR}")
endif()

# Which donor search MeshBlock::search() runs. Compile time rather than run
# time so that a build contains exactly one backend and the comparison between
# them cannot be skewed by dispatch left in the measured path.
string(TOLOWER "${TIOGA_SEARCH_BACKEND}" TIOGA_SEARCH_BACKEND_LC)
if (TIOGA_SEARCH_BACKEND_LC STREQUAL "cpu")
  add_definitions(-DTIOGA_SEARCH_BACKEND_CPU)
elseif (TIOGA_SEARCH_BACKEND_LC STREQUAL "adt_gpu")
  add_definitions(-DTIOGA_SEARCH_BACKEND_ADT_GPU)
elseif (TIOGA_SEARCH_BACKEND_LC STREQUAL "cubql")
  add_definitions(-DTIOGA_SEARCH_BACKEND_CUBQL)
elseif (TIOGA_SEARCH_BACKEND_LC STREQUAL "cubql_batch")
  add_definitions(-DTIOGA_SEARCH_BACKEND_CUBQL_BATCH)
else()
  message(FATAL_ERROR
    "TIOGA_SEARCH_BACKEND must be one of cpu, adt_gpu, cubql, cubql_batch "
    "(got '${TIOGA_SEARCH_BACKEND}')")
endif()

if (NOT TIOGA_SEARCH_BACKEND_LC STREQUAL "cpu")
  if (NOT TIOGA_ENABLE_CUDA)
    message(FATAL_ERROR
      "TIOGA_SEARCH_BACKEND=${TIOGA_SEARCH_BACKEND} needs -DTIOGA_ENABLE_CUDA=ON")
  endif()
  # The GPU backends search every query point, so they only pay off once the
  # host de-duplication pass is compiled out. With it on they stand down and
  # search() runs its own loop.
  if (TIOGA_ENABLE_UNIQUEID)
    message(WARNING
      "#tioga: TIOGA_SEARCH_BACKEND=${TIOGA_SEARCH_BACKEND} is inactive while "
      "TIOGA_ENABLE_UNIQUEID is on; the host search will be used. "
      "Configure with -DTIOGA_ENABLE_UNIQUEID=off to enable it.")
  endif()
endif()
message(STATUS "#tioga: search backend=${TIOGA_SEARCH_BACKEND_LC}, "
               "uniqueid=${TIOGA_ENABLE_UNIQUEID}")

# Always build libtioga
add_subdirectory(src)

# Optionally build driver exe and gridGen if the user requests it
if (BUILD_TIOGA_EXE)
  add_subdirectory(driver)
endif()

if (BUILD_GRIDGEN_EXE)
  add_subdirectory(gridGen)
endif()

# CMake installation configuration

install(EXPORT TIOGALibraries
  DESTINATION lib/cmake/TIOGA
  FILE TIOGALibraries.cmake)

# Create TIOGA config so that other codes can find TIOGA
include(CMakePackageConfigHelpers)
set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)
configure_package_config_file(
  cmake/TIOGAConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
  INSTALL_DESTINATION lib/cmake/TIOGA
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
  DESTINATION lib/cmake/TIOGA)
