cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(TIOGA CXX C Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" on)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" on)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" on)
option(TIOGA_HAS_NODEGID "Support node global IDs (default: on)" on)
option(TIOGA_ENABLE_UNIQUEID "Support evaluation of global unique IDs for nearbody mesh (default: on)" on)
option(TIOGA_ENABLE_TIMERS "Track timing information for TIOGA (default: off)" ON)
option(TIOGA_OUTPUT_STATS "Output statistics for TIOGA holecutting (default: off)" OFF)
option(TIOGA_ENABLE_CUDA "Enable CUDA support (default: off)" OFF)
option(TIOGA_ENABLE_HIP "Enable AMD HIP support (default: off)" OFF)
option(TIOGA_FAKE_GPU "Fake GPU pathways for testing (default: off)" OFF)

# CUDA specific options
set(TIOGA_CUDA_SM "70" CACHE STRING "CUDA arch option")
set(TIOGA_HIP_ARCH "gfx906,gfx908" CACHE STRING "HIP arch option")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if (TIOGA_ENABLE_CUDA)
  enable_language(CUDA)
  if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}. "
      "Please use CUDA toolkit version 9.0 or newer")
  endif()
  set(TIOGA_CUDA_FLAGS "-g -rdynamic --Wno-deprecated-gpu-targets -m64")
  set(TIOGA_CUDA_ARCH_FLAGS "-gencode=arch=compute_${TIOGA_CUDA_SM},code=sm_${TIOGA_CUDA_SM}")
  set(CMAKE_CUDA_STANDARD 11)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_EXTENSIONS OFF)
endif()

if (TIOGA_ENABLE_HIP)
  if (NOT DEFINED HIP_PATH)
    set(HIP_PATH $ENV{HIP_PATH})
  endif()
  set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" "${CMAKE_MODULE_PATH}")
  find_package(HIP REQUIRED)
  if ("${HIP_COMPILER}" STREQUAL "hcc")
    message(FATAL_ERROR "TIOGA does not support deprecated HCC compiler. Use newer ROCm")
  endif()
  message(STATUS "HIP: ${HIP_VERSION}, platform = ${HIP_PLATFORM}, compiler = ${HIP_COMPILER}")
  set(TIOGA_HIP_ARCH_FLAGS "-m64 --amdgpu-target=${TIOGA_HIP_ARCH}")
endif()

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fbounds-check -fbacktrace -fdefault-real-8 -ffree-line-length-none")
 set(CMAKE_CXX_FLAGS
   "${CMAKE_CXX_FLAGS} -g -rdynamic")
 set(CMAKE_C_FLAGS
   "${CMAKE_C_FLAGS} -g -rdynamic")
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
