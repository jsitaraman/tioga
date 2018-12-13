cmake_minimum_required(VERSION 2.8)

project(TIOGA CXX C Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" on)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" on)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" on)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Ensure C++11 standard is enabled
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "-g -rdynamic -std=c++0x" ${CMAKE_CXX_FLAGS})
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_FLAGS "-O3" ${CMAKE_CXX_FLAGS})
endif()

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fbounds-check -fbacktrace -fdefault-real-8 -ffree-line-length-none")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128")
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