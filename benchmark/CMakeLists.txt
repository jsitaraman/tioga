# Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
add_executable(bench_search.exe bench_search.C)
target_include_directories(bench_search.exe PRIVATE
  ${CMAKE_SOURCE_DIR}/src ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(bench_search.exe tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

install(TARGETS bench_search.exe
  RUNTIME DESTINATION bin)

add_executable(bench_multiblock.exe bench_multiblock.C)
target_include_directories(bench_multiblock.exe PRIVATE
  ${CMAKE_SOURCE_DIR}/src ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(bench_multiblock.exe tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

install(TARGETS bench_multiblock.exe
  RUNTIME DESTINATION bin)
