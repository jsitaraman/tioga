add_executable(bench_search.exe bench_search.C)
target_include_directories(bench_search.exe PRIVATE
  ${CMAKE_SOURCE_DIR}/src ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(bench_search.exe tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

install(TARGETS bench_search.exe
  RUNTIME DESTINATION bin)
