
set(GRIDGEN_SOURCES
  modules.f
  buildGrid.f
  clust.f
  partSphere.f90
  findfact.f90)

add_executable(buildGrid ${GRIDGEN_SOURCES})

target_compile_options(buildGrid PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(buildGrid PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(buildGrid
  tioga ${MPI_Fortran_LIBRARIES} ${CMAKE_DL_LIBS})
set_target_properties(buildGrid PROPERTIES LINKER_LANGUAGE Fortran)

install(TARGETS buildGrid
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
