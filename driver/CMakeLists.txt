
set(TIOGA_EXE_SOURCES
  modules.f90
  gridutils.f90)

add_executable(tioga.exe testTioga.f90 ${TIOGA_EXE_SOURCES})
add_executable(tioga_read.exe testTioga_dataRead.f90 ${TIOGA_EXE_SOURCES})

target_compile_definitions(tioga.exe PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(tioga.exe PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(tioga.exe
  tioga ${MPI_Fortran_LIBRARIES} ${CMAKE_DL_LIBS})
set_target_properties(tioga.exe PROPERTIES LINKER_LANGUAGE Fortran)

target_compile_definitions(tioga_read.exe PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(tioga_read.exe PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(tioga_read.exe
  tioga ${MPI_Fortran_LIBRARIES} ${CMAKE_DL_LIBS})
set_target_properties(tioga_read.exe PROPERTIES LINKER_LANGUAGE Fortran)

install(TARGETS tioga.exe
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(TARGETS tioga_read.exe
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
