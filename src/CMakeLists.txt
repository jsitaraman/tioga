

set(TIOGA_SOURCES
  # Fortran sources
  kaiser.f
  cellVolume.f90
  median.F90

  # C sources
  buildADTrecursion.c
  linklist.c
  math.c
  get_amr_index_xyz.c
  utils.c

  # CXX sources
  ADT.C
  CartBlock.C
  CartGrid.C
  MeshBlock.C
  bookKeeping.C
  cartOps.C
  checkContainment.C
  dataUpdate.C
  dataUpdateGPU.C
  exchangeAMRDonors.C
  exchangeBoxes.C
  exchangeDonors.C
  exchangeSearchData.C
  getCartReceptors.C
  highOrder.C
  holeMap.C
  parallelComm.C
  search.C
  searchADTrecursion.C
  tioga.C
  tiogaInterface.C
  )

file(GLOB TIOGA_HEADERS *.h)

if(BUILD_CUDA_SUPPORT)
  # Load the CUDA cmake variables
  find_package(CUDA REQUIRED)
  # Set architecture flags
  if((CUDA_VERSION_MAJOR EQUAL 8) OR (CUDA_VERSION_MAJOR GREATER 8))
    set(CUDA_OPTS "-arch=sm_30")
    list(APPEND CUDA_OPTS "-gencode=arch=compute_30,code=sm_30")
    list(APPEND CUDA_OPTS "-gencode=arch=compute_60,code=sm_60")
  elseif((CUDA_VERSION_MAJOR EQUAL 7) OR (CUDA_VERSION_MAJOR GREATER 7))
    set(CUDA_OPTS "-arch=sm_30")
    list(APPEND CUDA_OPTS "-gencode=arch=compute_30,code=sm_30")
  else()
    set(CUDA_OPTS "-arch=sm_20")
    list(APPEND CUDA_OPTS "-gencode=arch=compute_20,code=sm_20")
  endif()
  message(STATUS "CUDA target : ${CUDA_OPTS}")
  # Basic CUDA settings
  set(CUDA_SEPARABLE_COMPILATION ON CACHE INTERNAL "Compile GPU code seperately")
  set(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE INTERNAL "Do not propagate host flags")
  list(APPEND CUDA_OPTS "-O2")
  list(APPEND TIOGA_DEFINITIONS "-DUSE_CUDA")
  list(APPEND CUDA_OPTS -ccbin ${CMAKE_C_COMPILER} --compiler-options -fPIC)

  # Get all the CUDA files
  file(GLOB cu_sources "*.cu")
  message(STATUS "sources: ${cu_sources}")
  cuda_include_directories( ${MPI_C_INCLUDE_PATH} )
  # cuda_add_library(gpufunc STATIC ${cu_sources} OPTIONS ${CUDA_OPTS} ${TIOGA_DEFINITIONS})
  # set_target_properties(gpufunc PROPERTIES
  #   ARCHIVE_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/lib
  #   LIBRARY_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/lib
  #   )

  cuda_add_library(tioga STATIC ${TIOGA_SOURCES} ${cu_sources} OPTIONS ${CUDA_OPTS} ${TIOGA_DEFINITIONS})
  set_target_properties(tioga PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    LIBRARY_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    )
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  message(STATUS "cudalibs: ${CUDA_LIBRARIES}")
  #
  # Static libraries dont actually link to external libraries... BUT
  # if you tell cmake it does, cmake will remember and add those
  # libraries when the ".a" file is linked somewhere else. If the ".a"
  # file is linked externally (like form another project), just make
  # sure the shared object or executable also links these libraries:
  #
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} ${CUDA_LIBRARIES})

else()

  # add_library(tioga SHARED ${TIOGA_SOURCES})
  add_library(tioga ${TIOGA_SOURCES})
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

endif()



# set_target_properties(tioga PROPERTIES COMPILE_FLAGS "-std=c++11 -fPIC")

if(MPI_COMPILE_FLAGS)
  set_target_properties(tioga PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tioga PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# if(BUILD_CUDA_SUPPORT)
#   install(TARGETS tioga gpufunc
#     EXPORT "${CMAKE_PROJECT_NAME}Libraries"
#     RUNTIME DESTINATION bin
#     ARCHIVE DESTINATION lib
#     LIBRARY DESTINATION lib)
# else()
  install(TARGETS tioga
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
# endif()

install(FILES ${TIOGA_HEADERS}
  DESTINATION include)
