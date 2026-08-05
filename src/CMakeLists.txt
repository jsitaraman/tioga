

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
  exchangeAMRDonors.C
  exchangeBoxes.C
  exchangeDonors.C
  exchangeSearchData.C
  getCartReceptors.C
  highOrder.C
  holeMap.C
  linCartInterp.C
  parallelComm.C
  search.C
  searchADTrecursion.C
  tioga.C
  tiogaInterface.C
  )

if (TIOGA_ENABLE_CUDA)
  list(APPEND TIOGA_SOURCES searchGPU.cu)
else()
  list(APPEND TIOGA_SOURCES searchGPU_stub.C)
endif()

file(GLOB TIOGA_HEADERS *.h)

add_library(tioga ${TIOGA_SOURCES})
target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

if (TIOGA_ENABLE_CUDA)
  # cuBQL is header-only here: searchGPU.cu defines
  # CUBQL_GPU_BUILDER_IMPLEMENTATION and instantiates the float3 builders
  # it needs, so no cuBQL library has to be built or linked.
  target_include_directories(tioga PUBLIC ${TIOGA_CUBQL_DIR})
  set_source_files_properties(searchGPU.cu PROPERTIES LANGUAGE CUDA)
  set_target_properties(tioga PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_SEPARABLE_COMPILATION OFF
    POSITION_INDEPENDENT_CODE ON)
  target_compile_options(tioga PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3>)
  # nvcc needs the MPI headers for codetypes.h
  target_include_directories(tioga PRIVATE ${MPI_CXX_INCLUDE_DIRS})
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(tioga PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tioga PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS tioga
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES ${TIOGA_HEADERS}
  DESTINATION include)
