add_library(tioga)

target_sources(tioga PRIVATE
  # Fortran sources
  kaiser.f
  cellVolume.f90
  median.F90

  # C sources
  buildADTrecursion.c
  linklist.c
  math.c
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

target_include_directories(tioga PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

if (TIOGA_HAS_NODEGID)
  target_compile_definitions(tioga PUBLIC TIOGA_HAS_NODEGID)
endif()

if (TIOGA_ENABLE_TIMERS)
  target_compile_definitions(tioga PUBLIC TIOGA_ENABLE_TIMERS)
endif()

if (TIOGA_OUTPUT_STATS)
  target_compile_definitions(tioga PUBLIC TIOGA_OUTPUT_STATS)
endif()

target_link_libraries(tioga INTERFACE ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
if (TIOGA_ENABLE_ARBORX)
  target_link_libraries(tioga PRIVATE ArborX::ArborX)
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(tioga PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tioga PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if (TIOGA_ENABLE_CUDA)
  separate_arguments(TIOGA_CUDA_FLAGS)
  target_compile_definitions(tioga PUBLIC
    TIOGA_HAS_GPU TIOGA_HAS_CUDA)
  target_compile_options(tioga PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:${TIOGA_CUDA_FLAGS}>)
  set(CMAKE_CUDA_FLAGS ${TIOGA_CUDA_ARCH_FLAGS})
  set_target_properties(tioga PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON)
  get_target_property(TIOGALIB_CPP_SRC tioga SOURCES)
  list(FILTER TIOGALIB_CPP_SRC INCLUDE REGEX "\\.C")
  set_source_files_properties(${TIOGALIB_CPP_SRC} PROPERTIES LANGUAGE CUDA)
endif()

if (TIOGA_ENABLE_HIP)
  separate_arguments(TIOGA_HIP_ARCH_FLAGS)
  target_compile_definitions(tioga PUBLIC
    TIOGA_HAS_GPU TIOGA_HAS_HIP)
  target_compile_options(tioga PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${TIOGA_HIP_ARCH_FLAGS}>)
  target_compile_options(tioga PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-x hip>)
endif()

# Enable position independent code
set_target_properties(tioga PROPERTIES POSITION_INDEPENDENT_CODE ON)


file(GLOB TIOGA_HEADERS *.h)
install(TARGETS tioga
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES ${TIOGA_HEADERS}
  DESTINATION include)
