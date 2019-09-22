

find_package(PythonInterp)
find_package(PythonLibs)
if(NOT ${PYTHONLIBS_FOUND})
  find_package(PythonInterp 2.7 REQUIRED)
  find_package(PythonLibs 2.7 REQUIRED)
endif()
if(${PYTHON_VERSION_STRING} VERSION_LESS "3")
  set(PYTHON_2_COMPATIBILITY "from __future__ import absolute_import, division, print_function")
else()
  set(PYTHON_2_COMPATIBILITY "# Using Python 3")
endif()

# execute_process(COMMAND python -c "import mpi4py; print( mpi4py.get_include() )"
#   OUTPUT_VARIABLE mpi4py_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND python -c "import numpy; print( numpy.get_include() )"
  OUTPUT_VARIABLE numpy_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "mpi4py: ${mpi4py_inc}")
message(STATUS "numpy: ${numpy_inc}")

#
# Include Directories
#
set(INCLUDE_DIRS ${MPI_C_INCLUDE_PATH}
                 ${PYTHON_INCLUDE_DIRS}
		 ${mpi4py_inc}
		 ${numpy_inc}
		 ${TIOGA_SOURCE_DIR}/src)

file(GLOB PYTIOGA_SRC *.cpp)

if(BUILD_CUDA_SUPPORT)
  list(APPEND PYTIOGA_DEFINITIONS "-DUSE_CUDA")
endif()

message(STATUS "defs: ${PYTIOGA_DEFINITIONS}")

add_library(pytioga SHARED ${PYTIOGA_SRC})
if(BUILD_CUDA_SUPPORT)
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${PYTHON_LIBRARIES} ${CUDA_LIBRARIES})
else()
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${PYTHON_LIBRARIES})
endif()
target_include_directories(pytioga PUBLIC ${INCLUDE_DIRS})
target_compile_definitions(pytioga PRIVATE ${PYTIOGA_DEFINITIONS})
set_target_properties(pytioga PROPERTIES
  ARCHIVE_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  COMPILE_FLAG "-fPIC"
  PREFIX "" # no lib prefix
  )

