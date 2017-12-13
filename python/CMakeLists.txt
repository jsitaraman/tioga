
find_package(PythonLibs REQUIRED)
execute_process(COMMAND python -c "import mpi4py; print( mpi4py.get_include() )"
  OUTPUT_VARIABLE mpi4py_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "mpi4py: ${mpi4py_inc}")

#
# Include Directories
#
set(INCLUDE_DIRS ${MPI_C_INCLUDE_PATH}
                 ${PYTHON_INCLUDE_DIRS}
                 ${mpi4py_inc}
		 ${TIOGA_SOURCE_DIR}/src)

file(GLOB PYTIOGA_SRC *.cpp)

add_library(pytioga SHARED ${PYTIOGA_SRC})
target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${PYTHON_LIBRARIES})
target_include_directories(pytioga PUBLIC ${INCLUDE_DIRS})
set_target_properties(pytioga PROPERTIES
  ARCHIVE_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  COMPILE_FLAG "-fPIC"
  PREFIX "" # no lib prefix
  )

