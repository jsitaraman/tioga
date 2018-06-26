ifeq ($(CONFIG),)
include default.config
else
include $(CONFIG)
endif

MODULENAME=tioga
ifeq ($(strip $(INTEL)),YES)
	F90 = ifort
	CXX = icpc
	CC  = icc
else
	F90 = mpif90
	CXX = mpicxx
	CC  = mpicc
endif

ifeq ($(strip $(CUC)),)
	CUC = nvcc
endif
AR = ar -rvs
LD = $(CXX)
ld=$(CXX)

SWIG = $(SWIG_BIN)/swig

CFLAGS = -std=c++11 -fPIC -rdynamic
CUFLAGS = -std=c++11 --default-stream per-thread -Xcompiler -fPIC
FFLAGS = -fPIC  #-CB -traceback #-fbacktrace -fbounds-check
SFLAGS = -I$(strip $(PYTHON_INC_DIR))/ -I$(strip $(MPI4PY_INC_DIR))/ -I$(strip $(NUMPY_INC_DIR))/

# Intel compiler flags: 
# Floating-point exception; underflow gives 0.0: -fpe0

ifeq ($(strip $(CUDA)),YES)
	CFLAGS += -D_GPU
	CUFLAGS += -D_GPU
	SFLAGS += -D_GPU
	LIBS += -L$(strip $(CUDA_LIB_DIR))/ -lcudart -lcublas -lnvToolsExt -Wl,-rpath=$(strip $(CUDA_LIB_DIR))/
	# If compiling on Ubuntu 16.04 with default GCC, this might be needed:
	CUFLAGS += -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES
ifeq ($(strip $(CU_SM)),)
	CUFLAGS += -arch=sm_20
else
	CUFLAGS += -arch=sm_$(CU_SM)
endif
else
  CFLAGS += -D_CPU
  SFLAGS += -D_CPU
endif
INCS += -I$(strip $(CUDA_INC_DIR))/
INCS += -I$(strip $(MPI_INC_DIR))/

ifeq ($(strip $(DEBUG_LEVEL)),0)
	CFLAGS += -Ofast #-D_NVTX
	CUFLAGS += -O3 -use_fast_math #-D_NVTX
endif
ifeq ($(strip $(DEBUG_LEVEL)),1)
	CFLAGS += -g -O2 #-D_NVTX
	CUFLAGS += -g -O2 #-D_NVTX
endif
ifeq ($(strip $(DEBUG_LEVEL)),2)
	CFLAGS += -g -O0 -D_NVTX
	CUFLAGS += -g -G -O0 -D_NVTX
endif

WARN_ON = -Wall -Wextra -Wconversion
WARN_OFF = -Wno-narrowing -Wno-unused-result -Wno-narrowing -Wno-literal-suffix

ifeq ($(strip $(WARNINGS)),ON)
	CFLAGS += $(WARN_ON)
endif
ifeq ($(strip $(WARNINGS)),OFF)
	CFLAGS += $(WARN_OFF)
	CUFLAGS += -Xcompiler=-Wno-narrowing,-Wno-unused-result,-Wno-narrowing,-Wno-literal-suffix -Xcudafe "--diag_suppress=subscript_out_of_range"
endif

ifeq ($(strip $(OPENMP)),YES)
	CFLAGS += -fopenmp
	CUFLAGS += -Xcompiler=-fopenmp
endif

SRCDIR=$(CURDIR)/src
INCDIR=$(CURDIR)/include
BINDIR=$(CURDIR)/bin

CFLAGS += -I$(INCDIR)/
FFLAGS += -I$(INCDIR)/
CUFLAGS += -I$(INCDIR)/

INCLUDES = codetypes.h MeshBlock.h ADT.h tioga.h globals.h error.hpp points.hpp funcs.hpp
OBJF90 = $(BINDIR)/kaiser.o $(BINDIR)/median.o 
OBJECTS = $(BINDIR)/buildADTrecursion.o $(BINDIR)/searchADTrecursion.o $(BINDIR)/ADT.o\
	$(BINDIR)/MeshBlock.o $(BINDIR)/search.o $(BINDIR)/checkContainment.o $(BINDIR)/bookKeeping.o \
	$(BINDIR)/dataUpdate.o $(BINDIR)/math_funcs.o $(BINDIR)/utils.o $(BINDIR)/linklist.o\
	$(BINDIR)/tioga.o $(BINDIR)/holeMap.o $(BINDIR)/exchangeBoxes.o $(BINDIR)/exchangeSearchData.o $(BINDIR)/exchangeDonors.o\
	$(BINDIR)/parallelComm.o $(BINDIR)/highOrder.o \
	$(BINDIR)/cartOps.o $(BINDIR)/CartGrid.o $(BINDIR)/CartBlock.o $(BINDIR)/getCartReceptors.o $(BINDIR)/get_amr_index_xyz.o\
	$(BINDIR)/exchangeAMRDonors.o $(BINDIR)/funcs.o $(BINDIR)/points.o \
	$(BINDIR)/superMesh.o $(BINDIR)/tiogaInterface.o
OBJSWIG = $(BINDIR)/tioga_wrap.o

ifeq ($(strip $(CUDA)),YES)
	OBJECTS += $(BINDIR)/highOrder_kernels.o $(BINDIR)/dADT.o $(BINDIR)/dMeshBlock.o
endif

ifeq ($(strip $(INTEL)),YES)
	LDFLAGS= -L$(MPI_LIB_DIR)/ -lmpich_intel
	LIBS= -L$(MPI_LIB_DIR)/ -lmpich_intel
else
ifeq ($(CXX),g++)
	# If CXX is not mpicxx, need to manually add MPI library
	LIBS += -L$(MPI_LIB_DIR)/ -lmpi -Wl,-rpath=$(MPI_LIB_DIR)/
endif
endif

lib:	$(OBJECTS) $(OBJF90)
	$(AR) $(BINDIR)/lib$(MODULENAME).a $(OBJECTS) $(OBJF90)

shared:	$(OBJECTS) $(OBJF90)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) -fPIC -shared -o $(BINDIR)/lib$(MODULENAME).so -lc $(LIBS)

swig: CFLAGS += $(SFLAGS)
swig: $(OBJECTS) $(OBJF90) $(OBJSWIG)
	$(CXX) $(CFLAGS) $(SFLAGS) $(OBJECTS) $(OBJF90) $(OBJSWIG) -shared -fPIC -o $(BINDIR)/_tioga.so -lc $(LIBS)

convert: $(BINDIR)/convert_wrap.o
	$(CXX) $(CFLAGS) $(SFLAGS) $(BINDIR)/convert_wrap.o -shared -fPIC -o $(BINDIR)/_convert.so

default: $(OBJECTS) $(OBJF90)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(BINDIR)/$(MODULENAME).exe

$(BINDIR)/%_wrap.cpp: $(INCDIR)/%.i $(INCDIR)/tiogaInterface.h
	@mkdir -p bin
	$(SWIG) -c++ -python $(SFLAGS) -o $@ $<

$(BINDIR)/%_wrap.o: $(BINDIR)/%_wrap.cpp 
	@mkdir -p bin
	$(CXX) $(CFLAGS) $(INCS) $(SFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.cu  $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(CUC) $(INCS) $(CUFLAGS) $(FLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.C $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.c $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.F90 $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(F90) $(FFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.f90 $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(F90) $(FFLAGS) -c $< -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.f $(INCDIR)/*.h $(INCDIR)/*.hpp
	@mkdir -p bin
	$(F90) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean: 
	rm -rf $(BINDIR)/*.o $(BINDIR)/_$(MODULENAME).so $(BINDIR)/$(MODULENAME).py $(BINDIR)/lib$(MODULENAME).a $(BINDIR)/lib$(MODULENAME).so
