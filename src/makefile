ifeq ($(CONFIG),)
include default.config
else
include $(CONFIG)
endif

MODULENAME=tioga
F90= mpif90
#CC = mpicc
CXX= mpicxx
ifeq ($(strip $(CUC)),)
	CUC=nvcc
endif
AR = ar -rvs

SWIG = $(SWIG_BIN)/swig -c++ -python

CFLAGS = -std=c++11 -g -O3 -fPIC -rdynamic 
CUFLAGS = -std=c++11 --default-stream per-thread -Xcompiler -fPIC
FFLAGS = -fPIC  #-CB -traceback #-fbacktrace -fbounds-check
SFLAGS = -I$(PYTHON_INC_DIR)/ -I$(MPI4PY_INC_DIR)/

# Intel compiler flags: 
# Floating-point exception; underflow gives 0.0: -fpe0

ifeq ($(strip $(CUDA)),YES)
	CFLAGS += -D_GPU
	CUFLAGS += -D_GPU
	LIBS += -L$(strip $(CUDA_DIR))/lib64 -lcudart -lcublas -lnvToolsExt -Wl,-rpath=$(strip $(CUDA_DIR))/lib64
	INCS += -I$(strip $(CUDA_DIR))/include
ifeq ($(strip $(CU_SM)),)
	CUFLAGS += -arch=sm_20
else
	CUFLAGS += -arch=sm_$(CU_SM)
endif
endif

ifeq ($(strip $(DEBUG_LEVEL)),0)
	CFLAGS += -Ofast
	CUFLAGS += -O3 -use_fast_math
endif
ifeq ($(strip $(DEBUG_LEVEL)),1)
	CFLAGS += -g -O3 -D_NVTX
	CUFLAGS += -g -O3 -D_NVTX
endif
ifeq ($(strip $(DEBUG_LEVEL)),2)
	CFLAGS += -g -O0 -D_NVTX
	CUFLAGS += -g -O0 -D_NVTX
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

INCLUDES = codetypes.h MeshBlock.h ADT.h tioga.h globals.h error.hpp points.hpp funcs.hpp
OBJF90 = kaiser.o median.o 
OBJECTS = buildADTrecursion.o searchADTrecursion.o ADT.o\
	MeshBlock.o search.o checkContainment.o bookKeeping.o \
	dataUpdate.o math.o utils.o linklist.o\
	tioga.o holeMap.o exchangeBoxes.o exchangeSearchData.o exchangeDonors.o\
	parallelComm.o highOrder.o \
	cartOps.o CartGrid.o CartBlock.o getCartReceptors.o get_amr_index_xyz.o\
	exchangeAMRDonors.o funcs.o points.o \
	tiogaInterface.o
OBJSWIG = tioga_wrap.o

ifeq ($(strip $(CUDA)),YES)
	OBJECTS += highOrder_kernels.o
endif

LDFLAGS= -L/usr/local/intel/10.1.011/fce/lib /usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib  -lifcore  -limf -ldl

lib:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(AR) lib$(MODULENAME).a $(OBJECTS) $(OBJF90)

shared:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) -fPIC -shared -o lib$(MODULENAME).so -lc $(LIBS)

swig: CFLAGS += $(SFLAGS)
swig: $(OBJECTS) $(OBJF90) $(OBJSWIG)
	$(CXX) $(CFLAGS) $(SFLAGS) $(OBJECTS) $(OBJF90) $(OBJSWIG) -shared -fPIC -o _tioga.so -lc $(LIBS)

default: $(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -rf *.o _$(MODULENAME).so $(MODULENAME).py lib$(MODULENAME).a lib$(MODULENAME).so

%_wrap.cpp: %.i tiogaInterface.h
	$(SWIG) $(SFLAGS) -o $@ $<

%.o:%.cu
	$(CUC) $(INCS) $(CUFLAGS) $(FLAGS) -c $< -o $*.o
%.o:%.cpp
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $*.o
%.o:%.c
	$(CXX) $(INCS) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o
