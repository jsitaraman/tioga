MODULENAME=tioga
F90= mpif90
#CC = mpicc
CXX= mpicxx
AR = ar -rvs
SWIG = swig -c++ -python
CFLAGS = -std=c++11 -fPIC -O0 -g -rdynamic -Wno-literal-suffix #-Wall -Wextra#-fpe0
FFLAGS = -fPIC  #-CB -traceback #-fbacktrace -fbounds-check
SFLAGS = -I/usr/include/python2.7 -I/usr/local/lib/python2.7/dist-packages/mpi4py/include
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

LDFLAGS= -L/usr/local/intel/10.1.011/fce/lib /usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib  -lifcore  -limf -ldl

lib:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(AR) lib$(MODULENAME).a $(OBJECTS) $(OBJF90)

shared:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) -fPIC -shared -o lib$(MODULENAME).so -lc

swig: $(OBJECTS) $(OBJF90) $(OBJSWIG)
	$(CXX) $(CFLAGS) $(SFLAGS) $(OBJECTS) $(OBJF90) $(OBJSWIG) -shared -fPIC -o _tioga.so -lc

default: $(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -rf *.o lib$(MODULENAME).a lib$(MODULENAME).so

tioga_wrap.o: tioga.i
	$(SWIG) $(SFLAGS) -o tioga_wrap.cpp tioga.i
	$(CXX) -fPIC $(SFLAGS) tioga_wrap.cpp -c -o tioga_wrap.o

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.cpp
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.c
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o
