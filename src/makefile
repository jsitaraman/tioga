MODULENAME=tioga
F90= mpif90
CC = mpicc
CXX= mpicxx
AR = ar -rvs
CFLAGS = -fPIC -O2 -rdynamic -g -std=c++11# -g -Wall -Wextra#-fpe0
FFLAGS = -fPIC  #-CB -traceback #-fbacktrace -fbounds-check
INCLUDES = codetypes.h MeshBlock.h ADT.h tioga.h globals.h
OBJF90 = kaiser.o median.o cellVolume.o
OBJECTS = buildADTrecursion.o searchADTrecursion.o ADT.o\
	MeshBlock.o search.o checkContainment.o bookKeeping.o \
	dataUpdate.o math.o utils.o linklist.o\
	tioga.o holeMap.o exchangeBoxes.o exchangeSearchData.o exchangeDonors.o\
	parallelComm.o highOrder.o \
	cartOps.o CartGrid.o CartBlock.o getCartReceptors.o get_amr_index_xyz.o\
	exchangeAMRDonors.o\
	tiogaInterface.o

LDFLAGS= -L/usr/local/intel/10.1.011/fce/lib /usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib  -lifcore  -limf -ldl

lib:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(AR) lib$(MODULENAME).a $(OBJECTS) $(OBJF90)

shared:	$(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC)  -shared -o lib$(MODULENAME).so -lc

default: $(OBJECTS) $(OBJF90) $(INCLUDES)
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -r *.o lib$(MODULENAME).a lib$(MODULENAME).so

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o
