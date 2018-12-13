MODULENAME=tioga
F90= mpif90
CC = mpicc
CXX= mpicxx
AR = ar -rvs
FFLAGS=-fdefault-real-8 -fbounds-check -fbacktrace -g #-CB -traceback -r8
CFLAGS = -g -rdynamic #-fpe0
INCLUDES = 
OBJEXEC=
OBJECTS=
OBJF90 = modules.o gridutils.o testTioga.o
TIOGALIB=../build/libtioga.a

#LDFLAGS= -L../src/ -L/usr/local/intel/10.1.011/fce/lib -L/usr/local/openmpi/openmpi-1.4.3/x86_64/ib/intel10/lib  -lifcore  -limf -ldl  -ltioga -lstdc++

LDFLAGS= -L../build/ -ltioga -lstdc++ #-lmpi_cxx

default: $(OBJF90) $(LIB)
	$(F90) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -rf *.o *.mod *.exe 2>/dev/null

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
