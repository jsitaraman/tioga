# Generated automatically from Makefile.in by configure.
ALL: default
##### User configurable options #####

SHELL       = /bin/sh
ARCH        = LINUX
CC          = gcc
CLINKER     = $(CC)
F77         = ~/tests/openmpi/gnu/bin/mpif90
FLINKER     = $(F77)
F90         = ~/tests/openmpi/gnu/bin/mpif90
F90LINKER   = $(F90)      
CCC         = gcc
CCLINKER    = $(CCC)

MAKE        = make --no-print-directory
#
srcdir= $(PWD)

### End User configurable options ###

.SUFFIXES:	.cc .C .f90
CFLAGS	  = $(OPTFLAGS) 
CCFLAGS	  = $(CFLAGS)
#FFLAGS	  = $(OPTFLAGS) -r8 -traceback -CB
FFLAGS	  = $(OPTFLAGS) -fdefault-real-8 -ffree-line-length-132 -fbounds-check -g
LDFLAGS   =  
LIBS	  = -lm
ALL_EXECS = buildGrid 

ALL:	 buildGrid

default: $(EXECS)

OBJS    = modules.o buildGrid.o clust.o partSphere.o findfact.o
buildGrid: $(OBJS)
	$(FLINKER) $(FFLAGS) -o buildGrid $(OBJS) $(LDFLAGS)

clean:
	rm -f work.pc work.pcl *.mod
	rm -f *.o *~ PI* $(ALL_EXECS) upshot rdb.* startup.* core

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(F77) $(FFLAGS) -c $<
.f90.o:
	$(F90) -c $(FFLAGS) $<
.f:
	$(F77) -c $<
	$(FLINKER) -o $* $*.o
.c:	
	$(CC) -c $<
	$(CLINKER) -o $* $*.o
