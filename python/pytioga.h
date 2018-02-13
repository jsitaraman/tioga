#ifndef PYTIOGA_H
#define PYTIOGA_H
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "tioga.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // ignore warnings from numpy
#include <sys/time.h>

#define NTIMERS 5
#define TIMER_EXCHANGE 0
#define TIMER_CONNECT 1
#define TIMER_PROFILE 2

typedef struct {
  PyObject_HEAD
  tioga* tg;
  int initialized;
  int mpi_size, mpi_rank;
  MPI_Comm comm;
  double* xyz;
  int nv;
  int *iblank, *iblankcell;
  int nv2[4],ncell[4];
  int *connptr[4];
  struct timeval tv1, tv2;
  double timers[NTIMERS];
  void tick(){
    gettimeofday(&tv1, NULL);
  }
  double tock(){
    gettimeofday(&tv2, NULL);
    return ( (double)(tv2.tv_usec - tv1.tv_usec)/1000000 +
	     (double)(tv2.tv_sec  - tv1.tv_sec) );
  }
} PyTioga;

void PyTioga_dealloc(PyTioga* self);
int PyTioga_init(PyTioga* self, PyObject* args, PyObject *kwds);
PyObject* PyTioga_dummy(PyTioga* self);
PyObject* PyTioga_register_data(PyTioga* self, PyObject *args);
PyObject* PyTioga_connect(PyTioga* self);
PyObject* PyTioga_update(PyTioga* self, PyObject *args);
PyObject* PyTioga_test_interpolation(PyTioga* self, PyObject* args);


#endif
