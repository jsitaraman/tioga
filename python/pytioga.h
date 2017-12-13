#ifndef PYTIOGA_H
#define PYTIOGA_H
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "tioga.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // ignore warnings from numpy

typedef struct {
  PyObject_HEAD
  tioga* tg;
  int initialized;
  int mpi_size, mpi_rank;
  MPI_Comm comm;
  int nv2[4],ncell[4];
  int *connptr[4];
} PyTioga;

void PyTioga_dealloc(PyTioga* self);
int PyTioga_init(PyTioga* self, PyObject* args, PyObject *kwds);
PyObject* PyTioga_dummy(PyTioga* self);
PyObject* PyTioga_register_data(PyTioga* self, PyObject *args);
PyObject* PyTioga_preprocess(PyTioga* self);
PyObject* PyTioga_connect(PyTioga* self);
PyObject* PyTioga_update(PyTioga* self, PyObject *args);

// class PyTioga {

//   tioga *tg;

//   Py_buffer b_coord;
//   Py_buffer b_wallnode;
//   Py_buffer b_obcnode;
//   Py_buffer b_iblank;
//   Py_buffer b_tetConn;
//   Py_buffer b_pyraConn;
//   Py_buffer b_prismConn;
//   Py_buffer b_hexaConn;
//   Py_buffer b_qvars;

//   int initialized;

//   int mpi_size, mpi_rank;
//   int nv2[4],ncell[4];
//   int *connptr[4];
//   int *ndc4, *ndc5, *ndc6, *ndc8;

//   MPI_Comm comm;

//  public:
//   PyTioga();
//   PyTioga(boost::python::object pycomm);
//   ~PyTioga();

//   void setup(MPI_Comm comm);
//   void register_data(boost::python::dict data);
//   void preprocess_grids();
//   void perform_connectivity();
//   void register_solution(int block_id, boost::python::object qo);
//   void data_update(int nq);
//   void write_data(int nq);

// };

#endif
