#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // ignore warnings from numpy
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL tioga_ARRAY_API
#include <numpy/ndarrayobject.h>

template<typename T>
void numpy_to_array(PyObject* o, T** x, int *len){
  PyArrayObject *arr;
  arr = (PyArrayObject*)o;

  int ndim;
  npy_intp *dims;

  len[0] = 1;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  for(int i=0; i<ndim; i++){
    len[0] *= dims[i];
  }

  x[0] = (T*)PyArray_DATA(arr);

}

template void numpy_to_array<int>(PyObject* o, int** x, int *len);
template void numpy_to_array<double>(PyObject* o, double** x, int *len);
