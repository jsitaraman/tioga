#include <Python.h>


template<typename T> void numpy_to_array(PyObject* o, T** x, int *len);
