%module convert

// -----------------------------------------------------------------------------
// Header files required by any of the following C++ code
// -----------------------------------------------------------------------------
%header
%{
#include <numpy/arrayobject.h>
%}

%init
%{
  import_array();
%}

// -----------------------------------------------------------------------------
// Header files and other declarations to be parsed as SWIG input
// -----------------------------------------------------------------------------

// <-- Additional C++ declations [anything that would normally go in a header]

// -----------------------------------------------------------------------------
// Additional functions which have been declared, but not defined (including
// definition in other source files which will be linked in later)
// -----------------------------------------------------------------------------

// --------------------------------------------------------
// FUNCTIONS TO CONVERT POINTERS TO NUMPY ARRAYS (AND BACK)
// --------------------------------------------------------
%inline
%{
/* ------------ 1D Arrays ------------ */

  PyObject* ptrToArray(float* data, int n)
{
  npy_intp dims[1] = {n};
  return PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,(void*)data);
}

PyObject* ptrToArray(double* data, int n)
{
  npy_intp dims[1] = {n};
  return PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,(void*)data);
}

PyObject* ptrToArray(int* data, int n)
{
  npy_intp dims[1] = {n};
  return PyArray_SimpleNewFromData(1,dims,NPY_INT,(void*)data);
}

PyObject* ptrToArray(unsigned int* data, int n)
{
  npy_intp dims[1] = {n};
  return PyArray_SimpleNewFromData(1,dims,NPY_UINT,(void*)data);
}

/* ------------ 2D Arrays ------------ */

PyObject* ptrToArray(float* data, int n1, int n2)
{
  npy_intp dims[2] = {n1, n2};
  return PyArray_SimpleNewFromData(2,dims,NPY_FLOAT,(void*)data);
}

PyObject* ptrToArray(double* data, int n1, int n2)
{
  npy_intp dims[2] = {n1, n2};
  return PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,(void*)data);
}

PyObject* ptrToArray(int* data, int n1, int n2)
{
  npy_intp dims[2] = {n1, n2};
  return PyArray_SimpleNewFromData(2,dims,NPY_INT,(void*)data);
}

PyObject* ptrToArray(unsigned int* data, int n1, int n2)
{
  npy_intp dims[2] = {n1, n2};
  return PyArray_SimpleNewFromData(2,dims,NPY_UINT,(void*)data);
}

/* ------------ 3D Arrays ------------ */

PyObject* ptrToArray(float* data, int n1, int n2, int n3)
{
  npy_intp dims[3] = {n1, n2, n3};
  return PyArray_SimpleNewFromData(3,dims,NPY_FLOAT,(void*)data);
}

PyObject* ptrToArray(double* data, int n1, int n2, int n3)
{
  npy_intp dims[3] = {n1, n2, n3};
  return PyArray_SimpleNewFromData(3,dims,NPY_DOUBLE,(void*)data);
}

PyObject* ptrToArray(int* data, int n1, int n2, int n3)
{
  npy_intp dims[3] = {n1, n2, n3};
  return PyArray_SimpleNewFromData(3,dims,NPY_INT,(void*)data);
}

PyObject* ptrToArray(unsigned int* data, int n1, int n2, int n3)
{
  npy_intp dims[3] = {n1, n2, n3};
  return PyArray_SimpleNewFromData(3,dims,NPY_UINT,(void*)data);
}

/* ------------ 4D Arrays ------------ */

PyObject* ptrToArray(double* data, int n1, int n2, int n3, int n4)
{
  npy_intp dims[4] = {n1, n2, n3, n4};
  return PyArray_SimpleNewFromData(4,dims,NPY_DOUBLE,(void*)data);
}

/* ------------ Get Pointer from Numpy Array ------------ */

double* arrayToDblPtr(PyObject* arr)
{
  return (double *)(((PyArrayObject *)arr)->data);
}

float* arrayToFloatPtr(PyObject* arr)
{
  return (float *)(((PyArrayObject *)arr)->data);
}

int* arrayToIntPtr(PyObject* arr)
{
  return (int *)(((PyArrayObject *)arr)->data);
}

unsigned int* arrayToUintPtr(PyObject* arr)
{
  return (unsigned int *)(((PyArrayObject *)arr)->data);
}
%}


// -----------------------------------------------------------------------------
// Additional Python functions to add to module
// [can use any functions/variables declared above]
// -----------------------------------------------------------------------------

%pythoncode
%{

%}

