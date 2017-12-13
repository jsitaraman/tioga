#include "pytioga.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL tioga_ARRAY_API
#include <numpy/ndarrayobject.h>


//
// Python wrapping using the native Python-C API. This is tedious but
// it's only becuase it exposes all possible python functionality.
// Things like the PyTypeObject declaration can just be copy pasted
// from this template.
//
// Online Documentation
// Part 1 (basic extending)    https://docs.python.org/2/extending/extending.html
// Part 2 (structs)            https://docs.python.org/2/extending/newtypes.html
//
// For clarity, there should be module name and a python object name
// that are NOT the same. Here the module name is "pytioga" (imported
// from a "pytioga.so") and the class name is "PyTioga". So in python
// we would do the following:
//
// import pytioga                      
// connectivity = pytioga.PyTioga(...) # <-- PyTioga instantiation

// -----------------------------------------------------------------------------
//
// Methods for the PyTioga struct (Think class methods for the PyTioga
// class). METH_NOARGS means there are no arguments. METH_VARARGS
// means there are arguments passed from python (and have to be
// parsed, see pytioga.cpp for example).
//
static PyMethodDef PyTioga_methods[] = {
  {"dummy",            (PyCFunction)PyTioga_dummy,    METH_NOARGS,  "Dummy"},
  {"register_data",(PyCFunction)PyTioga_register_data,METH_VARARGS, "Register Data"},
  {"preprocess",   (PyCFunction)PyTioga_preprocess,   METH_NOARGS,  "Pre-process Grids"},
  {"connect",      (PyCFunction)PyTioga_connect,      METH_NOARGS,  "Connect Grids"},
  {"update",       (PyCFunction)PyTioga_update,       METH_VARARGS, "Update Data"},
  {NULL}  /* Sentinel */
};
// -----------------------------------------------------------------------------


// Methods for the pytioga MODULE (none)
static PyMethodDef pytioga_methods[] = {
  {NULL}  /* Sentinel */
};

// Declare a new python type with our methods
static PyTypeObject pytioga_PyTiogaType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "pytioga.PyTioga",           // tp_name 
  sizeof(PyTioga),             // tp_basicsize 
  0,                           // tp_itemsize 
  (destructor)PyTioga_dealloc, // tp_dealloc 
  0,                           // tp_print 
  0,                           // tp_getattr 
  0,                           // tp_setattr 
  0,                           // tp_compare 
  0,                           // tp_repr 
  0,                           // tp_as_number 
  0,                           // tp_as_sequence 
  0,                           // tp_as_mapping 
  0,                           // tp_hash 
  0,                           // tp_call 
  0,                           // tp_str 
  0,                           // tp_getattro 
  0,                           // tp_setattro 
  0,                           // tp_as_buffer 
  Py_TPFLAGS_DEFAULT |
  Py_TPFLAGS_BASETYPE,         // tp_flags 
  "PyTioga objects",           // tp_doc 
  0,                           // tp_traverse 
  0,                           // tp_clear 
  0,                           // tp_richcompare 
  0,                           // tp_weaklistoffset 
  0,                           // tp_iter 
  0,                           // tp_iternext 
  PyTioga_methods,             // tp_methods 
  0,                           // tp_members (no members for simplicity) 
  0,                           // tp_getset 
  0,                           // tp_base 
  0,                           // tp_dict 
  0,                           // tp_descr_get 
  0,                           // tp_descr_set 
  0,                           // tp_dictoffset 
  (initproc)PyTioga_init,      // tp_init 
  0,                           // tp_alloc 
};


// We're using c++ and python can only see C init functions
extern "C" {

#ifndef PyMODINIT_FUNC/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
// In python if we want to do "import pytioga", we have to declare a
// "void initpytioga()" function. This would be for a library called
// "pytioga.so". The naming here matters so make sure to follow this
// format. CMake likes to add "lib" as a prefix to libraries but you
// can tell it not to.
PyMODINIT_FUNC
initpytioga(void){

  PyObject* m;

  import_array(); // numpy stuff

  //
  // make sure we set the generic "new" function (cannot be done in
  // struct above). The "new" function is called before "init" and is
  // not relevant to what we're doing.
  pytioga_PyTiogaType.tp_new = PyType_GenericNew;

  // check PyTioga was created properly
  if (PyType_Ready(&pytioga_PyTiogaType) < 0)
    return;

  // Initialize the tioga module
  m = Py_InitModule3("pytioga", pytioga_methods,
  		     "Python Tioga module including PyTioga object");

  // Tell python we're using this type, dont clean it up from memory
  Py_INCREF(&pytioga_PyTiogaType);
  // Add our type to the module and connect the functions we wrote
  PyModule_AddObject(m, "PyTioga", (PyObject *)&pytioga_PyTiogaType);
 }

}


//// Old Boost Stuff:
//   class_<PyTioga>("PyTioga", init<object>())
//     .def(init<>())
//     .def("register_data", &PyTioga::register_data)
//     .def("preprocess_grids", &PyTioga::preprocess_grids)
//     .def("perform_connectivity", &PyTioga::perform_connectivity)
//     .def("register_solution", &PyTioga::register_solution)
//     .def("write_data", &PyTioga::write_data)
//     .def("data_update", &PyTioga::data_update);
