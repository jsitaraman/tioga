%module tioga

// -----------------------------------------------------------------------------
// Header files required by any of the following C++ code
// -----------------------------------------------------------------------------
%header
%{
#include <mpi.h>
#include "tiogaInterface.h"
%}

// -----------------------------------------------------------------------------
// Header files and other declarations to be parsed as SWIG input
// -----------------------------------------------------------------------------

// SWIG interface file for MPI, and typemap for MPI_Comm to Python Comm
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm,MPI_Comm);

%pythoncallback;
// Functions declared here will be able to act like (C-style) function pointers
void tioga_dataupdate_ab(int nvar, double* q_spts, int gradFlag);
void tioga_preprocess_grids_(void);
void tioga_performconnectivity_(void);
%nopythoncallback;

%ignore tioga_dataupdate_ab;
%ignore tioga_preprocess_grids_;
%ignore tioga_performconnectivity_;
%include "tiogaInterface.h"

// <-- Additional C++ declations [anything that would normally go in a header]

// -----------------------------------------------------------------------------
// Additional functions which have been declared, but not defined (including
// definition in other source files which will be linked in later)
// -----------------------------------------------------------------------------

%inline
%{
// <-- Additional C++ definitions [anything that would normally go in a .cpp]
%}

// -----------------------------------------------------------------------------
// Additional Python functions to add to module
// [can use any functions/variables declared above]
// -----------------------------------------------------------------------------

%pythoncode
%{
# Python functions here
%}
