#include "pytioga.h"
#include "tioga.h"
#include <mpi4py/mpi4py.h>
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL tioga_ARRAY_API
#include <numpy/ndarrayobject.h>
#include "utils.h"

void PyTioga_dealloc(PyTioga* self){
  //
  // Release hold of the python data buffers
  //

  printf("# TIOGA timers: %12.4f s (exchange) %12.f s (connect)", 
	 self->timers[TIMER_EXCHANGE], self->timers[TIMER_CONNECT]);

  delete [] self->tg;

  Py_TYPE(self)->tp_free((PyObject*)self);
}

int PyTioga_init(PyTioga* self, PyObject* args, PyObject *kwds){

  PyObject* pycom;
  int size, rank;
  MPI_Comm *comm;

  // we haven't done anything yet!
  self->initialized = 0;

  // Check mpi4py is loaded
  if(import_mpi4py()<0) return 0;

  // Parse the argument as a PyObject ("O")
  if(!PyArg_ParseTuple(args, "O", &pycom)){
    printf("Missing MPI Comm arg");
    return 1;
  }

  // Try to translate this into an MPI_Comm*
  comm = PyMPIComm_Get(pycom);
  if(comm == NULL){
    printf("Error with the python mpi communicator\n");
    return 1;
  }

  MPI_Comm_size(comm[0], &size);
  MPI_Comm_rank(comm[0], &rank);

  // note the use of "self" instead of "this"
  self->comm     = comm[0];
  self->mpi_size = size;
  self->mpi_rank = rank;

  self->tg = new tioga[1];

  self->tg->setCommunicator(comm[0],rank,size);

  self->timers[TIMER_CONNECT]  = 0.0;
  self->timers[TIMER_EXCHANGE] = 0.0;

  return 0;

} // end __init__

PyObject* PyTioga_dummy(PyTioga* self){

  printf("PyTioga dummy function\n");

  Py_INCREF(Py_None);
  return Py_None; // every python function returns something, even if it's "None"
}

bool check_python_data(PyObject* d){
  if(not PyDict_Contains(d,Py_BuildValue("s", "tetConn"         ))){ printf("missing tetConn key\n");          return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "pyraConn"        ))){ printf("missing pyraConn key\n");         return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "prismConn"       ))){ printf("missing prismConn key\n");        return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "hexaConn"        ))){ printf("missing hexaConn key\n");         return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "wallnode"        ))){ printf("missing wallnode key\n");         return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "obcnode"         ))){ printf("missing obcnode key\n");          return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "iblanking"       ))){ printf("missing iblanking key\n");        return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "bodyTag"         ))){ printf("missing bodyTag key\n");          return 1;}
  if(not PyDict_Contains(d,Py_BuildValue("s", "grid-coordinates"))){ printf("missing grid-coordinates key\n"); return 1;}
  return 0;
}

PyObject* PyTioga_register_data(PyTioga* self, PyObject *args){

  PyObject *dict;
  PyObject *list, *data; // temporary pointers for unpacking the dictionary

  int btag, bid, iblk, nv, nwbc, nobc;
  double* xyz;
  int *wbcnode, *obcnode, *iblank;
  int *tmpbtag;
  int *ndc4, *ndc5, *ndc6, *ndc8;
  int ntypes, count4, count5, count6, count8;
  int dbg, nq, len;

  if(!PyArg_ParseTuple(args, "O", &dict)){
    printf("Missing argument\n");
    return Py_BuildValue("i", 1);
  }

  if(!PyDict_Check(dict)){
    printf("That's not a dictionary!\n");
    return Py_BuildValue("i", 1);
  }

  if( check_python_data(dict) ){
    return Py_BuildValue("i", 1);
  }

  // With the python-c api, most things are "borrowed" references
  // meaning you don't have to decrement the pyobject reference when
  // you're done. Getting items from dictionary and lists are this
  // type of operation so that saves some trouble. When SETTING a list
  // item, though, you have to be careful. See here:
  // https://docs.python.org/3/extending/extending.html#ownership-rules

  data = PyList_GetItem(PyDict_GetItemString(dict, "tetConn"), 0);
  numpy_to_array(data, &ndc4, &count4);

  data = PyList_GetItem(PyDict_GetItemString(dict, "pyraConn"), 0);
  numpy_to_array(data, &ndc5, &count5);

  data = PyList_GetItem(PyDict_GetItemString(dict, "prismConn"), 0);
  numpy_to_array(data, &ndc6, &count6);

  data = PyList_GetItem(PyDict_GetItemString(dict, "hexaConn"), 0);
  numpy_to_array(data, &ndc8, &count8);

  data = PyList_GetItem(PyDict_GetItemString(dict, "wallnode"), 0);
  numpy_to_array(data, &wbcnode, &nwbc);

  data = PyList_GetItem(PyDict_GetItemString(dict, "obcnode"), 0);
  numpy_to_array(data, &obcnode, &nobc);

  data = PyList_GetItem(PyDict_GetItemString(dict, "grid-coordinates"), 0);
  numpy_to_array(data, &xyz, &len);

  data = PyList_GetItem(PyDict_GetItemString(dict, "iblanking"), 0);
  numpy_to_array(data, &iblank, &nv);

  if(nv != len/3){
    printf("Uh oh... nv != len(xyz)/3 \n");
  }

  data = PyList_GetItem(PyDict_GetItemString(dict, "bodyTag"), 0);
  numpy_to_array(data, &tmpbtag, &len);
  btag = tmpbtag[0];

  // Now we'll look through all the connectivity data and see which
  // one are / aren't used. Keep track of the number of types of
  // connectivity.
  ntypes = 0;
  
  // 4-vertex Cells
  if(count4 > 0){
    self->ncell[ntypes]   = count4/4;
    self->nv2[ntypes]     = 4;
    self->connptr[ntypes] = ndc4;
    ntypes++;
  }
  // 5-vertex Cells
  if(count5 > 0){
    self->ncell[ntypes]   = count5/5;
    self->nv2[ntypes]     = 5;
    self->connptr[ntypes] = ndc5;
    ntypes++;
  }
  // 6-vertex Cells
  if(count6 > 0){
    self->ncell[ntypes]   = count6/6;
    self->nv2[ntypes]     = 6;
    self->connptr[ntypes] = ndc6;
    ntypes++;
  }
  // 8-vertex Cells
  if(count8 > 0){
    self->ncell[ntypes]   = count8/8;
    self->nv2[ntypes]     = 8;
    self->connptr[ntypes] = ndc8;
    ntypes++;
  }

  if(ntypes == 0){
    printf("Error: using tioga without providing data\n");
    return Py_BuildValue("i", 1);
  } else {
    self->tg->registerGridData(btag,nv,xyz,iblank,nwbc,nobc,
			       wbcnode,obcnode,ntypes,self->nv2,self->ncell,self->connptr);
  } 

  return Py_BuildValue("i", 0);

}

PyObject* PyTioga_connect(PyTioga* self){

  self->tick();

  self->tg->profile();
  self->tg->performConnectivity();

  self->timers[TIMER_CONNECT] += self->tock();

  return Py_BuildValue("i", 0);
}

PyObject* PyTioga_update(PyTioga* self, PyObject* args){

  PyObject *data;
  int nq, len;
  int interptype = 0; // only row implemented here
  double *q;

  self->tick();

  if(!PyArg_ParseTuple(args, "Oi", &data, &nq)){
    printf("Missing argument\n");
    return Py_BuildValue("i", 1);
  }

  numpy_to_array(data, &q, &len);

  self->tg->dataUpdate(nq,q,interptype);

  self->timers[TIMER_EXCHANGE] += self->tock();

  return Py_BuildValue("i", 0);

}
