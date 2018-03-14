#include "pytioga.h"
#include "tioga.h"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL tioga_ARRAY_API
#include <numpy/ndarrayobject.h>
#include "utils.h"

void PyTioga_dealloc(PyTioga* self){
  //
  // Release hold of the python data buffers
  //

  printf("# TIOGA timers: %12.4f s (exchange) %12.4f s (connect) %12.4f s (profile)\n", 
	 self->timers[TIMER_EXCHANGE], self->timers[TIMER_CONNECT], self->timers[TIMER_PROFILE]);

  delete [] self->tg;

  Py_TYPE(self)->tp_free((PyObject*)self);
}

int PyTioga_init(PyTioga* self, PyObject* args, PyObject *kwds){

  int size, rank;
  int fcom;
  MPI_Comm comm;

  // we haven't done anything yet!
  self->initialized = 0;

  // Parse the argument as an int
  if(!PyArg_ParseTuple(args, "i", &fcom)){
    printf("Missing MPI Comm arg");
    return 1;
  }

  // note the use of "self" instead of "this"
  self->comm = MPI_Comm_f2c((MPI_Fint)fcom);
  
  MPI_Comm_size(self->comm, &size);
  MPI_Comm_rank(self->comm, &rank);

  self->mpi_size = size;
  self->mpi_rank = rank;

  self->tg = new tioga[1];

  self->iblank     = NULL;
  self->iblankcell = NULL;

  self->tg->setCommunicator(self->comm,rank,size);

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

  int btag, bid, iblk, nwbc, nobc;
  int *wbcnode, *obcnode;
  int *tmpbtag;
  int *ndc4, *ndc5, *ndc6, *ndc8;
  int ntypes, count4, count5, count6, count8;
  int dbg, nq, len, dummy;

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
  numpy_to_array(data, &(self->xyz), &len);

  if(PyDict_Contains(dict,Py_BuildValue("s", "iblanking-cell"))){
    data = PyList_GetItem(PyDict_GetItemString(dict, "iblanking-cell"), 0);
    numpy_to_array(data, &(self->iblankcell), &dummy);
  }

  data = PyList_GetItem(PyDict_GetItemString(dict, "iblanking"), 0);
  numpy_to_array(data, &(self->iblank), &(self->nv));

  if(self->nv != len/3){
    printf("Uh oh... nv != len(xyz)/3 \n");
  }

  data = PyList_GetItem(PyDict_GetItemString(dict, "bodyTag"), 0);
  numpy_to_array(data, &tmpbtag, &len);
  btag = tmpbtag[0];

  // printf("# pytioga (%2d) %2d nv = %5d\n", self->mpi_rank, btag,  nv);

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
    self->tg->registerGridData(btag,self->nv,self->xyz,self->iblank,nwbc,nobc,
			       wbcnode,obcnode,ntypes,self->nv2,self->ncell,self->connptr);
  } 

  return Py_BuildValue("i", 0);

}

PyObject* PyTioga_connect(PyTioga* self){

  self->tick();

  self->tg->profile();

  self->timers[TIMER_PROFILE] += self->tock();

  self->tick();

  self->tg->performConnectivity();

  self->timers[TIMER_CONNECT] += self->tock();

  if(self->iblankcell != NULL){
    self->tg->getiBlankCell(self->iblankcell);
  }

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

PyObject* PyTioga_test_interpolation(PyTioga* self, PyObject* args){

  int j, k, aux, nvartmp=1;
  int nv = self->nv;
  double *xyz = self->xyz;
  int *iblank_node = self->iblank;
  double l2norm, l2old;

  double *qtmp = (double *)malloc(sizeof(double)*1*nvartmp*nv);

  int myid, nproc;

  int count8  = self->ncell[0];
  
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  for (j = 0;  j < nv; j++){
    if(iblank_node[j] == 1){
      for (k = 0; k < nvartmp; k++)
	qtmp[nvartmp*j+k] = xyz[3*j]+xyz[3*j+1]+xyz[3*j+2];
    } else {
      for (k = 0; k < nvartmp; k++)
	qtmp[nvartmp*j+k] = -100.0;
    }
  }

  self->tg->dataUpdate(nvartmp,qtmp,0);

  //
  l2norm=0.;
  l2old=0.;
  aux = 1;
  for (j = 0;  j < nv; j++){
    if(iblank_node[j] == -1){
      for (k = 0; k < nvartmp; k++){
	if(qtmp[nvartmp*j+k] != qtmp[nvartmp*j+k]){
	  printf("pytioga.cpp: qtmp Nanned out..\n");
	  printf("Proc %d, k %d, cellID %d\n",myid,k,j);
	  MPI_Abort(MPI_COMM_WORLD,-124);
	}

	l2norm += ((qtmp[nvartmp*j+k] - (xyz[3*j]+xyz[3*j+1]+xyz[3*j+2]))*
		   (qtmp[nvartmp*j+k] - (xyz[3*j]+xyz[3*j+1]+xyz[3*j+2])));
         
	// errorcheck similar to cell based formulation but since
	// l2norm is added, ony change of norm at certain position is
	// relevant to find bad spots
	if((l2norm-l2old)>1.e-10){
	  printf("# l2norm is not zero: %e\n",l2norm);
	  printf("proc ID is [%d] \n",myid);
	  printf("node position : %e %e %e\n",xyz[3*j],xyz[3*j+1],xyz[3*j+2]);
	  printf("q value : %e\n",xyz[3*j]+xyz[3*j+1]+xyz[3*j+2]);
	  printf("q computed : %e\n\n",qtmp[nvartmp*j+k]);
	  l2old = l2norm;
	}
      }
    }
  }

  //
  // find maximum interpolation error from all processors
  //
  // MPI_Gather to proc 0
  //
  double *interp_error = (double *)malloc(sizeof(double)*nproc);
  MPI_Gather(&l2norm,1,MPI_DOUBLE,interp_error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //
  if(myid==0){
    int id_max_interp_error;
    double max_interp_error=0.0;
    for(j=0;j<nproc;j++){
      max_interp_error = max(max_interp_error,interp_error[j]);
      if(interp_error[j] > 1.0e-10){
        printf("\n interpol. error per proc= %e at proc %d \n", interp_error[j],j);
      }
    }
    
    printf("\n#pyTIOGA: ihigh [%d]. Maximum interpolation error: %e \n\n",0,max_interp_error);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  free(qtmp);

  return Py_BuildValue("i", 0);
}

