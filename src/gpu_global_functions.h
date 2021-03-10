#include "tioga_gpu.h"
#include "gpu_device_functions.h"
TIOGA_GPU_GLOBAL
void g_reset_iblanks(int *iblank, int nnodes)
{
#ifdef TIOGA_HAS_GPU
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nnodes) {
#else
  for(int idx=0;idx<nnodes;idx++) {
#endif
    iblank[idx]=1;
  }

}

TIOGA_GPU_GLOBAL
void g_adt_search(double *x, int **vconn,int *nc, int *nv, int ntypes,
             double *coord, double *adtExtents, int *adtIntegers, double *adtReals,
             int *elementList, int *donorId, double *xsearch, 
	     int ndim, int nelem, int nsearch)
{
#ifdef TIOGA_HAS_GPU
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < nsearch) {
#else
  for (int idx =0;idx<nsearch;idx++) {
#endif
    int cellIndex[2];
    int nchecks;
    donorId[idx]=-1;
    d_searchIntersections_containment(cellIndex,
                                     adtIntegers,
                                     adtReals,
                                     coord,
                                     &(xsearch[3*idx]),
                                     x, nc, nv, ntypes,
                                     elementList,
                                     vconn,
                                     nelem,
                                     ndim,
                                     &nchecks);
    if (cellIndex[0] > -1 && cellIndex[1]==0) donorId[idx]=cellIndex[0];
  }
}
