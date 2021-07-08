#include "tioga_gpu.h"
TIOGA_GPU_GLOBAL
void g_interp_data_cart(int *interpList_wcft,
			double *interpList_weights,
			int *interpList_inode,
			int ninterp,
			int nvar_cell,int nvar_node,
			double *realData,
			double *qcell,
			double *qnode)
{
  int nsend=(nvar_cell+nvar_node)*ninterp;
#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nsend) 
#else
  for(int idx=0;idx<nsend;idx++) 
#endif
   {
    int i=idx/(nvar_cell+nvar_node);
    int k=idx-i*(nvar_cell+nvar_node);
    realData[idx]=0;
    int nweights=interpList_wcft[i+1]-interpList_wcft[i];
    double *q=(k<nvar_cell) ? qcell: qnode;
    int offset=(k < nvar_cell) ? 0 : nweights;
    int nvar=(k < nvar_cell ) ? nvar_cell: nvar_node;
    int m1=interpList_wcft[i]+offset;
    int m2=m1+nweights;
    for(int m=m1;m<m2;m++)
      {
	int inode=interpList_inode[m];
	double weight=interpList_weights[m];
	realData[idx]+=q[inode*nvar+k]*weight;
      }
   }
}
