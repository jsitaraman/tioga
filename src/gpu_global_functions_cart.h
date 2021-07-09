#include "tioga_gpu.h"
TIOGA_GPU_GLOBAL
void g_interp_data_cart(int *interpList_wcft,
			double *interpList_weights,
			int *interpList_inode,
			int ninterp,
			int nvar_cell,int ncell_nf,
                        int nvar_node,int nnode_nf,
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
    int nweights=(interpList_wcft[i+1]-interpList_wcft[i])/2;
    double *q;
    int offset,nvar,ndof,koffset;
    if (k < nvar_cell) {
        q=qcell;
        offset=0;
        nvar=nvar_cell;
        ndof=ncell_nf;
        koffset=0;
    } else {
        q=qnode;
        offset=nweights;
        nvar=nvar_node;
        ndof=nnode_nf;
        koffset=nvar_cell;
    }
    int m1=interpList_wcft[i]+offset;
    int m2=m1+nweights;
    realData[idx]=0;
    for(int m=m1;m<m2;m++)
      {
	int inode=interpList_inode[m];
	double weight=interpList_weights[m];
	realData[idx]+=q[(k-koffset)*ndof+inode]*weight;
      }
   }
}
