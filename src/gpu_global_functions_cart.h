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
    int offset,ndof,koffset;
    if (k < nvar_cell) {
        q=qcell;
        offset=0;
        ndof=ncell_nf;
        koffset=0;
    } else {
        q=qnode;
        offset=nweights;
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

TIOGA_GPU_GLOBAL
void g_update_sol_cart(
  int ncell_nf,
  int num_updates,
  int* q_ind_full,
  int* q_ind_cell_nd,
  double* q_val,
  double *qcell,
  double *qnode)
{
#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < num_updates)
#else
  for(int idx=0; idx<num_updates; idx++)
#endif
  {
    if(q_ind_full[idx] >= ncell_nf) {
      qnode[q_ind_cell_nd[idx]] = q_val[idx];
    }
    else {
      qcell[q_ind_cell_nd[idx]] = q_val[idx];
    }
  }
}
