#include "tioga_gpu.h"
#include "gpu_device_functions.h"
//void g_reset_iblanks(int *iblank, int nnodes)
TIOGA_GPU_GLOBAL
void g_reset_iblanks(TIOGA::MeshBlockInfo * m_info)
{
  int *iblank=m_info->iblank_node.dptr;
  int nnodes=m_info->num_nodes;
#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nnodes) {
#else
  for(int idx=0;idx<nnodes;idx++) {
#endif
    iblank[idx]=1;
  }

}

  //void g_adt_search(double *x, int **vconn,int *nc, int *nv, int ntypes,
  //           double *coord, double *adtExtents, int *adtIntegers, double *adtReals,
  //           int *elementList, int *donorId, double *xsearch, 
  //	     int ndim, int nelem, int nsearch)
TIOGA_GPU_GLOBAL
void g_adt_search(TIOGA::MeshBlockInfo* m_info, 
		    double *coord, double *adtExtents, int *adtIntegers, double *adtReals,
		    int *elementList, int *donorId, double *xsearch,
		    int ndim, int nelem, int nsearch, int myid)
{
  double *x = m_info->xyz.dptr;
  int ntypes=m_info->num_vert_per_elem.sz;
  int* nv= m_info->num_vert_per_elem.dptr;
  int* nc= m_info->num_cells_per_elem.dptr;
  int **vconn;
  int *vconn_ptrs[4];
  int ncells=0;
  int maxchecks=0;
  int imax;
  /*
  FILE *fp=fopen("log.dat","w");
  for(int i=0;i<ntypes;i++) ncells+=nc[i];
  fprintf(fp,"%d %d %d %d\n",nc[0],nc[1],nc[2],nc[3]);
  fprintf(fp,"ncells=%d\n",ncells);
  fprintf(fp,"ntypes=%d\n",ntypes);
  fprintf(fp,"nsearch=%d\n",nsearch);
  fprintf(fp,"%d %d %d %d\n",nv[0],nv[1],nv[2],nv[3]);
  fclose(fp);
  */
  for (int i=0; i < TIOGA::MeshBlockInfo::max_vertex_types; i++) {
    vconn_ptrs[i] = m_info->vertex_conn[i].dptr;
  }
  vconn=&vconn_ptrs[0];

#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < nsearch) {
#else
  for (int idx =0;idx<nsearch;idx++) {
#endif
    int cellIndex[2];
    int nchecks=0;
    donorId[idx]=-1;
    //printf("%d %d\n",idx,nsearch);
    d_searchIntersections_containment(idx,
                                     cellIndex,
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
    if (maxchecks < nchecks) {
     imax=idx;
     maxchecks=nchecks;
    }
    if (cellIndex[0] > -1 && cellIndex[1]==0) donorId[idx]=cellIndex[0];
  }
 //printf("maxchecks= %d %d %d %f %f %f\n",myid,imax,maxchecks,
 //                   xsearch[3*imax],xsearch[3*imax+1],xsearch[3*imax+2]);
}

TIOGA_GPU_GLOBAL
void g_interp_data(int *interpList_wcft,
		   double *interpList_weights,
		   int *interpList_inode,
		   int ninterp,
		   int nvar,
		   double *realData,
		   TIOGA::MeshBlockInfo* m_info)
{
  double *qnode=m_info->qnode.dptr;
  int nsend=nvar*ninterp;
#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nsend) 
#else
  for(int idx=0;idx<nsend;idx++) 
#endif
  {
    int i=idx/nvar;
    int k=idx-i*nvar;
    realData[idx]=0;
    for(int m=interpList_wcft[i];m<interpList_wcft[i+1];m++)
      {
	int inode=interpList_inode[m];
	double weight=interpList_weights[m];
        double val=qnode[inode*nvar+k]*weight;
	realData[idx]+=val;
      }
   }
}

TIOGA_GPU_GLOBAL
void g_update_sol(
  int num_updates,
  int* q_ind,
  double* q_val,
  TIOGA::MeshBlockInfo* m_info)
{
  double* qnode = m_info->qnode.dptr;

#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < num_updates)
#else
  for(int idx=0; idx<num_updates; idx++)
#endif
  {
    qnode[q_ind[idx]] = q_val[idx];
  }
}
