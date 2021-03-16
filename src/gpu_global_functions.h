#include "tioga_gpu.h"
#include "gpu_device_functions.h"
#include<stdio.h>
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
void g_adt_search(double *x, TIOGA::MeshBlockInfo* m_info,int ntypes, int *nc, int *nv, int *vconn, 
		    double *coord, double *adtExtents, int *adtIntegers, double *adtReals,
		    int *elementList, int *donorId, double *xsearch,
		    int ndim, int nelem, int nsearch)
{

  int numverts[4][6]={3,3,3,3,0,0,4,3,3,3,3,0,3,4,4,4,3,0,4,4,4,4,4,4};
  int faceInfo[4][24]={1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0,
		       1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0,
		       1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0,
		       1,2,3,4,1,5,6,2,2,6,7,3,3,7,8,4,1,4,8,5,5,8,7,6};
  int nchecks;

  //double *x = m_info->xyz.dptr;
  /*
  int ntypes=m_info->num_vert_per_elem.sz;
  int* nv= m_info->num_vert_per_elem.dptr;
  int* nc= m_info->num_cells_per_elem.dptr;
  int **vconn;
  int *vconn_ptrs[4];
  int ncells=0;
  for(int i=0;i<ntypes;i++) ncells+=nc[i];

  for (int i=0; i < TIOGA::MeshBlockInfo::max_vertex_types; i++) {
    vconn_ptrs[i] = m_info->vertex_conn[i].dptr;
  }
  vconn=&vconn_ptrs[0];
  */

#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < nsearch) {
#else
  for (int idx =0;idx<nsearch;idx++) {
#endif
    //    int cellIndex[2];
    donorId[idx]=-1;
    d_searchIntersections_containment(&(donorId[idx]),
                                     adtIntegers,
                                     adtReals,
                                     coord,
                                     &(xsearch[3*idx]),
                                     x, nc, nv, ntypes,
                                     elementList,
                                     vconn,
                                     nelem,
                                     ndim,
				      &nchecks,
				      numverts,faceInfo);
    //    if (cellIndex[0] > -1 && cellIndex[1]==0) donorId[idx]=cellIndex[0];
  }
}
