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
		    int ndim, int nelem, int nsearch)
{
  double *x = m_info->xyz.dptr;
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

#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
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
