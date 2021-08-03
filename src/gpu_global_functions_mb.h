#ifndef GPU_GLOBAL_FUNCTIONS_MB_H
#define GPU_GLOBAL_FUNCTIONS_MB_H

#include "TiogaMeshInfo.h"
#include "tioga_gpu.h"

TIOGA_GPU_GLOBAL
void g_reset_iblanks(TIOGA::MeshBlockInfo * m_info);

TIOGA_GPU_GLOBAL
void g_adt_search(TIOGA::MeshBlockInfo* m_info,
  double *coord, double *adtExtents, int *adtIntegers, double *adtReals,
  int *elementList, int *donorId, double *xsearch,
  int ndim, int nelem, int nsearch, int myid);

TIOGA_GPU_GLOBAL
void g_interp_data(int *interpList_wcft,
  double *interpList_weights,
  int *interpList_inode,
  int ninterp,
  int nvar,
  double *realData,
  TIOGA::MeshBlockInfo* m_info);

TIOGA_GPU_GLOBAL
void g_update_sol(
  int num_updates,
  int* q_ind,
  double* q_val,
  TIOGA::MeshBlockInfo* m_info);

#endif /* GPU_GLOBAL_FUNCTIONS_MB_H */
