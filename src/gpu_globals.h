#include "tioga_gpu.h"

TIOGA_GPU_GLOBAL
void g_reset_iblanks(int *iblank, int nnodes)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nnodes) {
    iblank[idx]=1;
  }

}
