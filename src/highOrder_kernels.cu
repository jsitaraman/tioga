#include "cuda_runtime.h"
#include "error.hpp"

__global__
void interp_u(double* U_spts, double *U_out, int *donors, double *weights,
    int nFringe, int nSpts, int nVars, int estride, int sstride, int vstride)
{
  int fpt = (blockDim.x * blockIdx.x + threadIdx.x) / nVars;
  int var = (blockDim.x * blockIdx.x + threadIdx.x) % nVars;

  if (fpt >= nFringe || var >= nVars)
    return;

  int ind = nVars * fpt + var;
  int u_ind = donors[fpt] * estride + var * vstride;
  int w_ind = nSpts * fpt;

  U_out[ind] = 0.;

  for (int spt = 0; spt < nSpts; spt++)
  {
    U_out[ind] += weights[w_ind + spt] * U_spts[u_ind + spt * sstride];
  }
}

void interp_u_wrapper(double *U_spts, double *U_out, int *donors,
    double *weights, int nFringe, int nSpts, int nVars, int estride,
    int sstride, int vstride)
{
  int threads = 128;
  int blocks = (nVars * nFringe + threads - 1) / threads;

  interp_u<<<blocks, threads>>>(U_spts, U_out, donors, weights, nFringe, nSpts,
      nVars, estride, sstride, vstride);
}
