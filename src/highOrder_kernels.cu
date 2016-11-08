#include "cuda_runtime.h"
#include "error.hpp"

__global__
void interp_u(double* U_spts, double *U_out, int *donors, double *weights,
    int *out_inds, int nFringe, int nSpts, int nVars, int estride, int sstride,
    int vstride)
{
  const int fpt = (blockDim.x * blockIdx.x + threadIdx.x) / nVars;
  const int var = (blockDim.x * blockIdx.x + threadIdx.x) % nVars;

  if (fpt >= nFringe || var >= nVars)
    return;

  int ind = nVars * out_inds[fpt] + var;
  int u_ind = donors[fpt] * estride + var * vstride;
  int w_ind = nSpts * fpt;

  U_out[ind] = 0.;

  for (int spt = 0; spt < nSpts; spt++)
  {
    U_out[ind] += weights[w_ind+spt] * U_spts[u_ind + spt*sstride];
  }
}

void interp_u_wrapper(double *U_spts, double *U_out, int *donors,
    double *weights, int* out_inds, int nFringe, int nSpts, int nVars, int estride,
    int sstride, int vstride)
{
  unsigned int threads = 128;
  unsigned int blocks = (nVars * nFringe + threads - 1) / threads;

  interp_u<<<blocks, threads>>>(U_spts, U_out, donors, weights, out_inds,
      nFringe, nSpts, nVars, estride, sstride, vstride);

  check_error();
}


template <unsigned int nDims>
__global__
void interp_du(double* dU_spts, double *dU_out, int *donors, double *weights,
    int *out_inds, int nFringe, int nSpts, int nVars, int estride, int sstride,
    int vstride, int dstride)
{
  const int fpt = (blockDim.x * blockIdx.x + threadIdx.x) / nVars;
  const int var = (blockDim.x * blockIdx.x + threadIdx.x) % nVars;

  if (fpt >= nFringe || var >= nVars)
    return;

  int u_ind = donors[fpt] * estride + var * vstride;
  int w_ind = nSpts * fpt;

  for (int dim = 0; dim < nDims; dim++)
  {
    int ind = nVars * (dim + nDims * out_inds[fpt]) + var;
    dU_out[ind] = 0.;
    for (int spt = 0; spt < nSpts; spt++)
    {
      dU_out[ind] += weights[w_ind+spt] * dU_spts[u_ind + spt*sstride + dim*dstride];
    }
  }
}

void interp_du_wrapper(double *dU_spts, double *dU_out, int *donors,
    double *weights, int* out_inds, int nFringe, int nSpts, int nVars,
    int nDims, int estride, int sstride, int vstride, int dstride)
{
  unsigned int threads = 128;
  unsigned int blocks = (nVars * nFringe + threads - 1) / threads;

  if (nDims == 3)
    interp_du<3><<<blocks, threads>>>(dU_spts, dU_out, donors, weights,
        out_inds, nFringe, nSpts, nVars, estride, sstride, vstride, dstride);
  else
    FatalError("TIOGA support for 3D only currently!");

  check_error();
}
