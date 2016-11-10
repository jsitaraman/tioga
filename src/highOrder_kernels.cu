#include "cuda_runtime.h"
#include "error.hpp"

__global__
void interp_u(const double* __restrict__ U_spts, double *U_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const int* __restrict__ out_inds, int nFringe, int nSpts,
    int nVars, int estride, int sstride, int vstride)
{
  const int fpt = (blockDim.x * blockIdx.x + threadIdx.x) / nVars;
  const int var = (blockDim.x * blockIdx.x + threadIdx.x) % nVars;

  if (fpt >= nFringe)
    return;

  int ind = nVars * out_inds[fpt] + var;
  int u_ind = donors[fpt] * estride + var * vstride;
  int w_ind = nSpts * fpt;

  double sum = 0;

  for (int spt = 0; spt < nSpts; spt++)
    sum += weights[w_ind+spt] * U_spts[u_ind + spt*sstride];

  U_out[ind] = sum;
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
void interp_du(const double* __restrict__ dU_spts, double *dU_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const int* __restrict__ out_inds, int nFringe, int nSpts,
    int nVars, int estride, int sstride, int vstride, int dstride)
{
  const int fpt = (blockDim.x * blockIdx.x + threadIdx.x) / nVars;
  const int var = (blockDim.x * blockIdx.x + threadIdx.x) % nVars;

  if (fpt >= nFringe)
    return;

  int u_ind = donors[fpt] * estride + var * vstride;
  int w_ind = nSpts * fpt;

  double sum;

  for (int dim = 0; dim < nDims; dim++)
  {
    sum = 0.;
    for (int spt = 0; spt < nSpts; spt++)
      sum += weights[w_ind+spt] * dU_spts[u_ind + spt*sstride + dim*dstride];

    int ind = nVars * (dim + nDims * out_inds[fpt]) + var;
    dU_out[ind] = sum;
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
