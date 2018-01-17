#include "cuda_runtime.h"
#include "error.hpp"

template<unsigned int nVars>
__global__
void interp_u(const double* __restrict__ U_spts, double *U_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const int* __restrict__ out_inds, int nFringe, int nSpts,
    int estride, int sstride, int vstride)
{
  const int fpt = blockDim.x * blockIdx.x + threadIdx.x;

  if (fpt >= nFringe)
    return;

  int ind = nVars * out_inds[fpt];
  int u_ind = donors[fpt] * estride;
  int w_ind = nSpts * fpt;

  double sum[nVars] = {0.0};

  for (int spt = 0; spt < nSpts; spt++)
  {
    double wt = weights[w_ind+spt];
    for (int var = 0; var < nVars; var++)
      sum[var] += wt * U_spts[u_ind + spt*sstride + var*vstride];
  }

  for (int var = 0; var < nVars; var++)
    U_out[ind+var] = sum[var];
}

void interp_u_wrapper(double *U_spts, double *U_out, int *donors,
    double *weights, int* out_inds, int nFringe, int nSpts, int nVars, int estride,
    int sstride, int vstride, cudaStream_t stream_h)
{
  unsigned int threads = 128;
  unsigned int blocks = (nFringe + threads - 1) / threads;

  if (nVars == 1)
    interp_u<1><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, out_inds,
        nFringe, nSpts, estride, sstride, vstride);
  else if (nVars == 4)
    interp_u<4><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, out_inds,
        nFringe, nSpts, estride, sstride, vstride);
  else if (nVars == 5)
    interp_u<5><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, out_inds,
        nFringe, nSpts, estride, sstride, vstride);

  check_error();
}

template<unsigned int nVars>
__global__
void interp_u_types(const double* const*  U_spts, double *U_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const char* __restrict__ etypes, const int* __restrict__ wgt_inds,
    const int* __restrict__ out_inds, int nFringe,
    const int* __restrict__ nweights, const int* __restrict__ strides)
{
  const int fpt = blockDim.x * blockIdx.x + threadIdx.x;

  if (fpt >= nFringe)
    return;

  const int ind = nVars * out_inds[fpt];
  const int N = (int)etypes[fpt];
  const int u_ind = donors[fpt] * strides[3*N+0];
  const int w_ind = wgt_inds[fpt];
  const int nSpts = nweights[fpt];
  const int sstride = strides[3*N+1];
  const int vstride = strides[3*N+2];

  double sum[nVars] = {0.0};

  for (int spt = 0; spt < nSpts; spt++)
  {
    double wt = weights[w_ind+spt];
    for (int var = 0; var < nVars; var++)
      sum[var] += wt * U_spts[N][u_ind + spt*sstride + var*vstride];
  }

  for (int var = 0; var < nVars; var++)
    U_out[ind+var] = sum[var];
}

void interp_u_types_wrapper(double **U_spts, double *U_out, int *donors,
    double *weights, char *etypes, int* wgt_inds, int* out_inds, int nFringe,
    int* nSpts, int nVars, int *strides, cudaStream_t stream_h)
{
  unsigned int threads = 128;
  unsigned int blocks = (nFringe + threads - 1) / threads;

  if (nVars == 1)
    interp_u_types<1><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, etypes,
        wgt_inds, out_inds, nFringe, nSpts, strides);
  else if (nVars == 4)
    interp_u_types<4><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, etypes,
        wgt_inds, out_inds, nFringe, nSpts, strides);
  else if (nVars == 5)
    interp_u_types<5><<<blocks, threads, 0, stream_h>>>(U_spts, U_out, donors, weights, etypes,
        wgt_inds, out_inds, nFringe, nSpts, strides);

  check_error();
}


template <unsigned int nDims, unsigned int nVars>
__global__
void interp_du_types(const double* const*  dU_spts, double *dU_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const char* __restrict__ etypes, const int* __restrict__ wgt_inds,
    const int* __restrict__ out_inds, int nFringe,
    const int* __restrict__ nweights, const int* __restrict__ strides)
{
  const int fpt = blockDim.x * blockIdx.x + threadIdx.x;

  if (fpt >= nFringe)
    return;

  const int N = (int)etypes[fpt];
  const int u_ind = donors[fpt] * strides[4*N+0];
  const int w_ind = wgt_inds[fpt];
  const int nSpts = nweights[fpt];
  const int sstride = strides[4*N+1];
  const int vstride = strides[4*N+2];
  const int dstride = strides[4*N+3];

  double sum[nDims][nVars] = {0.0};

  for (int spt = 0; spt < nSpts; spt++)
  {
    double wgt = weights[w_ind + spt];
    for (int dim = 0; dim < nDims; dim++)
      for (int var = 0; var < nVars; var++)
        sum[dim][var] += wgt * dU_spts[N][u_ind + spt*sstride + dim*dstride + var * vstride];
  }

  for (int dim = 0; dim < nDims; dim++)
  {
    for (int var = 0; var < nVars; var++)
    {
      int ind = var + nVars * (dim + nDims * out_inds[fpt]);
      dU_out[ind] = sum[dim][var];
    }
  }
}

void interp_du_types_wrapper(double **dU_spts, double *dU_out, int *donors,
    double *weights, char *etypes, int* wgt_inds, int* out_inds, int nFringe,
    int* nSpts, int nVars, int nDims, int *strides, cudaStream_t stream_h)
{
  int threads = 128;
  int blocks = (nFringe + threads - 1) / threads;

  if (nDims == 3)
  {
    if (nVars == 1)
      interp_du_types<3,1><<<blocks, threads, 0, stream_h>>>(dU_spts, dU_out, donors,
          weights, etypes, wgt_inds, out_inds, nFringe, nSpts, strides);
    else if (nVars == 5)
      interp_du_types<3,5><<<blocks, threads, 0, stream_h>>>(dU_spts, dU_out, donors,
          weights, etypes, wgt_inds, out_inds, nFringe, nSpts, strides);
    else
      FatalError("3D nVars case not recognized (expecting 1 or 5)");
  }
  else
    FatalError("TIOGA support for 3D only currently!");

  check_error();
}

template <unsigned int nDims, unsigned int nVars>
__global__
void interp_du(const double* __restrict__ dU_spts, double *dU_out,
    const int* __restrict__ donors, const double* __restrict__ weights,
    const int* __restrict__ out_inds, int nFringe, int nSpts,
    int estride, int sstride, int vstride, int dstride)
{
  const int fpt = blockDim.x * blockIdx.x + threadIdx.x;

  if (fpt >= nFringe)
    return;

  int u_ind = donors[fpt] * estride;
  int w_ind = nSpts * fpt;

  double sum[nDims][nVars] = {0.0};

  for (int spt = 0; spt < nSpts; spt++)
  {
    double wgt = weights[w_ind + spt];
    for (int dim = 0; dim < nDims; dim++)
      for (int var = 0; var < nVars; var++)
        sum[dim][var] += wgt * dU_spts[u_ind + spt*sstride + dim*dstride + var * vstride];
  }

  for (int dim = 0; dim < nDims; dim++)
  {
    for (int var = 0; var < nVars; var++)
    {
      int ind = var + nVars * (dim + nDims * out_inds[fpt]);
      dU_out[ind] = sum[dim][var];
    }
  }
}

void interp_du_wrapper(double *dU_spts, double *dU_out, int *donors,
    double *weights, int* out_inds, int nFringe, int nSpts, int nVars,
    int nDims, int estride, int sstride, int vstride, int dstride, cudaStream_t stream_h)
{
  int threads = 128;
  int blocks = (nFringe + threads - 1) / threads;

  if (nDims == 3)
  {
    if (nVars == 1)
      interp_du<3,1><<<blocks, threads, 0, stream_h>>>(dU_spts, dU_out, donors, weights,
        out_inds, nFringe, nSpts, estride, sstride, vstride, dstride);
    else if (nVars == 5)
      interp_du<3,5><<<blocks, threads, 0, stream_h>>>(dU_spts, dU_out, donors, weights,
        out_inds, nFringe, nSpts, estride, sstride, vstride, dstride);
    else
      FatalError("3D nVars case not recognized (expecting 1 or 5)");
  }
  else
    FatalError("TIOGA support for 3D only currently!");

  check_error();
}
