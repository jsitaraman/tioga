#include "error.hpp"

void interp_u_wrapper(double *U_spts, double *U_out, int *donors,
    double *weights, int *out_inds, int nFringe, int nSpts, int nVars,
    int estride, int sstride, int vstride, cudaStream_t stream_h);


void interp_du_wrapper(double *dU_spts, double *dU_out, int *donors,
    double *weights, int *out_inds, int nFringe, int nSpts, int nVars,
    int nDims, int estride, int sstride, int vstride, int dstride,
    cudaStream_t stream_h);
