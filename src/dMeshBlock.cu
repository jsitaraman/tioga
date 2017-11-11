#include "dMeshBlock.h"
#include "funcs.hpp"

#include "device_functions.h"
#include "math.h"

#define MAX_UCHAR 255

/* --- Handy Vector Operation Macros --- */

#define NF1 16 // 16-32 depending on unstructured-ness of grid & desire for robustness
#define NF2  4 // 3-6 depending on unstructured-ness of grid & desire for robustness

#define CROSS(a, b, c) { \
  c[0] = a[1]*b[2] - a[2]*b[1]; \
  c[1] = a[2]*b[0] - a[0]*b[2]; \
  c[2] = a[0]*b[1] - a[1]*b[0]; }

#define CROSS4(a1, a2, b1, b2, c) { \
  c[0] = (a1[1]-a2[1])*(b1[2]-b2[2]) - (a1[2]-a2[2])*(b1[1]-b2[1]); \
  c[1] = (a1[2]-a2[2])*(b1[0]-b2[0]) - (a1[0]-a2[0])*(b1[2]-b2[2]); \
  c[2] = (a1[0]-a2[0])*(b1[1]-b2[1]) - (a1[1]-a2[1])*(b1[0]-b2[0]); }

#define DOT(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])

static
__device__ __forceinline__
float DOTCROSS4(const float* __restrict__ c,
                 const float* __restrict__ a1, const float* __restrict__ a2,
                 const float* __restrict__ b1, const float* __restrict__ b2)
{
  float d[3];
  CROSS4(a1,a2,b1,b2,d)
  return DOT(c,d);
}


/* --- Misc. Helpful CUDA kernels --- */

__device__ __forceinline__
void print_nodes(const double* pts, int id, int npts)
{
  int idx = threadIdx.x;
  for (int tid = 0; tid < 32; tid++)
  {
    if (idx == tid)
    {
      printf("Points%d = [",id);
      for (int i = 0; i < npts - 1; i++)
        printf("%f %f %f;\n",pts[3*i+0],pts[3*i+1],pts[3*i+2]);

      int I = npts-1;
      printf("%f %f %f];\n",pts[3*I+0],pts[3*I+1],pts[3*I+2]);
    }
  }
}

__device__ __forceinline__
void print_nodes(const float* pts, int id, int npts)
{
  int idx = threadIdx.x;
  for (int tid = 0; tid < 32; tid++)
  {
    if (idx == tid)
    {
      printf("Points%d = [",id);
      for (int i = 0; i < npts - 1; i++)
        printf("%f %f %f;\n",pts[3*i+0],pts[3*i+1],pts[3*i+2]);

      int I = npts-1;
      printf("%f %f %f];\n",pts[3*I+0],pts[3*I+1],pts[3*I+2]);
    }
  }
}

#define WARP_SZ 32

__device__
inline int lane_id(void) { return threadIdx.x % WARP_SZ; }

__device__
inline int warp_bcast(int v, int leader) { return __shfl(v, leader); }

__device__ __forceinline__
float warpAllReduceMin(float val)
{
  for (int mask = warpSize/2; mask > 0; mask /= 2)
    val = fminf(val, __shfl_xor(val, mask));
  return val;
}

__device__
int floatToOrderedInt(float floatVal)
{
  int intVal = __float_as_int(floatVal);

  return (intVal >= 0) ? intVal : intVal ^ 0x8FFFFFFF;
}

__device__
unsigned int floatToUint(float fval)
{
  unsigned int ival = __float_as_uint(fval);
  unsigned int mask = -int(ival >> 31) | 0x80000000;
  return ival ^ mask;
}

__device__
float uintToFloat(unsigned int ival)
{
  unsigned int mask = ((ival >> 31) - 1) | 0x80000000;
  return __uint_as_float(ival ^ mask);
}

__device__
float orderedIntToFloat(int intVal)
{
  return __int_as_float( (intVal >= 0) ? intVal : intVal ^ 0x8FFFFFFF );
}

__device__ float atomicMaxf(float* address, float val)
{
  //int *iaddr = (int*)address;
  int old = __float_as_int(*address);
  int assumed;
  while (val > __int_as_float(old))
  {
    assumed = old;
    old = atomicCAS((int*)address, assumed, __float_as_int(val));
  }
  return __int_as_float(old);
}

__device__ float atomicMinf(float* address, float val)
{
  int old = __float_as_int(*address);
  int assumed;
  while (val < __int_as_float(old))
  {
    assumed = old;
    old = atomicCAS((int*)address, assumed, __float_as_int(val));
  }
  return __int_as_float(old);
}

/*! Warp-aggregated atomic increment
 *  https://devblogs.nvidia.com/parallelforall/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/ */
__device__
int atomicAggInc(int *ctr)
{
  int mask = __ballot(1);
  // select the leader
  int leader = __ffs(mask) - 1;
  // leader does the update
  int res;
  if (lane_id() == leader)
    res = atomicAdd(ctr, __popc(mask));
  // brodcast result
  res = warp_bcast(res, leader);
  // each thread computes its own value
  return res + __popc(mask & ((1 << lane_id()) - 1));
}

/* ------ dMeshBlock Member Functions ------ */

void dMeshBlock::dataToDevice(int ndims, int nnodes, int ncells, int ncells_adt,
    int nsearch, int* nv, int* nc, int* eleList, double* eleBBox, int* isearch,
    double* xsearch, int rank)
{
  this->nnodes = nnodes;
  this->ncells = ncells;
  this->nc_adt = ncells_adt;

  this->nv = nv;
  this->nc = nc;

  nvert = nv[0];

  this->rank = rank;

  this->eleBBox.assign(eleBBox, ncells_adt*ndims*2);
  this->eleList.assign(eleList, ncells_adt);

  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);

  auto ijk2gmsh_h = tg_funcs::structured_to_gmsh_hex(nvert);
  ijk2gmsh.assign(ijk2gmsh_h.data(), ijk2gmsh_h.size());

  int nSide = std::cbrt(nvert);
  std::vector<double> xlist_h(nSide);
  std::vector<float> xlistf_h(nSide);
  double dxi = 2./(nSide-1);

  for (int i = 0; i < nSide; i++)
    xlist_h[i] = -1. + i*dxi;

  for (int i = 0; i < nSide; i++)
    xlistf_h[i] = xlist_h[i];

  xlist.assign(xlist_h.data(), xlist_h.size());
  xlistf.assign(xlistf_h.data(), xlistf_h.size());
}

void dMeshBlock::extraDataToDevice(int* vconn)
{
//  c2v.assign(vconn, nvert*ncells);
}

void dMeshBlock::assignHoleMap(bool hasWall, int* nx, int* sam, double* extents)
{
  if (hasWall)
  {
    int size = nx[0]*nx[1]*nx[2];

    std::vector<char> tmp_sam(size);
    for (int i = 0; i < size; i++)
      tmp_sam[i] = (char)sam[i];

    double dx[3];
    for (int d = 0; d < 3; d++)
      dx[d] = (extents[d+3] - extents[d]) / nx[d];

    hm_sam.assign(tmp_sam.data(), size);
    hm_extents.assign(extents, 6);
    hm_nx.assign(nx, 3);
    hm_dx.assign(dx, 3);
  }
  else
  {
    clearHoleMap();
  }
}

void dMeshBlock::clearHoleMap(void)
{
  int nx[3] = {0,0,0};
  double dx[3] = {0,0,0};
  double extents[6] = {0,0,0,0,0,0};

  hm_sam.resize(0);

  hm_nx.assign(nx, 3);
  hm_dx.assign(dx, 3);
  hm_extents.assign(extents, 6);
}

void dMeshBlock::updateADTData(int ncells_adt, int* eleList, double* eleBBox)
{
  this->eleBBox.assign(eleBBox, ncells_adt*nDims*2);
  this->eleList.assign(eleList, ncells_adt);
}

void dMeshBlock::updateSearchPoints(int nsearch, int *isearch, double *xsearch)
{
  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);
}

void dMeshBlock::setDeviceData(double* vx, double* ex, int* ibc, int* ibf)
{
  x = vx;
  iblank_cell = ibc;
  iblank_face = ibf;
  coord = ex;
}

void dMeshBlock::setTransform(double* mat, double* off, int ndim)
{
  if (ndim != nDims)
    ThrowException("dMeshBlock::set_transform: input ndim != nDims");

  rrot = true; /// WORKING ON ADT REBUILD - DISABLED RROT
  Rmat.assign(mat, ndim*ndim);
  offset.assign(off, ndim);
}

/* ---------------------------- Geometric Functions --------------------------- */

__device__ __forceinline__
int oobbCheck(float* vA, float* xA, float* dxA, float* vB, float* xB, float* dxB)
{
  double eps = 1e-10;

  // D = distance between centers
  float D[3];
  for (int i = 0; i < 3; i++)
    D[i] = xB[i] - xA[i];

  // C = scalar product of axes
  float c[3][3];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      c[i][j] = 0;
      for (int k = 0; k < 3; k++)
        c[i][j] = c[i][j] + vA[3*i+k]*vB[3*j+k];
    }
  }

  // separating axes based on the faces of box A
  for (int i = 0; i < 3; i++)
  {
    float r0 = dxA[i];
    float r1 = 0;
    float r = 0;
    for (int j = 0; j < 3; j++)
    {
      r1 += dxB[j] * fabs(c[i][j]);
      r += fabs(vA[3*i+j]) * D[j];
    }

    if (r > (r0+r1+eps))
      return 0;
  }

  // separating axes based on the faces of box B
  for (int i = 0; i < 3; i++)
  {
    float r1 = dxB[i];
    float r0 = 0;
    float r = 0;
    for (int j = 0; j < 3; j++)
    {
      r0 += dxA[j] * fabs(c[j][i]);
      r += fabs(vB[3*i+j]) * D[j];
    }
    if (r > (r0+r1+eps)) return 0;
  }

  // cross products
  for (int i = 0; i < 3; i++)
  {
    int i1 = (i+1)%3;
    int i2 = (i+2)%3;
    for (int j = 0; j < 3; j++)
    {
      int j1 = (j+1)%3;
      int j2 = (j+2)%3;

      float r0 = dxA[i1] * fabs(c[i2][j]) + dxA[i2] * fabs(c[i1][j]);
      float r1 = dxB[j1] * fabs(c[i][j2]) + dxB[j2] * fabs(c[i][j1]);

      float d2 = 0;
      float d1 = 0;
      for (int k = 0; k < 3; k++)
      {
        d2 += vA[3*i2+k]*D[k];
        d1 += vA[3*i1+k]*D[k];
      }

      float r = fabs( c[i1][j]*d2 - c[i2][j]*d1 );

      if (r > (r0+r1+eps)) {
        return 0;
      }
    }
  }

  // return 1 if no separation can be found [intersection]
  return 1;
}

/* From Demmel 1997, Applied Numerical Linear Algebra, pp. 232-235 */

//! Specific to a 3x3 matrix [used to get oriented bounding box axes]
__device__ __forceinline__
void jacobi_rotation(float A[3][3], float* __restrict__ J, int j, int k)
{
  if (abs(A[j][k]) > 1e-10)
  {
    float tau = (A[j][j] - A[k][k]) / (2*A[j][k]);
    float t = sign_d(tau) / (abs(tau) + sqrt(1+tau*tau));
    float c = 1. / (sqrt(1+t*t));
    float s = c*t;

    // Temporary matrix for matrix-matrix multiplication
    float B[6];

    // Apply R^T from left - rotates rows j,k by theta
    for (int d = 0; d < 3; d++)
    {
      B[d]   =  c*A[d][j] + s*A[d][k];
      B[d+3] = -s*A[d][j] + c*A[d][k];
    }

    // Update A to be result of R^T*A (Row i not changed)
    for (int d = 0; d < 3; d++)
    {
      A[j][d] = B[d];
      A[k][d] = B[d+3];
    }

    // Copy columns j,k of A to B
    for (int d = 0; d < 3; d++)
    {
      B[d]   = A[d][j];
      B[d+3] = A[d][k];
    }

    // Apply R from right - rotates cols j,k by theta
    for (int d = 0; d < 3; d++)
    {
      A[d][j] =  c*B[d] + s*B[d+3];
      A[d][k] = -s*B[d] + c*B[d+3];
    }

    // Copy columns j,k of J into tmp storage for multiplication
    for (int d = 0; d < 3; d++)
    {
      B[d]   = J[3*d+j];
      B[d+3] = J[3*d+k];
    }

    // Update J [eventual eigenvectors] [J = J*R]
    for (int d = 0; d < 3; d++)
    {
      J[3*d+j] =  c*B[d] + s*B[d+3];
      J[3*d+k] = -s*B[d] + c*B[d+3];
    }
  }
}

__device__ __forceinline__
float off_norm(float A[3][3])
{
  return sqrt(A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2]);
}

//! Get the eigenvectors of a 3x3 symmetric matrix [covariance matrix]
__device__ __forceinline__
void get_eigenvectors(float mat[3][3], float* __restrict__ evecs)
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < i; j++)
      evecs[3*i+j] = 0;

    evecs[3*i+i] = 1;

    for (int j = i+1; j < 3; j++)
      evecs[3*i+j] = 0;
  }

  //while (off_norm(mat) < 1.e-8)
  for (int i = 0; i < 4; i++)
  {
    jacobi_rotation(mat, evecs, 0, 1);
    jacobi_rotation(mat, evecs, 0, 2);
    jacobi_rotation(mat, evecs, 1, 2);
  }
}

/*! Get the 3x3 covariance matrix for a set of points
 *  WARNING: upon return 'pts' will be shifted to place centroid at origin */
__device__ __forceinline__
void get_covariance(float* __restrict__ pts, int npts, float xc[3], float cov[3][3])
{
  // Get the centroid of all points
  xc[0] = xc[1] = xc[2] = 0.f;

  for (int i = 0; i < npts; i++)
    for (int d = 0; d < 3; d++)
      xc[d] += pts[3*i+d];

  for (int d = 0; d < 3; d++)
    xc[d] /= (float)npts;

  // Shift the centroid to the origin
  for (int i = 0; i < npts; i++)
    for (int d = 0; d < 3; d++)
      pts[3*i+d] -= xc[d];

  // Compute the covarience [just the upper-triangular terms]
  for (int d1 = 0; d1 < 3; d1++)
  {
    for (int d2 = d1; d2 < 3; d2++)
    {
      cov[d1][d2] = 0.f;
      for (int i = 0; i < npts; i++)
        cov[d1][d2] += pts[3*i+d1]*pts[3*i+d2];
      cov[d1][d2] /= (float)npts;
    }
  }

  // Apply symmetry
  for (int d1 = 0; d1 < 3; d1++)
    for (int d2 = d1; d2 < 3; d2++)
      cov[d2][d1] = cov[d1][d2];
}

/*! Get the oriented bounding box of a set of points using the eigenvectors
 *  of the covariance matrix as the axis of the box
 *  First 9 values of obb are the new axes; next 6 are transformed bbox */
__device__ __forceinline__
void getOBB(float* __restrict__ pts, int npts, float* obb)
{
  float cov[3][3], xc[3], axes[9];
  get_covariance(pts, npts, xc, cov);

  get_eigenvectors(cov, axes);

  for (int i = 0; i < 9; i++)
    obb[i] = axes[i];

  for (int i = 0; i < 3; i++)
  {
    obb[9+i]   =  BIG_DOUBLE;
    obb[9+i+3] = -BIG_DOUBLE;
  }

  for (int i = 0; i < npts; i++)
  {
    float pt[3] = {0.f, 0.f, 0.f};
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        pt[j] += obb[3*k+j] * (pts[3*i+k]+xc[k]);

    for (int d = 0; d < 3; d++)
    {
      obb[9+d]   = fminf(obb[9+d],   pt[d]);
      obb[9+d+3] = fmaxf(obb[9+d+3], pt[d]);
    }
  }
}

/* ---------------------------- Direct Cut Method Functions --------------------------- */

static
__device__
float lineSegmentDistance(float *p1, float *p2, float *p3, float *p4, float *dx)
{
  // Get the line equations
  const float U[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
  const float V[3] = {p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2]};
  const float W[3] = {p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2]};
  const float uu = U[0]*U[0] + U[1]*U[1] + U[2]*U[2];
  const float vv = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
  const float uv = U[0]*V[0] + U[1]*V[1] + U[2]*V[2];

  const float uw = U[0]*W[0] + U[1]*W[1] + U[2]*W[2];
  const float vw = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];

  const float den = uu*vv - uv*uv;

  // NOTE: not finding exact minimum distance between the line segments in all
  // cases; plenty close enough for our purposes
  // (see http://geomalgorithms.com/a07-_distance.html for full algo)

  // Calculate line parameters (if nearly parallel, set one & calculate other)
  float s = (den < 1e-7f) ? 0.0f : (uv*vw - vv*uw) / den;
  float t = (den < 1e-7f) ? uw / uv: (uu*vw - uv*uw) / den;

  s = fmin(fmax(s, 0.f), 1.f);
  t = fmin(fmax(t, 0.f), 1.f);

  // vec = closest distance from segment 1 to segment 2
  for (int i = 0; i < 3; i++)
    dx[i] = t*V[i] - s*U[i] - W[i];

  float dist = 0.f;
  for (int i = 0; i < 3; i++)
    dist += dx[i]*dx[i];

  return sqrt(dist);
}

/*! Modified Moller triangle-triangle intersection algorithm
 *  Determines if triangles intersect, or returns an approximate minimum
 *  distance between them otherwise
 *  Also returns vector of minimum distance from T1 to T2 */
static
__device__
float triTriDistanceVec(float* __restrict__ T1, float* __restrict__ T2,
    float* __restrict__ minVec, float tol)
{
  float dist = 1e15f;
  float vec[3];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      int i2 = (i+1) % 3;
      int j2 = (j+1) % 3;
      float D = lineSegmentDistance(&T1[3*i], &T1[3*i2], &T2[3*j], &T2[3*j2], vec);

      if (D < dist)
      {
        for (int d = 0; d < 3; d++)
          minVec[d] = vec[d];
        dist = D;
      }
    }
  }

  // Pointers to points
  const float* V01 = T1;
  const float* V11 = T1+3;
  const float* V21 = T1+6;

  const float* V02 = T2;
  const float* V12 = T2+3;
  const float* V22 = T2+6;

  float N1[3], N2[3];

  // Plane for Triangle 1
  CROSS4(V11,V01, V21,V01, N1);

  float norm = NORM(N1);

  // Plane for Triangle 2
  for (int d = 0; d < 3; d++)
    N1[d] /= norm;

  float d1 = -DOT(N1,V01);

  CROSS4(V12,V02, V22,V02, N2);

  norm = NORM(N2);

  for (int d = 0; d < 3; d++)
    N2[d] /= norm;

  float d2 = -DOT(N2,V02);

  // Signed distances of T1's vertices to T2's plane
  float d01 = DOT(N2,V01) + d2;
  float d11 = DOT(N2,V11) + d2;
  float d21 = DOT(N2,V21) + d2;

  float d02 = DOT(N1,V02) + d1;
  float d12 = DOT(N1,V12) + d1;
  float d22 = DOT(N1,V22) + d1;

  // Round values near 0 to 0
  d01 = (fabs(d01) < 1e-10) ? 0 : d01;
  d11 = (fabs(d11) < 1e-10) ? 0 : d11;
  d21 = (fabs(d21) < 1e-10) ? 0 : d21;

  d02 = (fabs(d02) < 1e-10) ? 0 : d02;
  d12 = (fabs(d12) < 1e-10) ? 0 : d12;
  d22 = (fabs(d22) < 1e-10) ? 0 : d22;

  if (fabs(d01) + fabs(d11) + fabs(d21) < 3*tol ||
      fabs(d02) + fabs(d12) + fabs(d22) < 3*tol)
  {
    // Approximately coplanar; check if one triangle is inside the other /

    // Check if a point in T1 is inside T2
    bool inside = true;
    inside = inside && DOTCROSS4(N2, V12,V02, V01,V02) > 0;
    inside = inside && DOTCROSS4(N2, V02,V22, V01,V22) > 0;
    inside = inside && DOTCROSS4(N2, V22,V12, V01,V12) > 0;

    if (inside) return 0.;

    // Check if a point in T2 is inside T1
    inside = true;
    inside = inside && DOTCROSS4(N1, V11,V01, V02,V01) > 0;
    inside = inside && DOTCROSS4(N1, V01,V21, V02,V21) > 0;
    inside = inside && DOTCROSS4(N1, V21,V11, V02,V11) > 0;

    if (inside) return 0.;
  }

  bool noTouch = false;

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d01) == sgn(d11) && sgn(d01) == sgn(d21)) // && fabs(d01) > tol)
  {
    noTouch = true;

    // No intersection; check if projection of points provides closer distance
    if (fabs(d01) < dist)
    {
      float P01[3];
      for (int d = 0; d < 3; d++)
        P01[d] = V01[d] - N2[d]*d01;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P01,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P01,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P01,V12) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = -N2[i]*d01; // Vector from T1 to T2
        dist = fabs(d01);
      }
    }

    if (fabs(d11) < dist)
    {
      float P11[3];
      for (int d = 0; d < 3; d++)
        P11[d] = V11[d] - N2[d]*d11;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P11,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P11,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P11,V12) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = -N2[i]*d11;
        dist = fabs(d11);
      }
    }

    if (fabs(d21) < dist)
    {
      float P21[3];
      for (int d = 0; d < 3; d++)
        P21[d] = V21[d] - N2[d]*d21;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P21,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P21,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P21,V12) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = -N2[i]*d21;
        dist = fabs(d21);
      }
    }
  }

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d02) == sgn(d12) && sgn(d02) == sgn(d22)) // && fabs(d02) > tol)
  {
    noTouch = true;

    // No intersection; check if projection of points provides closer distance
    if (fabs(d02) < dist)
    {
      float P02[3];
      for (int d = 0; d < 3; d++)
        P02[d] = V02[d] - N1[d]*d02;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P02,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P02,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P02,V11) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = N1[i]*d02;
        dist = fabs(d02);
      }
    }

    if (fabs(d12) < dist)
    {
      float P12[3];
      for (int d = 0; d < 3; d++)
        P12[d] = V12[d] - N1[d]*d12;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P12,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P12,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P12,V11) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = N1[i]*d12;
        dist = fabs(d12);
      }
    }

    if (fabs(d22) < dist)
    {
      float P22[3];
      for (int d = 0; d < 3; d++)
        P22[d] = V22[d] - N1[d]*d22;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P22,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P22,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P22,V11) > 0;

      if (inside)
      {
        for (int i = 0; i < 3; i++)
          minVec[i] = N1[i]*d22;
        dist = fabs(d22);
      }
    }
  }

  // No intersection; return result from edge intersections & plane projections
  if (noTouch)
    return dist;

  // Compute intersection line
  float L[3];
  CROSS(N1, N2, L);
  norm = NORM(L);
  for (int d = 0; d < 3; d++)
    L[d] /= norm;

  float p0 = DOT(L,V01);
  float p1 = DOT(L,V11);
  float p2 = DOT(L,V21);

  float q0 = DOT(L,V02);
  float q1 = DOT(L,V12);
  float q2 = DOT(L,V22);

  // Figure out which point of each triangle is opposite the other two
  int npt1 = (sgn(d01) != sgn(d11)) ? ( (sgn(d11) == sgn(d21)) ? 0 : 1 ) : 2;
  int npt2 = (sgn(d02) != sgn(d12)) ? ( (sgn(d12) == sgn(d22)) ? 0 : 1 ) : 2;

  float s1, s2;
  switch (npt1)
  {
    case 0:
      s1 = p1 + (p0-p1) * (d11 / (d11-d01));
      s2 = p2 + (p0-p2) * (d21 / (d21-d01));
      break;
    case 1:
      s1 = p0 + (p1-p0) * (d01 / (d01-d11));
      s2 = p2 + (p1-p2) * (d21 / (d21-d11));
      break;
    case 2:
      s1 = p0 + (p2-p0) * (d01 / (d01-d21));
      s2 = p1 + (p2-p1) * (d11 / (d11-d21));
      break;
  }

  float t1, t2;
  switch (npt2)
  {
    case 0:
      t1 = q1 + (q0-q1) * (d12 / (d12-d02));
      t2 = q2 + (q0-q2) * (d22 / (d22-d02));
      break;
    case 1:
      t1 = q0 + (q1-q0) * (d02 / (d02-d12));
      t2 = q2 + (q1-q2) * (d22 / (d22-d12));
      break;
    case 2:
      t1 = q0 + (q2-q0) * (d02 / (d02-d22));
      t2 = q1 + (q2-q1) * (d12 / (d12-d22));
      break;
  }

  s1 = (fabs(s1) < 1e-10f) ? 0 : s1;
  s2 = (fabs(s2) < 1e-10f) ? 0 : s2;
  t1 = (fabs(t1) < 1e-10f) ? 0 : t1;
  t2 = (fabs(t2) < 1e-10f) ? 0 : t2;

  if (s1 > s2)
    swap(s1,s2);

  if (t1 > t2)
    swap(t1,t2);

  if (s2 < t1 || t2 < s1)
  {
    // No overlap; return min of dt*L and minDist
    float dt = fmin(fabs(t1-s2), fabs(s1-t2));
    float dl = 0;
    for (int d = 0; d < 3; d++)
      dl += (dt*L[d])*(dt*L[d]);
    dl = sqrt(dl);

    if (dl < dist)
    {
      dist = dl;
      for (int i = 0; i < 3; i++)
        minVec[i] = sgn(t1-s2)*dt*L[i]; // Ensure vec is T1 -> T2
    }

    return dist;
  }

  return 0.f;
}

/*! Modified Moller triangle-triangle intersection algorithm
 *  Determines if triangles intersect, or returns an approximate minimum
 *  distance between them otherwise */
static
__device__
float triTriDistance(float* __restrict__ T1, float* __restrict__ T2, float tol)
{
  float dist = 1e15f;
  float vec[3];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      int i2 = (i+1) % 3;
      int j2 = (j+1) % 3;
      float D = lineSegmentDistance(&T1[3*i], &T1[3*i2], &T2[3*j], &T2[3*j2], vec);

      if (D < dist)
      {
        dist = D;
      }
    }
  }

  // Pointers to points
  const float* V01 = T1;
  const float* V11 = T1+3;
  const float* V21 = T1+6;

  const float* V02 = T2;
  const float* V12 = T2+3;
  const float* V22 = T2+6;

  float N1[3], N2[3];

  // Plane for Triangle 1
  CROSS4(V11,V01, V21,V01, N1);

  float norm = NORM(N1);

  // Plane for Triangle 2
  for (int d = 0; d < 3; d++)
    N1[d] /= norm;

  float d1 = -DOT(N1,V01);

  CROSS4(V12,V02, V22,V02, N2);

  norm = NORM(N2);

  for (int d = 0; d < 3; d++)
    N2[d] /= norm;

  float d2 = -DOT(N2,V02);

  // Signed distances of T1's vertices to T2's plane
  float d01 = DOT(N2,V01) + d2;
  float d11 = DOT(N2,V11) + d2;
  float d21 = DOT(N2,V21) + d2;

  float d02 = DOT(N1,V02) + d1;
  float d12 = DOT(N1,V12) + d1;
  float d22 = DOT(N1,V22) + d1;

  // Round values near 0 to 0
  d01 = (fabs(d01) < 1e-7f) ? 0 : d01;
  d11 = (fabs(d11) < 1e-7f) ? 0 : d11;
  d21 = (fabs(d21) < 1e-7f) ? 0 : d21;

  d02 = (fabs(d02) < 1e-7f) ? 0 : d02;
  d12 = (fabs(d12) < 1e-7f) ? 0 : d12;
  d22 = (fabs(d22) < 1e-7f) ? 0 : d22;

  if (fabs(d01) + fabs(d11) + fabs(d21) < 3*tol ||
      fabs(d02) + fabs(d12) + fabs(d22) < 3*tol)
  {
    // Approximately coplanar; check if one triangle is inside the other /

    // Check if a point in T1 is inside T2
    bool inside = true;
    inside = inside && DOTCROSS4(N2, V12,V02, V01,V02) > 0;
    inside = inside && DOTCROSS4(N2, V02,V22, V01,V22) > 0;
    inside = inside && DOTCROSS4(N2, V22,V12, V01,V12) > 0;

    if (inside) return 0.;

    // Check if a point in T2 is inside T1
    inside = true;
    inside = inside && DOTCROSS4(N1, V11,V01, V02,V01) > 0;
    inside = inside && DOTCROSS4(N1, V01,V21, V02,V21) > 0;
    inside = inside && DOTCROSS4(N1, V21,V11, V02,V11) > 0;

    if (inside) return 0.;
  }

  bool noTouch = false;

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d01) == sgn(d11) && sgn(d01) == sgn(d21)) // && fabs(d01) > tol)
  {
    noTouch = true;

    // No intersection; check if projection of points provides closer distance
    if (fabs(d01) < dist)
    {
      float P01[3];
      for (int d = 0; d < 3; d++)
        P01[d] = V01[d] - N2[d]*d01;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P01,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P01,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P01,V12) > 0;

      if (inside)
      {
        dist = fabs(d01);
      }
    }

    if (fabs(d11) < dist)
    {
      float P11[3];
      for (int d = 0; d < 3; d++)
        P11[d] = V11[d] - N2[d]*d11;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P11,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P11,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P11,V12) > 0;

      if (inside)
      {
        dist = fabs(d11);
      }
    }

    if (fabs(d21) < dist)
    {
      float P21[3];
      for (int d = 0; d < 3; d++)
        P21[d] = V21[d] - N2[d]*d21;
      bool inside = true;
      inside = inside && DOTCROSS4(N2, V12,V02, P21,V02) > 0;
      inside = inside && DOTCROSS4(N2, V02,V22, P21,V22) > 0;
      inside = inside && DOTCROSS4(N2, V22,V12, P21,V12) > 0;

      if (inside)
      {
        dist = fabs(d21);
      }
    }
  }

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d02) == sgn(d12) && sgn(d02) == sgn(d22)) // && fabs(d02) > tol)
  {
    noTouch = true;

    // No intersection; check if projection of points provides closer distance
    if (fabs(d02) < dist)
    {
      float P02[3];
      for (int d = 0; d < 3; d++)
        P02[d] = V02[d] - N1[d]*d02;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P02,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P02,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P02,V11) > 0;

      if (inside)
      {
        dist = fabs(d02);
      }
    }

    if (fabs(d12) < dist)
    {
      float P12[3];
      for (int d = 0; d < 3; d++)
        P12[d] = V12[d] - N1[d]*d12;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P12,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P12,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P12,V11) > 0;

      if (inside)
      {
        dist = fabs(d12);
      }
    }

    if (fabs(d22) < dist)
    {
      float P22[3];
      for (int d = 0; d < 3; d++)
        P22[d] = V22[d] - N1[d]*d22;
      bool inside = true;
      inside = inside && DOTCROSS4(N1, V11,V01, P22,V01) > 0;
      inside = inside && DOTCROSS4(N1, V01,V21, P22,V21) > 0;
      inside = inside && DOTCROSS4(N1, V21,V11, P22,V11) > 0;

      if (inside)
      {
        dist = fabs(d22);
      }
    }
  }

  // No intersection; return result from edge intersections & plane projections
  if (noTouch)
    return dist;

  // Compute intersection line
  float L[3];
  CROSS(N1, N2, L);
  norm = NORM(L);
  for (int d = 0; d < 3; d++)
    L[d] /= norm;

  float p0 = DOT(L,V01);
  float p1 = DOT(L,V11);
  float p2 = DOT(L,V21);

  float q0 = DOT(L,V02);
  float q1 = DOT(L,V12);
  float q2 = DOT(L,V22);

  // Figure out which point of each triangle is opposite the other two
  int npt1 = (sgn(d01) != sgn(d11)) ? ( (sgn(d11) == sgn(d21)) ? 0 : 1 ) : 2;
  int npt2 = (sgn(d02) != sgn(d12)) ? ( (sgn(d12) == sgn(d22)) ? 0 : 1 ) : 2;

  float s1, s2;
  switch (npt1)
  {
    case 0:
      s1 = p1 + (p0-p1) * (d11 / (d11-d01));
      s2 = p2 + (p0-p2) * (d21 / (d21-d01));
      break;
    case 1:
      s1 = p0 + (p1-p0) * (d01 / (d01-d11));
      s2 = p2 + (p1-p2) * (d21 / (d21-d11));
      break;
    case 2:
      s1 = p0 + (p2-p0) * (d01 / (d01-d21));
      s2 = p1 + (p2-p1) * (d11 / (d11-d21));
      break;
  }

  float t1, t2;
  switch (npt2)
  {
    case 0:
      t1 = q1 + (q0-q1) * (d12 / (d12-d02));
      t2 = q2 + (q0-q2) * (d22 / (d22-d02));
      break;
    case 1:
      t1 = q0 + (q1-q0) * (d02 / (d02-d12));
      t2 = q2 + (q1-q2) * (d22 / (d22-d12));
      break;
    case 2:
      t1 = q0 + (q2-q0) * (d02 / (d02-d22));
      t2 = q1 + (q2-q1) * (d12 / (d12-d22));
      break;
  }

  s1 = (fabs(s1) < 1e-7f) ? 0 : s1;
  s2 = (fabs(s2) < 1e-7f) ? 0 : s2;
  t1 = (fabs(t1) < 1e-7f) ? 0 : t1;
  t2 = (fabs(t2) < 1e-7f) ? 0 : t2;

  if (s1 > s2)
    swap(s1,s2);

  if (t1 > t2)
    swap(t1,t2);

  if (s2 < t1 || t2 < s1)
  {
    // No overlap; return min of dt*L and minDist
    float dt = fmin(fabs(t1-s2), fabs(s1-t2));
    float dl = 0;
    for (int d = 0; d < 3; d++)
      dl += (dt*L[d])*(dt*L[d]);
    dl = sqrt(dl);

    if (dl < dist)
      dist = dl;

    return dist;
  }

  return 0.f;
}

static
__device__ __forceinline__
dPointf faceNormal(const float* xv)
{
  /* Assuming nodes of face ordered CCW such that right-hand rule gives
     * outward normal */

  // Triangle #1
  dPointf pt0 = dPointf(&xv[0]);
  dPointf pt1 = dPointf(&xv[3]);
  dPointf pt2 = dPointf(&xv[6]);
  dPointf norm1 = (pt1-pt0).cross(pt2-pt0);           // Face normal vector

  // Triangle #2
  pt1 = dPointf(&xv[9]);
  dPointf norm2 = (pt2-pt0).cross(pt1-pt0);

  // Average the two triangle's normals
  dPointf norm = 0.5*(norm1+norm2);

  return (norm / norm.norm());
}

//! Calculate the distance between a linear triangle and a curved face
template<int nSideF>
__device__
float intersectionCheckOne(dMeshBlock &mb, const float* __restrict__ fxv,
    float* __restrict__ minVec, float* TC)
{
  /* --- Prerequisites --- */

  const int sorderF = nSideF-1;

  float tol = 1e-7f;
  float TF[9];
  float minDist = BIG_DOUBLE;
  minVec[0] = minDist;
  minVec[1] = minDist;
  minVec[2] = minDist;

  /* Only 3 cases possible:
   * 1) Face entirely contained within element
   * 2) Face intersects with element's boundary
   * 3) Face and element do not intersect
   */

  for (int M = 0; M < sorderF; M++)
  {
    for (int N = 0; N < sorderF; N++)
    {
      int m0 = M + nSideF*N;
      int TriPtsF[2][3] = {{m0, m0+1, m0+nSideF+1}, {m0, m0+nSideF+1, m0+nSideF}};
      for (int m = 0; m < 2; m++)
        for (int n = 0; n < 3; n++)
          TriPtsF[m][n] = mb.ijk2gmsh_quad[TriPtsF[m][n]];

      // Intersection check between element face tris & cutting-face tris
      for (int j = 0; j < 2; j++)
      {
        for (int p = 0; p < 3; p++)
        {
          int ipt = TriPtsF[j][p];
          for (int d = 0; d < 3; d++)
            TF[3*p+d] = fxv[3*ipt+d];
        }

        float vec[3];
        float dist = triTriDistanceVec(TF, TC, vec, tol);

        if (dist < tol)
          return 0.;

        if (dist < minDist)
        {
          for (int d = 0; d < 3; d++)
            minVec[d] = vec[d];
          minDist = dist;
        }
      }
    }
  }

  return minDist;
}

//! Calc distance from near half of a hex to a face, assuming both linear
__device__ __forceinline__
float intersectionCheckLinear(const float* __restrict__ fxv,
    const float* __restrict__ exv, char &cornerOut)
{
  /* --- Prerequisites --- */

  // NOTE: Gmsh ordering  |  btm,top,left,right,front,back
  const char TriPts[12][3] = {{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,3,7},{0,7,4},
                       {1,2,6},{1,6,5},{0,4,5},{0,5,1},{3,2,6},{3,7,6}};

  float tol = 1e-7f;
  float TC[9], TF[9];
  float minDist = 1e15f;

  float xcf[3];
  cuda_funcs::getCentroid<3,4>(fxv,xcf);

  // Find nearest corner of element to face; check only that half of element
  int corner = -1;
  for (int i = 0; i < 8; i++)
  {
    float dist = 0.f;
    for (int d = 0; d < 3; d++)
      dist += (exv[3*i+d] - xcf[d]) * (exv[3*i+d] - xcf[d]);

    if (dist < minDist)
    {
      minDist = dist;
      corner = i;
    }
  }

  cornerOut = corner;

  // Faces 0 or 1, 2 or 3, and 4 or 5 (btm or top, L or R, etc.)
  const char fList[3] = {(char)(corner / 4), (char) (((corner + 1)%4) / 2 + 2), (char)(((corner%4) / 2) + 4)};

  // 3) Check those faces of element for intersection with face
  for (int F = 0; F < 3; F++)
  {
    char f = fList[F];
    // Get triangles for the sub-hex of the larger curved hex
    for (int i = 2*f; i < 2*f+2; i++)
    {
      for (int p = 0; p < 3; p++)
      {
        char ipt = TriPts[i][p];
        for (int d = 0; d < 3; d++)
          TC[3*p+d] = exv[3*ipt+d];
      }

      const char TriPtsF[2][3] = {{0, 1, 3}, {1, 2, 3}};

      // Intersection check between element face tris & cutting-face tris
      for (int j = 0; j < 2; j++)
      {
        for (int p = 0; p < 3; p++)
        {
          char ipt = TriPtsF[j][p];
          for (int d = 0; d < 3; d++)
            TF[3*p+d] = fxv[3*ipt+d];
        }

        float dist = triTriDistance(TF, TC, tol);

        if (dist < tol)
          return 0.;

        minDist = (dist < minDist) ? dist : minDist;
      }
    }
  }

  return minDist;
}

//! Perform initial ele-face distance sorting using oriented bounding boxes
__global__
void cuttingPass0B(dvec<float> eleBbox, dvec<float> eleXC, dvec<float> faceXC,
    dvec<float> cutFaces, int nEles, int nFaces, int nvertf, dvec<int> filt_faces, dvec<float> outDist)
{
  const int IC = blockIdx.x * blockDim.x + threadIdx.x;
  const int F = blockIdx.y * blockDim.y + threadIdx.y;

  if (IC >= nEles || F >= nFaces) return;

  const int nDims = 3;

  const int ff = filt_faces[F];
  const int stride = nDims*nvertf;

  /// Alternatively: compute distance to ele bbox for all 4 corner points
  float pts[12];
  for (int i = 0; i < 4; i++)
  {
    float pt[3];
    for (int j = 0; j < 3; j++)
      pt[j] = cutFaces[ff*stride+3*i+j];

    for (int j = 0; j < 3; j++)
    {
      pts[3*i+j] = 0;
      for (int k = 0; k < 3; k++)
      {
        pts[3*i+j] += eleBbox[16*IC+3*k+j] * pt[k];
      }
    }
  }

  float bboxF[6];
  cuda_funcs::getBoundingBox<3,4>(pts, bboxF);

  float dist = 0;
  for (int i = 0; i < 3; i++)
    dist += (eleXC[3*IC+i]-faceXC[3*F+i])*(eleXC[3*IC+i]-faceXC[3*F+i]);
  dist = sqrt(dist);

  outDist[nEles*F+IC] = .01f*dist + cuda_funcs::boundingBoxDist<3>(bboxF, &eleBbox[16*IC+9]);
}

__global__
void sortFaces0(dvec<int> faceList, int nEles, int nFaces, dvec<float> distList,
    dvec<int> outFaces)
{
  const int IC = blockDim.x * blockIdx.x + threadIdx.x;

  if (IC >= nEles) return;

  float dists[NF1];
  int faces[NF1];

  for (int i = 0; i < NF1; i++)
  {
    dists[i] = BIG_FLOAT;
    faces[i] = -1;
  }

  for (int F = 0; F < nFaces; F++)
  {
    int ind = NF1-1;

    float Dist = distList[nEles*F+IC];
    if (Dist > dists[ind]) continue;

    dists[ind] = Dist;
    faces[ind] = faceList[F];
    while (ind > 0 && dists[ind] < dists[ind-1])
    {
      swap(dists[ind], dists[ind-1]);
      swap(faces[ind], faces[ind-1]);
      ind--;
    }
  }

  for (int i = 0; i < NF1; i++)
  {
    outFaces[nEles*i+IC] = faces[i];
  }
}

__global__
void cuttingPass1(dMeshBlock mb, dvec<int> filt_eles, int nEles,
    dvec<float> cutFaces, int nvertf, int nFaces, dvec<int> checkFaces,
    dvec<char> outCorner, dvec<float> outDist)
{
  const int IC = blockIdx.x * blockDim.x + threadIdx.x;
  const int F = blockIdx.y * blockDim.y + threadIdx.y;

  if (IC >= nEles || F >= nFaces) return;

  const int nDims = 3;

  const int ic = filt_eles[IC];  // Get filtered cell ID
  const int ff = checkFaces[nEles*F+IC];
  const int stride = nDims*nvertf;

  // Load up the cell nodes into an array
  float xv[8*nDims];
  for (int i = 0; i < 8*nDims; i++)
  {
    int d = i % 3;
    int v = i / 3;
    xv[i] = mb.coord[ic+mb.ncells*(d+nDims*v)]; /// NOTE: 'row-major' ZEFR layout
  }

  // Each thread will check against 1 face
  outDist[nEles*F+IC] = intersectionCheckLinear(&cutFaces[ff*stride], xv, outCorner[nEles*F+IC]);
}

__global__
void sortFaces(int nEles, int nFaces, dvec<float> distList,
    dvec<char> cornerList, dvec<int> outFaces, dvec<char> outCorners)
{
  const int IC = blockDim.x * blockIdx.x + threadIdx.x;

  if (IC >= nEles) return;

  float dists[NF2];
  int faces[NF2];
  char corners[NF2];

  for (int i = 0; i < NF2; i++)
  {
    dists[i] = BIG_FLOAT;
    faces[i] = -1;
    corners[i] = -1;
  }

  for (int F = 0; F < nFaces; F++)
  {
    int ind = NF2-1;

    float Dist = distList[nEles*F+IC];
    if (Dist > dists[ind]) continue;

    dists[ind] = Dist;
    faces[ind] = outFaces[nEles*F+IC];
    corners[ind] = cornerList[nEles*F+IC];
    while (ind > 0 && dists[ind] < dists[ind-1])
    {
      swap(dists[ind], dists[ind-1]);
      swap(faces[ind], faces[ind-1]);
      swap(corners[ind], corners[ind-1]);
      ind--;
    }
  }

  for (int i = 0; i < NF2; i++)
  {
    outFaces[nEles*i+IC] = faces[i];
    outCorners[nEles*i+IC] = corners[i];
  }
}

__global__
void getFaceCentroids(dvec<float> cutFaces, dvec<int> faceList, int nFaces,
    int nvertf, dvec<float> faceCentroid)
{
  int F = blockIdx.x * blockDim.x + threadIdx.x;

  if (F >= nFaces) return;

  int ff = faceList[F];

  cuda_funcs::getCentroid<3,4>(&cutFaces[3*nvertf*ff], &faceCentroid[3*F]);
}

__global__
void getElementCentroids(dMeshBlock mb, dvec<int> eleList, int nEles,
    dvec<float> eleCentroid)
{
  int IC = blockIdx.x * blockDim.x + threadIdx.x;

  if (IC >= nEles) return;

  int ic = eleList[IC];

  float xv[8*3];  // Only concerning ourselves with linear portion of ele
  for (int i = 0; i < 8; i++)
    for (int d = 0; d < 3; d++)
      xv[3*i+d] = mb.coord[ic+mb.ncells*(d+3*i)];

  cuda_funcs::getCentroid<3,8>(xv,&eleCentroid[3*IC]);
}

__global__
void getElementBoundingBoxes(dMeshBlock mb, dvec<int> eleList, int nEles,
    dvec<float> eleBbox)
{
  int IC = blockIdx.x * blockDim.x + threadIdx.x;

  if (IC >= nEles) return;

  int ic = eleList[IC];

  float xv[8*3];  // Only concerning ourselves with linear portion of ele
  for (int i = 0; i < 8; i++)
    for (int d = 0; d < 3; d++)
      xv[3*i+d] = mb.coord[ic+mb.ncells*(d+3*i)];

  cuda_funcs::getBoundingBox<3,8>(xv,&eleBbox[6*IC]); /// TODO: swap rows/cols in eleBbox for better coalesced access
}

__global__
void getElementOrientedBoundingBoxes(dMeshBlock mb, dvec<int> eleList, int nEles,
    dvec<float> eleBbox)
{
  int IC = blockIdx.x * blockDim.x + threadIdx.x;

  if (IC >= nEles) return;

  int ic = eleList[IC];

  float xv[8*3];  // Only concerning ourselves with linear portion of ele
  for (int i = 0; i < 8; i++)
    for (int d = 0; d < 3; d++)
      xv[3*i+d] = mb.coord[ic+mb.ncells*(d+3*i)];

  getOBB(xv, 8, &eleBbox[16*IC]);

  // Use last entry for 'href' [average length of oriented bounding box]
  eleBbox[16*IC+15] = ( (eleBbox[16*IC+12] - eleBbox[16*IC+9]) +
                        (eleBbox[16*IC+13] - eleBbox[16*IC+10]) +
                        (eleBbox[16*IC+14] - eleBbox[16*IC+11]) ) / 3.;
}

template<int nSideF>
__global__
void cuttingPass2(dMeshBlock mb, dvec<float> cutFaces, dvec<int> checkFaces,
    dvec<int> list, int nEles, int nSideC, dvec<float> eleBbox,
    dvec<char> corners, dvec<float> outDist, dvec<float> outVec)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int FID = threadIdx.y;

  const int nDims = 3;
  const int nvertf = nSideF*nSideF;
  const int stride = nDims*nvertf;
  const int sOrderC = nSideC - 1;

  const int nQuadFace = sOrderC*sOrderC;
  const int nTriFace = 2*nQuadFace;
  const int nTri = 3*nTriFace;

  /* --- Get our specific element & sub-triangle of element --- */

  const int IC = tid / (3 * nTriFace);
  const int T = tid % (3 * nTriFace);

  const unsigned char q = (T % (nTriFace)) / 2;
  const unsigned char t = (T % (nTriFace)) % 2;
  char F = T / nTriFace;

  if (IC >= nEles) return;

  const int ic = list[IC];  // Get filtered cell ID
  const int ff = checkFaces[nEles*FID+IC];

  if (ff < 0)
    return;

  /* ---- Check against our reduced list of faces ---- */

  /// Or TODO on storing href in global memory instead...
  const float href = eleBbox[16*IC+15];
  const float dtol = .2*href;

  // Only checking half the element's faces; figure out which ones
  const char corner = corners[nEles*FID+IC];

  switch (F)
  {
    case 0:
      F = corner / 4;
      break;
    case 1:
      F = ((corner + 1)%4) / 2 + 2;
      break;
    case 2:
      F = (((corner%4) / 2) + 4);
      break;
  }

  /* ---- Get our specific triangle ---- */

  // NOTE: Structured ordering  |  btm,top,left,right,front,back
  const char TriPts[12][3] = {{0,1,3},{0,3,2},{4,7,5},{4,6,7},{0,2,6},{0,6,4},
    {1,3,7},{1,7,5},{0,4,5},{0,5,1},{2,3,7},{2,6,7}};

  int I, J, K;
  switch (F)
  {
    case 0: // Bottom
      I = q / sOrderC;
      J = q % sOrderC;
      K = 0;
      break;
    case 1: // Top
      I = q / sOrderC;
      J = q % sOrderC;
      K = sOrderC - 1;
      break;
    case 2: // Left
      I = 0;
      J = q / sOrderC;
      K = q % sOrderC;
      break;
    case 3: // Right
      I = sOrderC - 1;
      J = q / sOrderC;
      K = q % sOrderC;
      break;
    case 4: // Front
      I = q / sOrderC;
      J = 0;
      K = q % sOrderC;
      break;
    case 5: // Back
      I = q / sOrderC;
      J = sOrderC - 1;
      K = q % sOrderC;
      break;
  }

  int i0 = I+nSideC*(J+nSideC*K);
  int j0 = i0 + nSideC*nSideC;
  int lin2curv[8] = {i0, i0+1, i0+nSideC, i0+nSideC+1, j0, j0+1, j0+nSideC, j0+nSideC+1};
  for (int i = 0; i < 8; i++)
    lin2curv[i] = mb.ijk2gmsh[lin2curv[i]];

  float TC[9];
  for (int p = 0; p < 3; p++)
  {
    int ipt = lin2curv[TriPts[2*F+t][p]];
    for (int d = 0; d < 3; d++)
      TC[3*p+d] = mb.coord[ic+mb.ncells*(d+nDims*ipt)]; /// NOTE: 'row-major' ZEFR layout
  }

  // Find distance from face to cell
  /// NOTE: ignoring case of face entirely inside cell, since any valid grid
  /// will also have a different face which intersects its boundary
  dPointf vec;
  float myDist = intersectionCheckOne<nSideF>(mb, &cutFaces[ff*stride], &vec[0], TC);
  vec /= vec.norm();

  if (myDist < dtol) // They intersect
  {
    myDist = 0.;
  }

  // Write out results to global memory for future reduction
  outDist[T+nTri*(IC+nEles*FID)] = myDist;
  for (int i = 0; i < 3; i++)
  {
    //outNorm[IC+nEles*(i+3*FID)] = myNorm[i];
    outVec[T+nTri*(IC+nEles*(FID+NF2*i))] = vec[i];
  }
}

__global__
void getMinDist(dvec<float> dists, dvec<float> vecs, int nEles, int nTri)
{
  const int IC = blockDim.x * blockIdx.x + threadIdx.x;
  const int F = threadIdx.y;

  if (IC >= nEles) return;

  // Find minimum tri-face distance for this face/element
  float minDist = BIG_FLOAT;
  float myVec[3] = {0.0f};
  for (int i = 0; i < nTri; i++)
  {
    float dist = dists[i+nTri*(IC+nEles*F)];

    if (dist < minDist)
    {
      minDist = dist;
      for (int d = 0; d < 3; d++)
        myVec[d] = vecs[i+nTri*(IC+nEles*(F+NF2*d))];
    }
  }

  // NOTE: Assuming NF2 always <= nTri
  dists[F+nTri*IC] = minDist;
  for (int d = 0; d < 3; d++)
    vecs[F+nTri*(IC+nEles*NF2*d)] = myVec[d];
}

__global__
void getFinalFlag(dvec<int> eleList, dvec<int> checkFaces,
    dvec<float> cutFaces, dvec<float> eleBbox, int nEles, int nvertf, int nTri,
    dvec<int> cutFlag, dvec<float> dists, dvec<float> vecs, int cutType)
{
  const int IC = blockDim.x * blockIdx.x + threadIdx.x;

  if (IC >= nEles) return;

  const int ic = eleList[IC];

  const float dtol = .1*eleBbox[16*IC+15];

  // Find nearest face distance for this element

  int nMin = 0;
  char myFlag = DC_UNASSIGNED;
  float minDist = BIG_FLOAT;
  float myDot = 0.f;
  dPointf myNorm;

  for (int F = 0; F < NF2; F++)
  {
    int ff = checkFaces[nEles*F+IC];
    dPointf norm = faceNormal(&cutFaces[ff*nvertf*3]);
    if (cutType == 0) norm *= -1;

    float dist = dists[F+nTri*IC];

    dPointf vec;
    for (int d = 0; d < 3; d++)
      vec[d] = vecs[F+nTri*(IC+nEles*NF2*d)];

    if (dist < dtol) // They intersect
    {
      myFlag = DC_CUT;
      minDist = 0.;
    }
    else if (myFlag == DC_UNASSIGNED || dist < (minDist - .02*dtol))
    {
      // Unflagged cell, or have a closer face to use
      minDist = dist;
      myNorm = norm;
      myDot = myNorm*vec;

      nMin = 1;

      if (myDot < 0)
        myFlag = DC_HOLE; // outwards normal = inside cutting surface
      else
        myFlag = DC_NORMAL;
    }
    else if (fabs(dist - minDist) <= .02*dtol)
    {
      // Approx. same dist. to two faces; avg. their normals to decide
      minDist = dist;
      for (int d = 0; d < 3; d++)
        myNorm[d] = (nMin*myNorm[d] + norm[d]) / (nMin + 1.);
      nMin++;

      myDot = myNorm*vec;

      if (myDot < 0)
        myFlag = DC_HOLE; // outwards normal = inside cutting surface
      else
        myFlag = DC_NORMAL;
    }
  }

  if (myFlag == DC_CUT)
    myFlag = (cutType != 0) ? DC_HOLE : DC_NORMAL;

  // Write out final result
  cutFlag[ic] = myFlag;
}

/*! Remove all elements which do not intersect with cut group's bbox from
 *  consideration (obviously do not intersect) */
template<int nvert>
__global__
void filterElements(dMeshBlock mb, dvec<double> cut_bbox, dvec<int> filt,
    dvec<int> cutFlag, dvec<int> nfilt, dvec<float> bboxOut, int cutType)
{
  const unsigned int ic = blockIdx.x * blockDim.x + threadIdx.x;

  // Initialize nfilt to 0; will be atomically added to at end
  if (ic == 0)
  {
    nfilt[0] = 0;
    for (int d = 0; d < 3; d++)
    {
      bboxOut[d]   =  1.e10f;
      bboxOut[d+3] = -1.e10f;
    }
  }

  __shared__ float bboxF[6];

  for (int i = threadIdx.x; i < 6; i += blockDim.x)
    bboxF[i] = (float)cut_bbox[i];

  __syncthreads();

  if (ic >= mb.ncells) return;

  // Set all cell flags initially to DC_NORMAL (filtered cells will remain 'NORMAL')
  cutFlag[ic] = DC_NORMAL;

  float href = .005f/3.f*(bboxF[3]-bboxF[0]+bboxF[4]-bboxF[1]+bboxF[5]-bboxF[2]);

  // Get element nodes
  float xv[nvert*3];
  for (int i = 0; i < nvert; i++)
    for (int d = 0; d < 3; d++)
      xv[3*i+d] = (float)mb.coord[ic+mb.ncells*(d+3*i)];

  // Get element bounding box
  float bboxC[6], xc[3];
  cuda_funcs::getBoundingBox<3,nvert>(xv, bboxC);

  bool checkH = false;

  if (cutType == 1)
  {
    // Wall boundary - set as hole if centroid inside wall
    cuda_funcs::getCentroid<3,nvert>(xv, xc);
    char tag = cuda_funcs::checkHoleMap(xc, mb.hm_sam.data(), mb.hm_nx.data(), mb.hm_extents.data());
    checkH = (tag != 1);
  }
  else
  {
    // Overset boundary - only set as hole if entirely inside hole map
    for (int i = 0; i < 8; i++)
    {
      for (int d = 0; d < 3; d++)
        xc[d] = xv[3*i+d];

      if (mb.rrot) // Transform xc to hole map's coordinate system
      {
        double x2[3] = {0.,0.,0.};
        for (int d1 = 0; d1 < 3; d1++)
          for (int d2 = 0; d2 < 3; d2++)
            x2[d1] += mb.Rmat[d1+3*d2]*(xc[d2]-mb.offset[d2]); //! TODO: include Rmat from other grid

        for (int d = 0; d < 3; d++)
          xc[d] = x2[d];
      }

      char tag = cuda_funcs::checkHoleMap(xc, mb.hm_sam.data(), mb.hm_nx.data(), mb.hm_extents.data());
      checkH = checkH || (tag != 1);
    }
  }

  bool checkB = cuda_funcs::boundingBoxCheck<3>(bboxC, bboxF, href);

  // If filtering element due to being completely inside hole region, tag as hole
  if (!checkH)
    cutFlag[ic] = DC_HOLE;

  if (checkH && checkB)
  {
    int ind = atomicAggInc(&nfilt[0]);
    filt[ind] = ic;
    for (int d = 0; d < 3; d++)
    {
      atomicMinf(&bboxOut[d], bboxC[d]);
      atomicMaxf(&bboxOut[d+3], bboxC[d+3]);
    }
  }
}

/*! Remove all cutting faces which do not intersect this rank's reduced bbox
 *  from consideration (obviously do not intersect) */
template<int nvertf>
__global__
void filterFaces(dvec<float> ele_bbox, int nCut,
    dvec<float> cutFaces, dvec<int> filt, dvec<int> nfilt, dvec<float> bboxOut)
{
  const unsigned int ff = blockIdx.x * blockDim.x + threadIdx.x;

  // Initialize nfilt to 0; will be atomically added to at end
  if (ff == 0)
  {
    nfilt[1] = 0;
    for (int d = 0; d < 3; d++)
    {
      bboxOut[d]   =  1e10f;
      bboxOut[d+3] = -1e10f;
    }
  }

  __shared__ float bboxE[6]; // Global bounding box of all filtered elements

  for (int i = threadIdx.x; i < 6; i += blockDim.x)
    bboxE[i] = ele_bbox[i];

  __syncthreads();

  if (ff >= nCut) return;

  float href = .01f/3.f*(bboxE[3]-bboxE[0]+bboxE[4]-bboxE[1]+bboxE[5]-bboxE[2]);

  // Get face nodes
  float fxv[nvertf*3];
  for (int i = 0; i < nvertf; i++)
    for (int d = 0; d < 3; d++)
      fxv[3*i+d] = (float)cutFaces[(ff*nvertf+i)*3+d];

  // Get face bounding box
  float bboxF[6];
  cuda_funcs::getBoundingBox<3,nvertf>(fxv, bboxF);

  bool checkB = cuda_funcs::boundingBoxCheck<3>(bboxF, bboxE, href);

  if (checkB)
  {
    int ind = atomicAggInc(&nfilt[1]);
    filt[ind] = ff;
    for (int d = 0; d < 3; d++)
    {
       atomicMinf(&bboxOut[d], bboxF[d]);
       atomicMaxf(&bboxOut[d+3], bboxF[d+3]);
    }
  }
}

void dMeshBlock::directCut(double* cutFaces_h, int nCut, int nvertf, double *cutBbox_h, int* cutFlag, int cutType)
{
  // Setup cutMap TODO: create initialization elsewhere?
  cutFlag_d.resize(ncells);
  filt_eles.resize(ncells);
  filt_faces.resize(nCut);

  std::vector<float> cutFaces_hf(nCut*nvertf*nDims);
  for (int i = 0; i < nvertf*nDims*nCut; i++)
    cutFaces_hf[i] = (float)cutFaces_h[i];

  dvec<float> cutFaces;
  cutFaces.assign(cutFaces_hf.data(), nCut*nvertf*nDims);

  dvec<double> cutBbox_d;
  cutBbox_d.assign(cutBbox_h, 2*nDims);

  if (ijk2gmsh_quad.size() != nvertf)
  {
    auto ijk2gmsh_quad_h = tg_funcs::structured_to_gmsh_quad(nvertf);
    ijk2gmsh_quad.assign(ijk2gmsh_quad_h.data(), nvertf);
  }

  /* Filter elements based upon cutting-surface bounding box & Cartesian approx. rep. */

  hvec<int> nfilt_h;
  dvec<int> nfilt_d;
  nfilt_h.resize(2);
  nfilt_h[0] = 0;  nfilt_h[1] = 0;
  nfilt_d.assign(nfilt_h.data(), nfilt_h.size());

  ele_bbox.resize(6);
  face_bbox.resize(6);

  int threads = 128;
  int blocks = (ncells + threads - 1) / threads;

  switch(nvert)
  {
    case 8:
      filterElements<8><<<blocks, threads, 6*sizeof(float)>>>(*this, cutBbox_d, filt_eles, cutFlag_d, nfilt_d, ele_bbox, cutType);
      break;
    case 27:
      filterElements<27><<<blocks, threads, 6*sizeof(float)>>>(*this, cutBbox_d, filt_eles, cutFlag_d, nfilt_d, ele_bbox, cutType);
      break;
    case 64:
      filterElements<64><<<blocks, threads, 6*sizeof(float)>>>(*this, cutBbox_d, filt_eles, cutFlag_d, nfilt_d, ele_bbox, cutType);
      break;
    default:
      printf("nvert = %d\n",nvert);
      ThrowException("nvert case not implemented for filterElements on device");
  }
  check_error();

  if (nCut == 0)
  {
    cuda_copy_d2h(cutFlag_d.data(), cutFlag, ncells);

    cutFaces.free_data();
    cutBbox_d.free_data();
    nfilt_h.free_data();
    nfilt_d.free_data();

    return;
  }

  /* Filter cutting faces by intersection with the filtered elements' bounding box */

  blocks = (nCut + threads - 1) / threads;

  switch(nvertf)
  {
    case 4:
      filterFaces<4><<<blocks, threads, 6*sizeof(float)>>>(ele_bbox, nCut, cutFaces, filt_faces, nfilt_d, face_bbox);
      break;
    case 9:
      filterFaces<9><<<blocks, threads, 6*sizeof(float)>>>(ele_bbox, nCut, cutFaces, filt_faces, nfilt_d, face_bbox);
      break;
    case 16:
      filterFaces<16><<<blocks, threads, 6*sizeof(float)>>>(ele_bbox, nCut, cutFaces, filt_faces, nfilt_d, face_bbox);
      break;
    default:
      printf("nvertf = %d\n",nvertf);
      ThrowException("nvertf case not implemented for filterFaces on device");
  }
  check_error();

  nfilt_h.assign(nfilt_d.data(), 2);
  int nfiltC = nfilt_h[0];
  int nfiltF = nfilt_h[1];

  /* Perform the Direct Cut algorithm on the filtered list of elements & faces */

  int nSideC = std::cbrt(nvert);
  int nTri = 3*2*(nSideC-1)*(nSideC-1);

  if (nTri < NF2) ThrowException("Incompatible nTri / NF2!");

  int nCheck1 = min(NF1, nfiltF);

  dvec<float> eleBbox, cfDist, cfNorm, cfVec;
  dvec<char> cfCorner;
  cfDist.resize(nfiltC*nfiltF);
  cfCorner.resize(nfiltC*NF1);

  dvec<int> checkFaces;
  dvec<char> corners;
  checkFaces.resize(NF1*nfiltC);
  corners.resize(NF2*nfiltC);

  dvec<float> eleXC, faceXC;

  if (nfiltC > 0 && nfiltF > 0)
  {
    /* Pass 0: 'Coarse-Grained' Check using Bounding Boxes
     * Pass 1: 'Medium-Grained' Check using Linear Parts of Elements/Faces
     * Pass 2: 'Finest-Grained' Direct Cut Check
     */

    // Tell CUDA to offer more space for registers over shared memory
    cudaFuncSetCacheConfig(getElementBoundingBoxes, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(sortFaces0, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(sortFaces, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(cuttingPass0B, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(cuttingPass1, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(cuttingPass2<2>, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(cuttingPass2<3>, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(cuttingPass2<4>, cudaFuncCachePreferL1);

    eleBbox.resize(16*nfiltC);

    eleXC.resize(nfiltC*3);
    faceXC.resize(nfiltC*3);

    int ThreadsB = 128;
    int BlocksB = (nfiltC + ThreadsB - 1) / ThreadsB;

    getElementOrientedBoundingBoxes<<<BlocksB,ThreadsB>>>(*this, filt_eles, nfiltC, eleBbox);

    getElementCentroids<<<BlocksB,ThreadsB>>>(*this, filt_eles, nfiltC, eleXC);
    int ThreadsF = 128;
    int BlocksF = (nfiltF + ThreadsF - 1) / ThreadsF;

    getFaceCentroids<<<BlocksF,ThreadsF>>>(cutFaces, filt_faces, nfiltF, nvertf, faceXC);
    check_error();

    // Have each filtered element calculate a rough distance to each filtered face

    dim3 Threads0(16,8);
    dim3 Blocks0( (nfiltC + Threads0.x - 1) / Threads0.x,
                  (nfiltF + Threads0.y - 1) / Threads0.y );

    cuttingPass0B<<<Blocks0,Threads0>>>(eleBbox, eleXC, faceXC, cutFaces, nfiltC, nfiltF, nvertf, filt_faces, cfDist);

    int ThreadsS0 = 128;
    int BlocksS0 = (nfiltC + ThreadsS0 - 1) / ThreadsS0;

    sortFaces0<<<BlocksS0, ThreadsS0>>>(filt_faces, nfiltC, nfiltF, cfDist, checkFaces);
    check_error();

    // Have each filtered element calculate a rough distance to the reduced face list

    dim3 Threads1(32,4);
    dim3 Blocks1( (nfiltC + Threads1.x - 1) / Threads1.x,
                  (nCheck1 + Threads1.y - 1) / Threads1.y );

    cuttingPass1<<<Blocks1,Threads1>>>(*this, filt_eles, nfiltC, cutFaces,
        nvertf, nCheck1, checkFaces, cfCorner, cfDist);

    // Sort the distance lists for each element [down to NF2 faces to check in detail]

    int ThreadsS1 = 128;
    int BlocksS1 = (nfiltC + ThreadsS1 - 1) / ThreadsS1;

    sortFaces<<<BlocksS1, ThreadsS1>>>(nfiltC, nCheck1, cfDist, cfCorner, checkFaces, corners);
    check_error();

    cfNorm.resize(3*nfiltC*nTri*NF2);
    cfVec.resize(3*nfiltC*nTri*NF2);
    cfDist.resize(nfiltC*nTri*NF2);


    dim3 t3(32,NF2);
    int b3 = (nfiltC*nTri + t3.x - 1) / t3.x;

    switch(nvertf)
    {
      case 4:
        cuttingPass2<2><<<b3, t3>>>(*this, cutFaces, checkFaces, filt_eles,
            nfiltC, nSideC, eleBbox, corners, cfDist, cfVec);
        break;

      case 9:
        cuttingPass2<3><<<b3, t3>>>(*this, cutFaces, checkFaces, filt_eles, nfiltC, nSideC, eleBbox, corners, cfDist, cfVec);
        break;

      case 16:
        cuttingPass2<4><<<b3, t3>>>(*this, cutFaces, checkFaces, filt_eles, nfiltC, nSideC, eleBbox, corners, cfDist, cfVec);
        break;

      default:
        printf("nvertFace = %d, nCut = %d\n",nvertf,nCut);
        ThrowException("nvertFace case not implemented for directCut on device");
    }
    check_error();

    int BlocksM = (nfiltC + t3.x - 1) / t3.x;
    getMinDist<<<BlocksM, t3>>>(cfDist, cfVec, nfiltC, nTri);

    getFinalFlag<<<BlocksM, 128>>>(filt_eles, checkFaces, cutFaces, eleBbox,
        nfiltC, nvertf, nTri, cutFlag_d, cfDist, cfVec, cutType);
    check_error();
  }

  cuda_copy_d2h(cutFlag_d.data(), cutFlag, ncells);

  // Free all data allocated in this function

  cfDist.free_data();
  cfCorner.free_data();
  cfNorm.free_data();
  cfVec.free_data();
  eleBbox.free_data();
  eleXC.free_data();
  faceXC.free_data();

  nfilt_d.free_data();
  nfilt_h.free_data();

  corners.free_data();
  checkFaces.free_data();

  cutFaces.free_data();
  cutBbox_d.free_data();
}
