#ifndef CUDA_FUNCS_H
#define CUDA_FUNCS_H

#include "cuda_runtime.h"
#include "error.hpp"

#ifndef BIG_DOUBLE
#define BIG_DOUBLE 1.0e15
#endif

#ifndef DC_CUT
#define DC_HOLE 0
#define DC_UNASSIGNED 1
#define DC_CUT 2
#define DC_NORMAL 3
#endif

template<typename T>
void cuda_malloc(T* &data_d, size_t size)
{
  cudaMalloc((void**)&data_d, size*sizeof(T));
  check_error();
}

//! Allocate pinned memory on host for asynchronous copying
template<typename T>
void cuda_malloc_pinned(T* &data_h, size_t size)
{
  cudaMallocHost((void**)&data_h, size*sizeof(T));
  check_error();
}

template <typename T>
void cuda_copy_h2d(T* data_d, const T* data_h, size_t size)
{
  cudaMemcpy(data_d, data_h, size*sizeof(T), cudaMemcpyHostToDevice);
  check_error();
}

//! Asynchronous Copy (Must use pinned memory)
template <typename T>
void cuda_copy_h2d(T* data_d, const T* data_h, size_t size, cudaStream_t stream)
{
  cudaMemcpyAsync(data_d, data_h, size*sizeof(T), cudaMemcpyHostToDevice, stream);
  check_error();
}

template <typename T>
void cuda_copy_d2h(const T* data_d, T* data_h, size_t size)
{
  cudaMemcpy(data_h, data_d, size*sizeof(T), cudaMemcpyDeviceToHost);
  check_error();
}

//! Asynchronous Copy (Must use pinned memory)
template <typename T>
void cuda_copy_d2h(const T* data_d, T* data_h, size_t size, cudaStream_t stream)
{
  cudaMemcpyAsync(data_h, data_d, size*sizeof(T), cudaMemcpyDeviceToHost, stream);
  check_error();
}

template <typename T>
void cuda_free(T* &data_d)
{
  cudaFree(data_d);
  data_d = NULL;
  check_error();
}

template <typename T>
void cuda_free_pinned(T* &data_h)
{
  cudaFreeHost(data_h);
  data_h = NULL;
  check_error();
}

template <typename T>
__host__ __device__ __forceinline__
void swap_d(T& a, T& b)
{
  T c(a);
  a = b; b = c;
}

template<typename T>
class dvec
{
private:
  int size_ = 0;
  int max_size_ = 0;
  T *data_ = NULL;
  bool allocated = false;

public:
  __host__ __device__
  dvec(void) { }

  __host__ __device__
  ~dvec(void);

  __host__ __device__
  int size(void) const { return size_; }

  __host__ __device__
  int capacity(void) const { return max_size_; }

  __host__ __device__
  T* data(void) { return data_; }

  __host__
  void resize(int size);

  __host__
  void assign(T* data_h, int size, cudaStream_t *stream = NULL);

  __host__
  void free_data(void);

  __device__
  T operator[](int ind) const;

  __device__
  T& operator[](int ind);
};

template<typename T>
__host__ __device__
dvec<T>::~dvec(void)
{
#ifndef __CUDA_ARCH__
//  free_data();
#endif
}

template<typename T>
__device__
T dvec<T>::operator[](int ind) const
{
  return data_[ind];
}

template<typename T>
__device__
T& dvec<T>::operator[](int ind)
{
  return data_[ind];
}

template<typename T>
__host__
void dvec<T>::resize(int size)
{
  if (allocated)
  {
    if (size > max_size_)
    {
      free_data();
      cuda_malloc(data_, size);
      max_size_ = size;
      size_ = size;
      allocated = true;
    }
    else
    {
      size_ = size;
    }
  }
  else
  {
    cuda_malloc(data_, size);
    size_ = size;
    max_size_ = size;
    allocated = true;
  }
}

template<typename T>
__host__
void dvec<T>::assign(T* data_h, int size, cudaStream_t *stream)
{
  resize(size);

  if (stream ==  NULL)
    cuda_copy_h2d(data_, data_h, size);
  else
    cuda_copy_h2d(data_, data_h, size, *stream);
}

template<typename T>
__host__
void dvec<T>::free_data(void)
{
  if (allocated)
  {
    cuda_free(data_);
    size_ = max_size_ = 0;
  }
}

template<typename T>
class hvec
{
private:
  int size_ = 0;
  int max_size_ = 0;
  T *data_ = NULL;
  bool allocated = false;

public:
  hvec(void) { }

  ~hvec(void);

  int size(void) const { return size_; }

  int capacity(void) const { return max_size_; }

  T* data(void) { return data_; }

  void resize(int size);

  void assign(T* data_d, int size, cudaStream_t *stream = NULL);

  void free_data(void);

  T operator[](int ind) const;
  T& operator[](int ind);
};

template<typename T>
hvec<T>::~hvec(void)
{
  free_data();
}

template<typename T>
T hvec<T>::operator[](int ind) const
{
  return data_[ind];
}

template<typename T>
T& hvec<T>::operator[](int ind)
{
  return data_[ind];
}

template<typename T>
void hvec<T>::resize(int size)
{
  if (allocated)
  {
    if (size > max_size_)
    {
      free_data();
      cuda_malloc_pinned(data_, size);
      max_size_ = size;
      size_ = size;
      allocated = true;
    }
    else
    {
      size_ = size;
    }
  }
  else
  {
    cuda_malloc_pinned(data_, size);
    size_ = size;
    max_size_ = size;
    allocated = true;
  }
}

template<typename T>
void hvec<T>::assign(T* data_d, int size, cudaStream_t *stream)
{
  resize(size);

  if (stream ==  NULL)
    cuda_copy_d2h(data_d, data_, size);
  else
    cuda_copy_d2h(data_d, data_, size, *stream);
}

template<typename T>
void hvec<T>::free_data(void)
{
  if (allocated)
  {
    cuda_free_pinned(data_);
    size_ = max_size_ = 0;
  }

  allocated = false;
}

struct dPoint
{
private:
  double _v[3] = {0.,0.,0.};

public:
  __device__ dPoint(void) {}

  __device__
  dPoint(double x, double y, double z) {
    _v[0] = x; _v[1] = y; _v[2] = z;
  }

  __device__
  dPoint(const double* pt) {
    _v[0] = pt[0]; _v[1] = pt[1]; _v[2] = pt[2];
  }

  __device__
  double& operator[](int ind) {
    return _v[ind];
  }

  __device__
  dPoint operator=(double* a) {
    struct dPoint pt;
    for (int i = 0; i < 3; i++)
      pt[i] = a[i];
    return pt;
  }

  __device__
  dPoint operator-(dPoint b) {
    struct dPoint c;
    for (int i = 0; i < 3; i++)
      c[i] = _v[i] - b[i];
    return c;
  }

  __device__
  dPoint operator+(dPoint b) {
    struct dPoint c;
    for (int i = 0; i < 3; i++)
      c[i] = _v[i] + b[i];
    return c;
  }

  __device__
  dPoint operator/(dPoint b) {
    struct dPoint c;
    for (int i = 0; i < 3; i++)
      c[i] = _v[i] / b[i];
    return c;
  }

  __device__
  dPoint& operator+=(dPoint b) {
    for (int i = 0; i < 3; i++)
      _v[i] += b[i];
    return *this;
  }

  __device__
  dPoint& operator-=(dPoint b) {
    for (int i = 0; i < 3; i++)
      _v[i] -= b[i];
    return *this;
  }

  __device__
  dPoint& operator+=(double* b) {
    for (int i = 0; i < 3; i++)
      _v[i] += b[i];
    return *this;
  }

  __device__
  dPoint& operator-=(double* b) {
    for (int i = 0; i < 3; i++)
      _v[i] -= b[i];
    return *this;
  }

  __device__
  dPoint& operator/=(double a) {
    for (int i = 0; i < 3; i++)
      _v[i] /= a;
    return *this;
  }

  __device__
  dPoint& operator*=(double a) {
    for (int i = 0; i < 3; i++)
      _v[i] *= a;
    return *this;
  }

  __device__
  double operator*(dPoint b) {
    return _v[0]*b[0] + _v[1]*b[1] + _v[2]*b[2];
  }

  __device__
  dPoint operator*(double b) {
    struct dPoint c;
    for (int i = 0; i < 3; i++)
      c[i] = _v[i] * b;
    return c;
  }

  __device__
  dPoint operator/(double b) {
    struct dPoint c;
    for (int i = 0; i < 3; i++)
      c[i] = _v[i] / b;
    return c;
  }

  __device__
  void abs(void) {
    for (int i = 0; i < 3; i++)
      _v[i] = fabs(_v[i]);
  }

  __device__
  double norm(void) {
    return sqrt(_v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]);
  }

  __device__
  dPoint cross(dPoint b) {
    dPoint c;
    c[0] = _v[1]*b[2] - _v[2]*b[1];
    c[1] = _v[2]*b[0] - _v[0]*b[2];
    c[2] = _v[0]*b[1] - _v[1]*b[0];
    return c;
  }
};

static
__device__
dPoint operator*(double a, dPoint b)
{
  dPoint c;
  for (int i = 0; i < 3; i++)
    c[i] = a * b[i];
  return c;
}

/* ------ Misc. Helper Functions ------ */

namespace cuda_funcs
{

static
__device__ __forceinline__
double det_3x3_part(const double* mat, int a, int b, int c)
{
  return mat[a] * (mat[3+b] * mat[6+c] - mat[3+c] * mat[6+b]);
}

static
__device__ __forceinline__
double det_3x3(const double* mat)
{
  return det_3x3_part(mat,0,1,2) - det_3x3_part(mat,1,0,2)
      + det_3x3_part(mat,2,0,1);
}

static
__device__ __forceinline__
void adjoint_3x3(const double* __restrict__ mat, double* __restrict__ adj)
{
  double a11 = mat[0], a12 = mat[1], a13 = mat[2];
  double a21 = mat[3], a22 = mat[4], a23 = mat[5];
  double a31 = mat[6], a32 = mat[7], a33 = mat[8];

  adj[0] = a22*a33 - a23*a32;
  adj[1] = a13*a32 - a12*a33;
  adj[2] = a12*a23 - a13*a22;

  adj[3] = a23*a31 - a21*a33;
  adj[4] = a11*a33 - a13*a31;
  adj[5] = a13*a21 - a11*a23;

  adj[6] = a21*a32 - a22*a31;
  adj[7] = a12*a31 - a11*a32;
  adj[8] = a11*a22 - a12*a21;
}

/*! Evaluates the Lagrange function corresponding to the specified mode on xiGrid at location xi.
 *
 * \param xiGrid The grid of interpolation points. Sorted in domain [-1,1].
 * \param mode Mode of the Lagrange function. Defined such that function is 1 at xiGrid(mode)
 * zero at other grid points.
 * \param xi  Point of evaluation in domain [-1,1].
 *
 * \return Value of Lagrange function at xi.
 */
__device__ __forceinline__
double Lagrange_gpu(const double* __restrict__ xiGrid, unsigned int npts,
                    double xi, unsigned int mode)
{
  double val = 1.0;

  for (unsigned int i = 0; i < mode; i++)
    val *= (xi - xiGrid[i])/(xiGrid[mode] - xiGrid[i]);

  for (unsigned int i = mode + 1; i < npts; i++)
    val *= (xi - xiGrid[i])/(xiGrid[mode] - xiGrid[i]);

  return val;
}

/*! Evaluates the first derivative of the Lagrange function corresponding to the specified mode
 *  on xiGrid at location xi.
 *
 * \param xiGrid The grid of interpolation points. Sorted in domain [-1,1].
 * \param mode Mode of the Lagrange function. Defined such that function is 1 at xiGrid[mode]
 * zero at other grid points.
 * \param xi  Point of evaluation in domain [-1,1].
 *
 * \return Value of first derivative of the Lagrange function at xi.
 */
__device__ __forceinline__
double dLagrange_gpu(const double* __restrict__ xiGrid, unsigned int npts,
                     double xi, unsigned int mode)
{
  double val = 0.0;

  /* Compute normalization constant */
  double den = 1.0;
  for (unsigned int i = 0; i < mode; i++)
    den *= (xiGrid[mode] - xiGrid[i]);

  for (unsigned int i = mode+1; i < npts; i++)
    den *= (xiGrid[mode] - xiGrid[i]);

  /* Compute sum of products */
  for (unsigned int j = 0; j < npts; j++)
  {
    if (j == mode)
      continue;

    double term = 1.0;
    for (unsigned int i = 0; i < npts; i++)
      if (i != mode and i != j)
        term *= (xi - xiGrid[i]);

    val += term;
  }

  return val/den;
}

__device__ __forceinline__
char checkHoleMap(const double* __restrict__ pt, const char* __restrict__ sam,
    const int* __restrict__ nx, const double* __restrict__ extents)
{
  double dx[3];
  int ix[3];

  for (int i = 0; i < 3; i++)
    dx[i] = (extents[i+3]-extents[i]) / nx[i];

  for (int i = 0; i < 3; i++)
  {
    ix[i] = (pt[i]-extents[i])/dx[i];
    if (ix[i] < 0 || ix[i] > nx[i]-1) return 0;
  }

  return sam[ix[0] + nx[0]*(ix[1] + nx[1]*ix[2])];
}


__device__ __forceinline__
char checkHoleMap(const float* __restrict__ pt, const char* __restrict__ sam,
    const int* __restrict__ nx, const double* __restrict__ extents)
{
  double dx[3];
  int ix[3];

  for (int i = 0; i < 3; i++)
    dx[i] = (extents[i+3]-extents[i]) / nx[i];

  for (int i = 0; i < 3; i++)
  {
    ix[i] = (pt[i]-extents[i])/dx[i];
    if (ix[i] < 0 || ix[i] > nx[i]-1) return 0;
  }

  return sam[ix[0] + nx[0]*(ix[1] + nx[1]*ix[2])];
}

template<int nDims, int nPts>
__device__ __forceinline__
void getBoundingBox(const float* __restrict__ pts, float *bbox)
{
  for (int i = 0; i < nDims; i++)
  {
    bbox[i]       =  BIG_DOUBLE;
    bbox[nDims+i] = -BIG_DOUBLE;
  }

  for (int i = 0; i < nPts; i++)
  {
    for (int dim = 0; dim < nDims; dim++)
    {
      bbox[dim]       = fminf(bbox[dim],      pts[i*nDims+dim]);
      bbox[nDims+dim] = fmaxf(bbox[nDims+dim],pts[i*nDims+dim]);
    }
  }
}

template<int nDims, int nPts>
__device__ __forceinline__
void getBoundingBox(const double* __restrict__ pts, double *bbox)
{
  for (int i = 0; i < nDims; i++)
  {
    bbox[i]       =  INFINITY;
    bbox[nDims+i] = -INFINITY;
  }

  for (int i = 0; i < nPts; i++)
  {
    for (int dim = 0; dim < nDims; dim++)
    {
      bbox[dim]       = fmin(bbox[dim],      pts[i*nDims+dim]);
      bbox[nDims+dim] = fmax(bbox[nDims+dim],pts[i*nDims+dim]);
    }
  }
}

template<int nDims>
__device__ __forceinline__
bool boundingBoxCheck(const float* __restrict__ bbox1,
                      const float* __restrict__ bbox2, float tol)
{
  bool check = true;
  for (int i = 0; i < nDims; i++)
  {
    check = check && (bbox1[i+nDims] >= bbox2[i] - tol);
    check = check && (bbox2[i+nDims] >= bbox1[i] - tol);
  }
  return check;
}

template<int nDims>
__device__ __forceinline__
bool boundingBoxCheck(const double* __restrict__ bbox1,
                      const double* __restrict__ bbox2, double tol)
{
  bool check = true;
  for (int i = 0; i < nDims; i++)
  {
    check = check && (bbox1[i+nDims] >= bbox2[i] - tol);
    check = check && (bbox2[i+nDims] >= bbox1[i] - tol);
  }
  return check;
}

template<int nDims>
__device__ __forceinline__
float boundingBoxDist(const float* __restrict__ bbox1,
                      const float* __restrict__ bbox2)
{
  float dist = 0.f;
  for (int i = 0; i < nDims; i++)
  {
    if ( (bbox1[i+nDims] < bbox2[i]) || (bbox2[i+nDims] < bbox1[i]) )
    {
      float d = fmaxf(bbox2[i]-bbox1[i+nDims], bbox1[i]-bbox2[i+nDims]);
      dist += d*d;
    }
  }

  return sqrt(dist);
}

template<int nDims>
__device__ __forceinline__
double boundingBoxDist(const double* __restrict__ bbox1,
                       const double* __restrict__ bbox2)
{
  double dist = 0.;
  for (int i = 0; i < nDims; i++)
  {
    if ( (bbox1[i+nDims] < bbox2[i]) || (bbox2[i+nDims] < bbox1[i]) )
    {
      double d = fmax(bbox2[i]-bbox1[i+nDims], bbox1[i]-bbox2[i+nDims]);
      dist += d*d;
    }
  }

  return sqrt(dist);
}

template<int nDims, int nPts>
__device__
void getCentroid(const float* __restrict__ pts, float* __restrict__ xc)
{
  for (int d = 0; d < nDims; d++)
    xc[d] = 0.f;

  for (int i = 0; i < nPts; i++)
    for (int d = 0; d < nDims; d++)
      xc[d] += pts[nDims*i+d]/nPts;
}

template<int nDims, int nPts>
__device__
void getCentroid(const double* __restrict__ pts, double* __restrict__ xc)
{
  for (int d = 0; d < nDims; d++)
    xc[d] = 0;

  for (int i = 0; i < nPts; i++)
    for (int d = 0; d < nDims; d++)
      xc[d] += pts[nDims*i+d]/nPts;
}

/*! Generate an approximate Object-Oriented Bounding Box for a hexahedron 
 *  OOBB stored as body axes, and min/max points of box in body frame (15 floats) */
template<int nDims, int nPts>
__device__
void getOOBB(const double* __restrict__ pts, float* __restrict__ oobb, bool PRINT)
{
  // List of edges in 8-node hex: xi-edges, eta-edges, zeta-edges
  const char edges[12][2] = { {0,1}, {3,2}, {4,5}, {7,6}, {0,3}, {1,2}, {4,7}, {5,6}, {0,4}, {1,5}, {3,7}, {2,6} };

  // 1) Find the longest edge in each quasi-reference direction
  int maxE[3];
  //char maxE[3];
  float length[3] = {0.f,0.f,0.f}; /// TODO: Use floats in more places throughout Direct Cut
  for (int j = 0; j < 3; j++)
  {
    for (int i = 0; i < 4; i++)
    {
      float dist = 0.f;
      int P0 = (int)edges[4*j+i][0];
      int P1 = (int)edges[4*j+i][1];
      for (int d = 0; d < nDims; d++)
        dist += (pts[3*P1+d] - pts[3*P0+d]) * (pts[3*P1+d] - pts[3*P0+d]);

      if (dist > length[j])
      {
        length[j] = dist;
        //maxE[j] = (char)i;
        maxE[j] = i;
      }
    }
  }

  // 1.5) Sort directions by length [descending]
  int Dims[3] = {0,1,2};
  if (length[1] > length[0])
  {
    swap_d(length[1], length[0]);
    swap_d(maxE[1], maxE[0]);
    swap_d(Dims[1], Dims[0]);
  }

  if (length[2] > length[1])
  {
    swap_d(length[1], length[2]);
    swap_d(maxE[1], maxE[2]);
    swap_d(Dims[1], Dims[2]);

    if (length[1] > length[0])
    {
      swap_d(length[1], length[0]);
      swap_d(maxE[1], maxE[0]);
      swap_d(Dims[1], Dims[0]);
    }
  }
   
  for (int d = 0; d < 3; d++)
    length[d] = sqrt(length[d]);

  // 2) Setup new coordinate axes
  float R[9];
  for (int m = 0; m < 3; m++)
  {
    int D = (int)Dims[m];
    int I = (int)maxE[m];
    int P0 = edges[4*D+I][0];
    int P1 = edges[4*D+I][1];
    for (int d = 0; d < 3; d++)
      R[3*m+d] = (pts[3*P1+d] - pts[3*P0+d]) / length[D];
  }

  // 2.5) Orthoganlize & normalize our approxmiate basis
  
  // Direction 2
  float dot = 0.f;
  for (int d = 0; d < 3; d++)
    dot += R[d]*R[3+d];
  for (int d = 0; d < 3; d++)
    R[3+d] -= dot*R[d];

  dot = 0.f;
  for (int d = 0; d < 3; d++)
    dot += R[3+d]*R[3+d];
  dot = sqrt(dot);
  for (int d = 0; d < 3; d++)
    R[3+d] /= dot;


  // Direction 3
  dot = 0.f;
  for (int d = 0; d < 3; d++)
    dot += R[d]*R[6+d];
  for (int d = 0; d < 3; d++)
    R[6+d] -= dot*R[d];

  dot = 0.f;
  for (int d = 0; d < 3; d++)
    dot += R[3+d]*R[6+d];
  for (int d = 0; d < 3; d++)
    R[6+d] -= dot*R[3+d];

  dot = 0.f;
  for (int d = 0; d < 3; d++)
    dot += R[6+d]*R[6+d];
  dot = sqrt(dot);
  for (int d = 0; d < 3; d++)
    R[6+d] /= dot;

  // 3) Find min/max pts in new axes
  float minPt[3] = { BIG_DOUBLE,  BIG_DOUBLE,  BIG_DOUBLE};
  float maxPt[3] = {-BIG_DOUBLE, -BIG_DOUBLE, -BIG_DOUBLE};
  for (int i = 0; i < nPts; i++)
  {
    float pt[3] = {0.f, 0.f, 0.f};
    for (int m = 0; m < 3; m++)
      for (int n = 0; n < 3; n++)
        pt[m] += R[3*m+n] * pts[3*i+n];

    for (int d = 0; d < 3; d++)
    {
      minPt[d] = fminf(minPt[d], pt[d]);
      maxPt[d] = fmaxf(maxPt[d], pt[d]);
    }
  }

  // 4) Store the OOBB
  for (int i = 0; i < 9; i++)
    oobb[i] = R[i];
  for (int i = 0; i < 3; i++)
    oobb[9+i] = minPt[i];
  for (int i = 0; i < 3; i++)
    oobb[12+i] = maxPt[i];
}

} // namespace cuda_funcs

#endif
