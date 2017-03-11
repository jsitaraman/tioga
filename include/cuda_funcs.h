#ifndef CUDA_FUNCS_H
#define CUDA_FUNCS_H

#include "cuda_runtime.h"
#include "error.hpp"

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
  int size(void) { return size_; }

  __host__ __device__
  int capacity(void) { return max_size_; }

  __host__ __device__
  T* data(void) { return data_; }

  __host__
  void resize(int size);

  __host__
  void assign(T* data_h, int size, cudaStream_t *stream = NULL);

  __host__
  void free_data(void);

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
__device__
T& dvec<T>::operator [](int ind)
{
  return data_[ind];
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

  int size(void) { return size_; }

  int capacity(void) { return max_size_; }

  T* data(void) { return data_; }

  void resize(int size);

  void assign(T* data_d, int size, cudaStream_t *stream = NULL);

  void free_data(void);

  T& operator[](int ind) { return data_[ind]; }
};

template<typename T>
hvec<T>::~hvec(void)
{
  free_data();
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
double Lagrange_gpu(double* xiGrid, unsigned int npts, double xi, unsigned int mode)
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
double dLagrange_gpu(double* xiGrid, unsigned int npts, double xi, unsigned int mode)
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

} // namespace cuda_funcs

#endif
