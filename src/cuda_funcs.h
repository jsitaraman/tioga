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
  check_error();
}

template <typename T>
void cuda_free_pinned(T* &data_h)
{
  cudaFreeHost(data_h);
  check_error();
}
