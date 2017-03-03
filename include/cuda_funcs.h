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
  dvec(void) { }

  ~dvec(void);

  int size(void) { return size_; }

  int capacity(void) { return max_size_; }

  T* data(void) { return data_; }

  void resize(int size);

  void assign(T* data_h, int size, cudaStream_t *stream = NULL);

  void free_data(void);
};

template<typename T>
dvec<T>::~dvec(void)
{
  free_data();
}

template<typename T>
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
void dvec<T>::assign(T* data_h, int size, cudaStream_t *stream)
{
  resize(size);

  if (stream ==  NULL)
    cuda_copy_h2d(data_, data_h, size);
  else
    cuda_copy_h2d(data_, data_h, size, *stream);
}

template<typename T>
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

#endif
