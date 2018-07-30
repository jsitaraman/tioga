#include "cuda_functions.h"

#define REAL_GPU

template<typename T>
GPUvec<T>::GPUvec(int p, int nv){
  this->owndata = true;
  this->pts     = p;
  this->nvar    = nv;
#ifdef REAL_GPU
  HANDLE_ERROR( cudaMalloc((void**) &this->data, nv*p*sizeof(T)) );
#else
  this->data = new T[p*nv];
#endif
}

template<typename T>
GPUvec<T>::GPUvec(int p, int nv, T* d){
  this->owndata = false;
  this->data    = d;
  this->pts     = p;
  this->nvar    = nv;
}

template<typename T>
GPUvec<T>::~GPUvec(){
  if(this->owndata){
    // printf("deallocating gpuvec\n");
#ifdef REAL_GPU
    HANDLE_ERROR( cudaFree(this->data) );
#else
    delete[] this->data;
#endif
  }
  this->data = NULL;
}

template<typename T>
void GPUvec<T>::to_gpu(T* q){
#ifdef REAL_GPU
  HANDLE_ERROR( cudaMemcpy(data, q, nvar*pts*sizeof(T), 
			   cudaMemcpyHostToDevice) );
#else
  memcpy(data, q, nvar*pts*sizeof(T));
#endif
}

template<typename T>
void GPUvec<T>::to_cpu(T* q){
#ifdef REAL_GPU
  HANDLE_ERROR( cudaMemcpy(q, data, nvar*pts*sizeof(T),
			   cudaMemcpyDeviceToHost) );
#else
  memcpy(q, data, nvar*pts*sizeof(T));
#endif
}

// explicit class type declaration
template class GPUvec<double>;
template class GPUvec<int>;
