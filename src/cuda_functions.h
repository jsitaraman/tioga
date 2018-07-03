#ifndef TIOGA_CUDA_FUNCTIONS_H
#define TIOGA_CUDA_FUNCTIONS_H
#include <stdio.h>
#include "codetypes.h"

#ifdef __NVCC__
//#ifdef __CUDACC__
static void HandleError(cudaError_t err, const char *file, int line){
  if(err != cudaSuccess){
    printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__))
#endif


template<typename T>
class GPUvec{
  bool owndata;
public:
  GPUvec(int p, int nv);
  GPUvec(int p, int nv, T* d);
  ~GPUvec(void);
  void to_gpu(T* q);
  void to_cpu(T* q);
  T* data;
  int pts;
  int nvar;
};

void tmp_update(int pts, double* q);

void freeGPUInterpList(INTERPLIST* d_interplist);
void allocGPUInterpList(INTERPLIST** d_interplist, int ninterp, INTERPLIST* interplist);


void interpolateVectorGPU(GPUvec<double> *vec, int nints, int nreals, int ninterp,
			  int** intData, double** realData, INTERPLIST* interplist);

void updateSolnGPU(int nrecv, PACKET *rcvPack, GPUvec<double> *vec);

#endif
