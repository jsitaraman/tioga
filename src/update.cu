#include <stdio.h>

__global__ void vec_mult(double* v1, double* v2, double *result, int pts){
  
  int i   = blockDim.x * blockIdx.x + threadIdx.x;

  if( i < pts ){
    result[i] = v1[i]*v2[i];
  }

}

void tmp_update(int pts, double* q){

  printf("This is the tmp_update function from update.cu\n");
}
