#include <stdio.h>
#include "cuda_functions.h"

__global__ void vec_mult(double* v1, double* v2, double *result, int pts){
  
  int i   = blockDim.x * blockIdx.x + threadIdx.x;

  if( i < pts ){
    result[i] = v1[i]*v2[i];
  }

}

void tmp_update(int pts, double* q){

  printf("This is the tmp_update function from update.cu\n");
}

void freeGPUInterpList(INTERPLIST* d_interpList){

  INTERPLIST first;

  // we allocated all the memory for inode and weights in the first
  // entry of the interplist on the GPU. All other entries in the
  // interpList just point to offsets from that memory.
  cudaMemcpy(&first, d_interpList, sizeof(INTERPLIST), cudaMemcpyDeviceToHost);
  cudaFree(first.inode);
  cudaFree(first.weights);
  cudaFree(d_interpList);

}

void allocGPUInterpList(INTERPLIST** d_interplist, int ninterp, INTERPLIST* interplist){

  int *d_inode, *h_inode;
  double *d_weights, *h_weights;
  INTERPLIST* tmplist;
  int count;

  // copy the list on the cpu into a temporary location
  tmplist = (INTERPLIST*)malloc(ninterp*sizeof(INTERPLIST));
  memcpy(tmplist, interplist, ninterp*sizeof(INTERPLIST));

  count = 0;
  for(int i=0; i<ninterp; i++){
    count+=tmplist[i].nweights;
  }

  // make space on the GPU
  HANDLE_ERROR( cudaMalloc((void**)d_interplist, ninterp*sizeof(INTERPLIST)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_inode,   count*sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_weights, count*sizeof(double)) );
  h_inode   = new int[count];
  h_weights = new double[count];

  count = 0;
  for(int i=0; i<ninterp; i++){
    tmplist[i].inode   = &d_inode[count];
    tmplist[i].weights = &d_weights[count];
    for(int k=0; k<interplist[i].nweights; k++){
      h_inode[count+k]   = interplist[i].inode[k];
      h_weights[count+k] = interplist[i].weights[k];
    }
    count += tmplist[i].nweights;
  }

  HANDLE_ERROR( cudaMemcpy(d_inode,   h_inode,   count*sizeof(int), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_weights, h_weights, count*sizeof(double), cudaMemcpyHostToDevice) );

  // now copy the whole list to the GPU:
  HANDLE_ERROR( cudaMemcpy((*d_interplist), tmplist, ninterp*sizeof(INTERPLIST),
			   cudaMemcpyHostToDevice) );

  // and cleanup
  free(tmplist);
  delete[] h_inode;
  delete[] h_weights;
}

__global__ void interp_vec(double* q, int nvar, int ninterp, INTERPLIST* interplist, 
		      int* intData, double* realData, int* cntr){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int idx, m, k, inode;
  double weight;
  if(i < ninterp && !interplist[i].cancel){

    // Atomic operations on global memory are faster now with CUDA
    // 9. No need to do fancy warp stuff:
    // https://devblogs.nvidia.com/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/
    // atomicAdd returns the old number and increments the pointer contents
    idx = atomicAdd(cntr, 1); 

    m = 0;
    // assignment, not +=
    inode=interplist[i].inode[m];
    weight=interplist[i].weights[m];
    for(k=0;k<nvar;k++){
      realData[idx*nvar+k]=q[inode*nvar+k]*weight; 
    }
    // now we use the += operator
    for(m=1;m<interplist[i].nweights;m++){
      inode=interplist[i].inode[m];
      weight=interplist[i].weights[m];
      for(k=0;k<nvar;k++){
      	realData[idx*nvar+k]+=q[inode*nvar+k]*weight;
      }
    }

    intData[idx*2+0]=interplist[i].receptorInfo[0];
    intData[idx*2+1]=interplist[i].receptorInfo[1];

  } // end if

}

__global__ void update_vec(double* q, int nvar, int nupdate, int* idata, double* ddata){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int m;

  if(i<nupdate){
    for(m=0; m<nvar; m++){
      // printf("_g_ q[%d] = q[%d]\n", idata[i]*nvar+m, i*nvar+m);
      q[idata[i]*nvar+m] = ddata[i*nvar+m];
    }
  }

}

__global__ void idbg_kernel(int* vec, int n){
  for(int i=0; i<n; i++){
    printf("_dbggpu__ %2d : %7d\n", vec[i]);
  }
}
__global__ void ddbg_kernel(double* vec, int n){
  for(int i=0; i<n; i++){
    printf("_dbggpu__ %2d : %16.8e\n", vec[i]);
  }
}

void interpolateVectorGPU(GPUvec<double> *vec, int nints, int nreals, int ninterp,
			  int** intData, double** realData, INTERPLIST* interpList){

  int* cntr;
  dim3 blocks(1,1,1), threads(256,1,1);
  blocks.x = (ninterp-1)/threads.x+1;

  int *d_intData;
  double *d_realData;
  HANDLE_ERROR( cudaMalloc((void**)&d_intData, 2*nints*sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_realData, nreals*sizeof(double)) );
  HANDLE_ERROR( cudaMalloc((void**)&cntr, sizeof(int)) );
  HANDLE_ERROR( cudaMemset(cntr, 0, sizeof(int)) );

  interp_vec<<<blocks,threads>>>(vec->data,vec->nvar,ninterp,interpList,d_intData,d_realData,cntr);

  // Allocate on CPU
  if((*intData) == NULL){
    (*intData)=(int *)malloc(sizeof(int)*2*nints);
    (*realData)=(double *)malloc(sizeof(double)*nreals);
  } else {
    printf("skipping alloc, assuming ok\n");
  }

  // Transfer to CPU
  HANDLE_ERROR( cudaMemcpy( (*intData),  d_intData, 2*nints*sizeof(int)   , cudaMemcpyDeviceToHost) );
  HANDLE_ERROR( cudaMemcpy( (*realData), d_realData, nreals*sizeof(double), cudaMemcpyDeviceToHost) );

  HANDLE_ERROR( cudaFree(d_intData) );
  HANDLE_ERROR( cudaFree(d_realData) );
  HANDLE_ERROR( cudaFree(cntr) );
}


void updateSolnGPU(int nrecv, PACKET *rcvPack, GPUvec<double> *vec){

  int i, k, m, count;
  int *h_idata, *d_idata;
  double *h_ddata, *d_ddata;
  dim3 blocks(1,1,1), threads(256,1,1);
  int nvar = vec->nvar;

  count = 0;
  for(k=0; k<nrecv; k++){
    count += rcvPack[k].nints;
  }
  if(count == 0) return;

  blocks.x = (count-1)/threads.x+1;

  h_idata = new int[count];
  h_ddata = new double[count*nvar];

  HANDLE_ERROR( cudaMalloc((void**)&d_idata, count*sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_ddata, count*nvar*sizeof(double)) );

  count = 0;
  for(k=0; k<nrecv; k++){
    for(i=0; i<rcvPack[k].nints; i++){
      h_idata[count+i] = rcvPack[k].intData[i];
      for(m=0; m<nvar; m++){
  	h_ddata[(count+i)*nvar+m] = rcvPack[k].realData[i*nvar+m];
      }
    }
    count += rcvPack[k].nints;
  }

  HANDLE_ERROR( cudaMemcpy(d_idata, h_idata, count*sizeof(int), 
  			   cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_ddata, h_ddata, count*nvar*sizeof(double), 
  			   cudaMemcpyHostToDevice) );

  update_vec<<<blocks,threads>>>(vec->data, nvar, count, d_idata, d_ddata);

  HANDLE_ERROR( cudaFree(d_idata) );
  HANDLE_ERROR( cudaFree(d_ddata) );
  delete[] h_idata;
  delete[] h_ddata;

}
