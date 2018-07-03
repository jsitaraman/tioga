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

  int *d_inode;
  double *d_weights;
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

  // the first item in the list points to the full array
  tmplist[0].inode   = d_inode;
  tmplist[0].weights = d_weights;
  // each item after that points to an offset of the array
  count = tmplist[0].nweights;
  for(int i=1; i<ninterp; i++){
    tmplist[i].inode   = &d_inode[count];
    tmplist[i].weights = &d_weights[count];
    count += tmplist[i].nweights;
  }

  // now copy the whole list to the GPU:
  HANDLE_ERROR( cudaMemcpy((*d_interplist), tmplist, ninterp*sizeof(INTERPLIST),
			   cudaMemcpyHostToDevice) );

  // and cleanup
  free(tmplist);
}

//__global__ interp_vec(double* q, int nvar, INTERPLIST* interplist

void interpolateVectorGPU(GPUvec<double> *vec, int nints, int nreals, int ninterp,
			  int** intData, double** realData, INTERPLIST* interpList){

  int i;
  int k,m,inode;
  double weight;
  double *qq, *q;
  int icount,dcount;
  int nvar = vec->nvar;

  // int *d_intData;
  // double *d_realData;
  // HANDLE_ERROR( cudaMalloc((void**)&d_intData, 2*nints*sizeof(int)) );
  // HANDLE_ERROR( cudaMalloc((void**)&d_realData, nreals*sizeof(double)) );

  q  = new double[vec->nvar*vec->pts];
  qq = new double[nvar];
  vec->to_cpu(q);

  //
  (*intData)=(int *)malloc(sizeof(int)*2*nints);
  (*realData)=(double *)malloc(sizeof(double)*nreals);
  icount=dcount=0;
  //
  for(i=0;i<ninterp;i++)
    {
      if (!interpList[i].cancel)
	{
	  for(k=0;k<nvar;k++) qq[k]=0;
	  for(m=0;m<interpList[i].nweights;m++)
	    {
	      inode=interpList[i].inode[m];
	      weight=interpList[i].weights[m];
	      if (weight < -TOL || weight > 1.0+TOL) {
		traced(weight);
		printf("warning: weights are not convex 1\n");
	      }
	      for(k=0;k<nvar;k++)
		qq[k]+=q[inode*nvar+k]*weight;
	    }
	  (*intData)[icount++]=interpList[i].receptorInfo[0];
	  (*intData)[icount++]=interpList[i].receptorInfo[1];
	  for(k=0;k<nvar;k++)
	    (*realData)[dcount++]=qq[k];
	}
    }

  delete qq;
  delete q;

  // HANDLE_ERROR( cudaFree(d_intData) );
  // HANDLE_ERROR( cudaFree(d_realData) );
}
