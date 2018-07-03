//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#ifdef USE_CUDA
#include "MeshBlock.h"
#include <string.h>
#include "cuda_functions.h"

void MeshBlock::getInterpolatedSolution(int *nints,int *nreals,int **intData,double **realData,
					GPUvec<double> *vec)
{
  int i;
  int k,m,inode;
  double weight;
  int icount,dcount;
  int nvar  = vec->nvar;

  // If this is the first time doing a dataUpdate since connecting, we
  // need to allocate the interpList on the GPU.
  if(d_interpList == NULL){
    allocGPUInterpList(&d_interpList, ninterp, interpList);
  }

  (*nints)=(*nreals)=0;
  for(i=0;i<ninterp;i++)
    {
      if (!interpList[i].cancel)
  	{
  	  (*nints)++;
  	  (*nreals)=(*nreals)+nvar;
  	}
    }
  if ((*nints)==0) return;

  interpolateVectorGPU(vec, (*nints), (*nreals), ninterp, intData, realData, d_interpList);

}

#endif
	   
    
      
