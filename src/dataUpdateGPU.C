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

void MeshBlock::getInterpolatedSolutionGPU(int *nints,int *nreals,
					   int **intData,double **realData,double *q,
					   int nvar, int interptype)
{
  int i;
  int k,m,inode;
  double weight;
  double *qq;
  int icount,dcount;

  //
  qq=(double *)malloc(sizeof(double)*nvar);
  //
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
  //
  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
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
}
	
  
#endif
	   
    
      
