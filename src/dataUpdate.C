// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"
#include "MeshBlock.h"
#include <string.h>

#define ROW 0
#define COLUMN 1

void MeshBlock::getInterpolatedSolution(int *nints,int *nreals,int **intData,double **realData,double *q,
					int nvar, int interptype)
{
  int i;
  int k,m,inode;
  double weight;
  double *qq = NULL;
  int icount,dcount;
  //
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
  qq=(double *)malloc(sizeof(double)*nvar);
  (*intData)=(int *)malloc(sizeof(int)*3*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  icount=dcount=0;
  //
  if (interptype==ROW)
    {    
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
                    TRACED(weight);
                    printf("warning: weights are not convex 1\n");
                   }
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[inode*nvar+k]*weight;
		}
	      (*intData)[icount++]=interpList[i].receptorInfo[0];
	      (*intData)[icount++]=interpList[i].receptorInfo[1];
	      (*intData)[icount++]=interpList[i].receptorInfo[2];
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
    }
  else if (interptype==COLUMN)
    {
      for(i=0;i<ninterp;i++)
	{
	  if (!interpList[i].cancel)
	    {
	      for(k=0;k<nvar;k++) qq[k]=0;
	      for(m=0;m<interpList[i].nweights;m++)
		{
		  inode=interpList[i].inode[m];
		  weight=interpList[i].weights[m];
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[k*nnodes+inode]*weight;
		}
	      (*intData)[icount++]=interpList[i].receptorInfo[0];
	      (*intData)[icount++]=interpList[i].receptorInfo[1];
	      (*intData)[icount++]=interpList[i].receptorInfo[2];
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
    }

  if (qq) TIOGA_FREE(qq);
}
	
void MeshBlock::updateSolnData(int inode,double *qvar,double *q,int nvar,int interptype)
{
  int k;

  if (interptype==ROW)
    {
      for(k=0;k<nvar;k++)
	q[inode*nvar+k]=qvar[k];
    }
  if (interptype==COLUMN)
    {
      for(k=0;k<nvar;k++)
	q[nnodes*k+inode]=qvar[k];
    }
}

void MeshBlock::getDonorCount(int *dcount,int *fcount)
{
  int i;  
  *dcount=0;
  *fcount=0;
  for(i=0;i<ninterp;i++)
    {
      if (!interpList[i].cancel) 
	{
	  (*dcount)++;
	  (*fcount)+=(interpList[i].nweights+1);
	}
    }
}

void MeshBlock::getDonorInfo(int *receptors,int *indices,double *frac)
{
  int i,j,k,m;
  int dcount=0;

  j=0;
  k=0;
  for(i=0;i<ninterp;i++)
    {
      if (!interpList[i].cancel) 
	{
	  for(m=0;m<interpList[i].nweights+1;m++)
	    {
	      indices[j]=interpList[i].inode[m];
	      frac[j]=interpList[i].weights[m];
	      j++;
	    }
	  receptors[k++]=interpList[i].receptorInfo[0];
	  receptors[k++]=interpList[i].receptorInfo[1];
	  receptors[k++]=interpList[i].receptorInfo[2];
	  receptors[k++]=interpList[i].nweights;
	}
    }
}

void MeshBlock::getReceptorInfo(int *receptors)
{
  int k=0;
  for (int i=0; i<ninterp; i++) {
    if (interpList[i].cancel) continue;

    receptors[k++] = interpList[i].receptorInfo[0];
    receptors[k++] = interpList[i].receptorInfo[1];
    receptors[k++] = interpList[i].receptorInfo[2];

    int donID = interpList[i].inode[interpList[i].nweights];

    // Copy the contents of uint64_t (8 bytes) into 2 4-byte locations in the
    // array
    memcpy(&receptors[k], &cellGID[donID], sizeof(uint64_t));
    k += 2;
  }
}
