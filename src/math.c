// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "codetypes.h"

void solvec(double **a,double *b,int *iflag,int n)
{
  int i,j,k,l,flag,temp1;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;

  
  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  temp=b[k];
		  b[k]=b[i];
		  b[i]=temp;
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      b[k]+=fact*b[i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      sum=0;
      for(j=i+1;j<n;j++)
	sum+=a[i][j]*b[j];
      b[i]=(b[i]-sum)/a[i][i];
    }
  *iflag=1;
  return;

}

void newtonSolve(double f[7][3],double *u1,double *v1,double *w1)
{
  int i,j,k;
  int iter,itmax,isolflag;
  double u,v,w;
  double uv,wu,vw,uvw,norm,convergenceLimit;
  double *rhs;
  double **lhs;
  double alph;
  //
  lhs=(double **)malloc(sizeof(double)*3);
  for(i=0;i<3;i++)
    lhs[i]=(double *)malloc(sizeof(double)*3);
  rhs=(double *)malloc(sizeof(double)*3);
  //
  itmax=500;
  convergenceLimit=1e-14;
  alph=1.0;
  isolflag=1.0;
  //
  u=v=w=0.5;
  //
  for(iter=0;iter<itmax;iter++)
    {
      uv=u*v;
      vw=v*w;
      wu=w*u;
      uvw=u*v*w;
      
      for(j=0;j<3;j++)
	rhs[j]=f[0][j]+f[1][j]*u+f[2][j]*v+f[3][j]*w+
	  f[4][j]*uv + f[5][j]*vw + f[6][j]*wu +
	  f[7][j]*uvw;
      
      norm=rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2];
      if (sqrt(norm) <= convergenceLimit) break;

      for(j=0;j<3;j++)
	{
	  lhs[j][0]=f[1][j]+f[4][j]*v+f[6][j]*w+f[7][j]*vw;
	  lhs[j][1]=f[2][j]+f[5][j]*w+f[4][j]*u+f[7][j]*wu;
	  lhs[j][2]=f[3][j]+f[6][j]*u+f[5][j]*v+f[7][j]*uv;
	}      
      
      solvec(lhs,rhs,&isolflag,3);
      if (isolflag==0) break;
      
      u-=(rhs[0]*alph);
      v-=(rhs[1]*alph);
      w-=(rhs[2]*alph);
    }
  if (iter==itmax) {u=2.0;v=w=0.;}
  if (isolflag==0) {
    u=2.0;
    v=w=0.;
  }
  *u1=u;
  *v1=v;
  *w1=w;
  for(i=0;i<3;i++) TIOGA_FREE(lhs[i]);
  TIOGA_FREE(lhs);
  TIOGA_FREE(rhs);
  return;
}


void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert)
{
  int i,j,k,isolflag;
  double **lhs;
  double *rhs;
  double f[8][3];
  double u,v,w;
  double oneminusU,oneminusV,oneminusW,oneminusUV;
 
  switch(nvert)
    {
    case 4:
      //
      // tetrahedron
      //
      lhs=(double **)malloc(sizeof(double)*3);
      for(i=0;i<3;i++)
	lhs[i]=(double *)malloc(sizeof(double)*3);
      rhs=(double *)malloc(sizeof(double)*3);       
      for(k=0;k<3;k++)
	{
	  for(j=0;j<3;j++)
	    lhs[j][k]=xv[k][j]-xv[3][j];
	  rhs[k]=xp[k]-xv[3][k];
	}
      //
      // invert the 3x3 matrix
      //
      solvec(lhs,rhs,&isolflag,3);
      //
      // check if the solution is not degenerate
      //
      if (isolflag) 
	{
	  for(k=0;k<3;k++) frac[k]=rhs[k];
	  frac[3]=1.-frac[0]-frac[1]-frac[2];
	}
      else
	{
	  frac[0]=1.0;
	  frac[1]=frac[2]=frac[3]=0;
	}
      for(i=0;i<3;i++) TIOGA_FREE(lhs[i]);
      TIOGA_FREE(lhs);
      TIOGA_FREE(rhs);
      break;
    case 5:
      //
      // pyramid
      //
      for(j=0;j<3;j++)
	{
	  f[0][j]=xv[0][j]-xp[j];
	  f[1][j]=xv[1][j]-xv[0][j];
	  f[2][j]=xv[3][j]-xv[0][j];
	  f[3][j]=xv[4][j]-xv[0][j];
	  //
	  f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
	  f[5][j]=xv[0][j]-xv[3][j];
	  f[6][j]=xv[0][j]-xv[1][j];
	  f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j];
	}
      //
      newtonSolve(f,&u,&v,&w);
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=w;
      //
      break;
    case 6:
      //
      // prizm
      //
      for(j=0;j<3;j++)
	{
	  f[0][j]=xv[0][j]-xp[j];
	  f[1][j]=xv[1][j]-xv[0][j];
	  f[2][j]=xv[2][j]-xv[0][j];
	  f[3][j]=xv[3][j]-xv[0][j];
	  //
	  f[4][j]=0;
	  f[5][j]=xv[0][j]-xv[2][j]-xv[3][j]+xv[5][j];
	  f[6][j]=xv[0][j]-xv[1][j]-xv[3][j]+xv[4][j];
	  f[7][j]=0.;
	}
      //
      newtonSolve(f,&u,&v,&w);
      //
      oneminusUV=1.0-u-v;
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusUV*oneminusW;
      frac[1]=u*oneminusW;
      frac[2]=v*oneminusW;
      frac[3]=oneminusUV*w;
      frac[4]=u*w;
      frac[5]=v*w;
      //
      break;
    case 8:
      //
      // hexahedra
      //
      for(j=0;j<3;j++)
	{
	  f[0][j]=xv[0][j]-xp[j];
	  f[1][j]=xv[1][j]-xv[0][j];
	  f[2][j]=xv[3][j]-xv[0][j];
	  f[3][j]=xv[4][j]-xv[0][j];
	  //
	  f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
	  f[5][j]=xv[0][j]-xv[3][j]+xv[7][j]-xv[4][j];
	  f[6][j]=xv[0][j]-xv[1][j]+xv[5][j]-xv[4][j];
	  f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j]+
	    xv[4][j]-xv[5][j]+xv[6][j]-xv[7][j];
	}
      //
      newtonSolve(f,&u,&v,&w);
      //
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=oneminusU*oneminusV*w;
      frac[5]=u*oneminusV*w;
      frac[6]=u*v*w;
      frac[7]=oneminusU*v*w;     
      //
      break;
    default:
      printf("Interpolation not implemented for polyhedra with %d vertices\n",nvert);
      break;
    }
}

void cellvolume_(double*, double[][3], int[][6], int[][24], int*, int*);

double computeCellVolume(double xv[8][3],int nvert)
{
 double vol;
 int itype;
 int nfaces;
 int numverts[4][6]={3,3,3,3,0,0,4,3,3,3,3,0,3,4,4,4,3,0,4,4,4,4,4,4};
 int faceInfo[4][24]={1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0,
                       1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0,
                       1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0,
                       1,2,3,4,1,5,6,2,2,6,7,3,3,7,8,4,1,4,8,5,5,8,7,6};
 switch(nvert)
   {
   case 4:
     itype=0;
     nfaces=4;
     break;
   case 5:
     itype=1;
     nfaces=5;
     break;
   case 6:
     itype=2;
     nfaces=5;
     break;
   case 8:
     itype=3;
     nfaces=6;
     break;
   }
     
 cellvolume_(&vol,xv,&numverts[itype],&faceInfo[itype],&nfaces,&nvert);
 return vol;
}
double tdot_product(double a[3],double b[3],double c[3])
{
  int k;
  double dp=0.0;
  for(k=0;k<3;k++) dp+=((a[k]-c[k])*(b[k]-c[k]));
  return dp;
}
