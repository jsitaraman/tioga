
/* This file is part of the Tioga software library */

/* Tioga  is a tool for overset grid assembly on parallel distributed systems */
/* Copyright (C) 2015 Jay Sitaraman */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA */
#include "codetypes.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void kaiser_wrap_(double *,int *,int *,double *,double *,double *,int *);

/***
 ** find oriented bounding box for a given set of points
*/
void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes)
{
  int i,j,k,m,i3;
  double *aa;
  double *eigenv;
  double trace,sume;
  int ier;
  double xd[3];
  double xmin[3],xmax[3];
  int nrows,ncols;
  printf("nnodes = %d, xend: %f\n",nnodes,x[nnodes*3-1]);
  //
  xc[0]=xc[1]=xc[2]=0;
  //
  // find centroid coordinates (not true centroid)
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      xc[0]+=x[i3];
      xc[1]+=x[i3+1];
      xc[2]+=x[i3+2];
    }
  //
  xc[0]/=nnodes;
  xc[1]/=nnodes;
  xc[2]/=nnodes;
  if (nnodes <4) 
    {
      vec[0][0]=vec[1][1]=vec[2][2]=1;
      vec[0][1]=vec[1][0]=vec[1][2]=vec[2][1]=0;
      vec[0][2]=vec[2][0]=0;
      if (nnodes==1) 
	{
	  dxc[0]=1e-3;	
	  dxc[1]=1e-3;
	  dxc[2]=1e-3;
	  return;
	}
      else if (nnodes==2)
	{
	  dxc[0]=max(1e-3,fabs(x[3]-x[0]))*0.5;
	  dxc[1]=max(1e-3,fabs(x[4]-x[1]))*0.5;
	  dxc[2]=max(1e-3,fabs(x[5]-x[2]))*0.5;
          return;
	}
      else
        {
         for(i=0;i<nnodes;i++)
          {
           i3=3*i;
           for(j=0;j<3;j++)
            dxc[j]=max(1e-3,fabs(x[i3+j]-x[0]));
          }
	 return;
        }
    }     
  //
  // find co-variance matrix
  // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
  //
  aa=(double *) malloc(sizeof(double)*9);
  eigenv=(double *)malloc(sizeof(double)*3);
  //
  for(i=0;i<9;i++) aa[i]=0;
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      aa[0]+=((x[i3]-xc[0])*(x[i3]-xc[0]));
      aa[4]+=((x[i3+1]-xc[1])*(x[i3+1]-xc[1]));
      aa[8]+=((x[i3+2]-xc[2])*(x[i3+2]-xc[2]));
      aa[3]+=((x[i3]-xc[0])*(x[i3+1]-xc[1]));
      aa[6]+=((x[i3]-xc[0])*(x[i3+2]-xc[2]));
      aa[7]+=((x[i3+1]-xc[1])*(x[i3+2]-xc[2]));
    }
  aa[1]=aa[3];
  aa[2]=aa[6];
  aa[5]=aa[7];
  //
  // use kaisers method to estimate
  // eigen values and vectors of the covariance matrix
  //
  nrows=3;
  ncols=3;  
  kaiser_wrap_(aa,&nrows,&ncols,eigenv,&trace,&sume,&ier);
  //
  // copy the eigen vector basis on to vec
  //
  m=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
       vec[i][j]=aa[m++];
      }
  //
  // find min and max bounds in the bounding box
  // vector basis
  //
  for(j=0;j<3;j++)
    {
      xmax[j]=-BIGVALUE;
      xmin[j]=BIGVALUE;
    }
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      //
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(x[i3+k]-xc[k])*vec[j][k];
      //
      for(j=0;j<3;j++)
	{
	  xmax[j]=max(xmax[j],xd[j]);
	  xmin[j]=min(xmin[j],xd[j]);
	}
    }
  //
  // find the extents of the box
  // and coordinates of the center w.r.t. xc
  // increase extents by 1% for tolerance
  //
  for(j=0;j<3;j++)
    {
      dxc[j]=(xmax[j]-xmin[j])*0.5*1.01;
      xd[j]=(xmax[j]+xmin[j])*0.5;
    }
  //
  // find the center of the box in 
  // actual cartesian coordinates
  //
  for(j=0;j<3;j++)
    {
    for(k=0;k<3;k++)
      xc[j]+=(xd[k]*vec[k][j]);
    }
  //
  free(aa);
  free(eigenv);
}

/** Check if a point is inside the provided hole map (return hole status) */
int checkHoleMap(double *x,int *nx,int *sam,double *extents)
{
  double dx[3];
  int ix[3];

  for (int i = 0; i < 3; i++)
    dx[i] = (extents[i+3]-extents[i]) / nx[i];

  for (int i = 0; i < 3; i++)
  {
    ix[i] = (x[i]-extents[i])/dx[i];
    if (ix[i] < 0 || ix[i] > nx[i]-1) return 0;
  }

  int mm = (ix[2]*nx[1] + ix[1])*nx[0] + ix[0];
  return sam[mm];
}

/**
 fill a given hole map using iterative
 flood fill from outside the marked boundary.
 boundary is marked by "2"
*/      
void fillHoleMap(int *holeMap, int ix[3],int isym)
{
  int m;
  int ii,jj,kk,mm;
  int ns2;
  int i,j,k;
  int ipaint,npaint,nneig;
  int mk[6];
  //
  // now start from outside and paint the
  // the exterior
  //
  ns2=ix[0]*ix[1];
  //
  for(kk=0;kk<ix[2];kk+=(ix[2]-1))
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii++)
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj+=(ix[1]-1))
      for(ii=0;ii<ix[0];ii++)
        {
  	  mm=kk*ns2+jj*ix[0]+ii;
         holeMap[mm]=1;
  	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii+=(ix[0]-1))
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}  
  npaint=ns2*ix[2];
  while(npaint > 0) 
    {
      npaint=0;
      for(k=1;k<ix[2]-1;k++)
        for(j=1;j<ix[1]-1;j++)
          for(i=1;i<ix[0]-1;i++)
            {
              m=k*ns2+j*ix[0]+i;
              if (holeMap[m]==0)
                {
                  ipaint=0;
                  if (isym==1) 
                   {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[2]=m-ix[0];
		     mk[3]=m+ix[0];
		     mk[4]=m-1;
		     mk[5]=m+1;
		     nneig=4;
		   }
		  else if (isym==2)
		    {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[4]=m-ix[0];
		     mk[5]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else if (isym==3)
		    {
		     mk[4]=m-ns2;
		     mk[5]=m+ns2;
		     mk[0]=m-ix[0];
		     mk[1]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else
		    {
		      mk[0]=m-ns2;
		      mk[1]=m+ns2;
		      mk[2]=m-ix[0];
		      mk[3]=m+ix[0];
		      mk[4]=m-1;
		      mk[5]=m+1;
		      nneig=6;
		    }		 
                  for (kk=0;kk<nneig && ipaint==0;kk++)
		    {
		      ipaint=(ipaint || holeMap[mk[kk]]==1);
		    }
                  if (ipaint > 0)
                    {
                      holeMap[m]=1;
                      npaint++;
                    }
                }
            }
    }
  for(i=0;i<ix[2]*ix[1]*ix[0];i++) 
   { 
    holeMap[i]=(holeMap[i] ==0 || holeMap[i]==2);
   }
}

int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
		       double vB[3][3],double xB[3],double dxB[3])
{
  int iflag;
  int i,j,k;
  int i1,i2,j1,j2;
  double r,r0,r1;
  double d1,d2;
  double eps=1e-12;
  double D[3];
  double c[3][3];
  //
  // D=distance between centers
  // C=scalar product of axes
  //
  for(i=0;i<3;i++) D[i]=xB[i]-xA[i];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	c[i][j]=0;
	for(k=0;k<3;k++)
	  c[i][j]=c[i][j]+vA[i][k]*vB[j][k];
      }
  //
  // separating axes based on the faces of box A
  //
  for(i=0;i<3;i++)
    {
      r0=dxA[i];
      r1=0;
      r=0;
      for(j=0;j<3;j++) 
	{
	  r1+=dxB[j]*fabs(c[i][j]);
	  r+=fabs(vA[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // separating axes based on the faces of box B
  //
  for(i=0;i<3;i++)
    {
      r1=dxB[i];
      r0=0;
      r=0;
      for(j=0;j<3;j++) 
	{
	  r0+=dxA[j]*fabs(c[j][i]);
	  r+=fabs(vB[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // cross products
  //
  for(i=0;i<3;i++)
    {
      i1=(i+1)%3;
      i2=(i+2)%3;
      for(j=0;j<3;j++)
	{
	  j1=(j+1)%3;
	  j2=(j+2)%3;
	  
	  r0=dxA[i1]*fabs(c[i2][j])+dxA[i2]*fabs(c[i1][j]);
	  r1=dxB[j1]*fabs(c[i][j2])+dxB[j2]*fabs(c[i][j1]);
	  
	  d2=0;
	  d1=0;
	  for(k=0;k<3;k++)
	    {
	      d2+=vA[i2][k]*D[k];
	      d1+=vA[i1][k]*D[k];
	    }
	  
	  r=fabs(c[i1][j]*d2-c[i2][j]*d1);
	  
	  if (r > (r0+r1+eps)) {
	    return 0;
	  }
	}
    }
  //
  // return zero if no separation can be found
  //
  return 1;
}
    

	  
void writebbox(OBB *obb,int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"qbox%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,
	  1);

  for(l=0;l<2;l++)
    {
      il=2*(l%2)-1;
      for(k=0;k<2;k++)
	{
	  ik=2*(k%2)-1;
	  for(j=0;j<2;j++)
	    {
	      ij=2*(j%2)-1;
	      xx[0]=xx[1]=xx[2]=0;
	      for(m=0;m<3;m++)
		xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
		  +ik*obb->vec[1][m]*obb->dxc[1]
		  +il*obb->vec[2][m]*obb->dxc[2];	      
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	    }
	}
    }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fclose(fp);
}
    
      
    
void writePoints(double *x,int nsearch,int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int i;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"points%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  for(i=0;i<nsearch;i++)
    fprintf(fp,"%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);
  fclose(fp);
}

#ifdef __cplusplus
}
#endif
