
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
#include <cmath>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

extern void kaiser_wrap_(double *,int *,int *,double *,double *,double *,int *);

/***
 ** find oriented bounding box for a given set of points
*/
void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes)
{
  double trace,sume;
  int ier;
  double xd[3];
  double xmin[3],xmax[3];
  
  xc[0]=xc[1]=xc[2]=0;
  
  // find centroid coordinates (not true centroid)
  for (int i = 0; i < nnodes; i++)
  {
    int i3 = 3*i;
    xc[0] += x[i3];
    xc[1] += x[i3+1];
    xc[2] += x[i3+2];
  }
  
  xc[0] /= nnodes;
  xc[1] /= nnodes;
  xc[2] /= nnodes;

  if (nnodes < 4) 
  {
    vec[0][0] = vec[1][1] = vec[2][2] = 1;
    vec[0][1] = vec[1][0] = vec[1][2] = vec[2][1] = 0;
    vec[0][2] = vec[2][0] = 0;

    if (nnodes == 1) 
    {
      dxc[0] = 1e-3;	
      dxc[1] = 1e-3;
      dxc[2] = 1e-3;
      return;
    }
    else if (nnodes == 2)
    {
      dxc[0] = max(1e-3, fabs(x[3]-x[0])) * 0.5;
      dxc[1] = max(1e-3, fabs(x[4]-x[1])) * 0.5;
      dxc[2] = max(1e-3, fabs(x[5]-x[2])) * 0.5;
      return;
    }
    else
    {
      for(int i = 0; i < nnodes; i++)
      {
        int i3 = 3*i;
        for (int j = 0; j < 3; j++)
          dxc[j] = max(1e-3,fabs(x[i3+j]-x[0]));
      }
      return;
    }
  }     
  
  // find co-variance matrix
  // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
  double aa[9], eigenv[3];
  
  for (int i = 0; i < 9; i++) aa[i] = 0;
  
  for (int i = 0; i < nnodes; i++)
  {
    int i3 = 3*i;
    aa[0] += ((x[i3]-xc[0]) * (x[i3]-xc[0]));
    aa[4] += ((x[i3+1]-xc[1]) * (x[i3+1]-xc[1]));
    aa[8] += ((x[i3+2]-xc[2]) * (x[i3+2]-xc[2]));
    aa[3] += ((x[i3]-xc[0]) * (x[i3+1]-xc[1]));
    aa[6] += ((x[i3]-xc[0]) * (x[i3+2]-xc[2]));
    aa[7] += ((x[i3+1]-xc[1]) * (x[i3+2]-xc[2]));
  }
  aa[1] = aa[3];
  aa[2] = aa[6];
  aa[5] = aa[7];

  // use kaisers method to estimate
  // eigen values and vectors of the covariance matrix
  int nrows = 3;
  int ncols = 3;  
  kaiser_wrap_(&aa[0],&nrows,&ncols,&eigenv[0],&trace,&sume,&ier);

  // copy the eigen vector basis on to vec
  int m = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      vec[i][j] = aa[m++];
    }
  
  // find min and max bounds in the bounding box
  // vector basis
  for (int j = 0; j < 3; j++)
  {
    xmax[j] = -BIGVALUE;
    xmin[j] =  BIGVALUE;
  }
  for (int i = 0; i < nnodes; i++)
  {
    int i3 = 3*i;
    for(int j = 0; j < 3; j++) xd[j] = 0;

    for (int j = 0; j < 3; j++)
    	for (int k = 0; k < 3; k++)
	      xd[j]+=(x[i3+k]-xc[k])*vec[j][k];
      
    for (int j = 0; j < 3; j++)
	  {
	    xmax[j] = max(xmax[j],xd[j]);
	    xmin[j] = min(xmin[j],xd[j]);
	  }
  }

  // find the extents of the box
  // and coordinates of the center w.r.t. xc
  // increase extents by 5% for tolerance - NECESSARY FOR MOVING GRIDS!!
  for (int j = 0; j < 3; j++)
  {
    dxc[j] = (xmax[j] - xmin[j]) * 0.5 * 1.05;
    xd[j] = (xmax[j] + xmin[j]) * 0.5;
  }
  
  // find the center of the box in 
  // actual cartesian coordinates
  for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
      xc[j] += (xd[k]*vec[k][j]);
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
  int mk[6];

  // now start from outside and paint the
  // the exterior
  int ns2 = ix[0]*ix[1];

  for (int k = 0; k < ix[2]; k += (ix[2]-1))
    for (int j = 0; j < ix[1]; j++)
      for (int i = 0; i < ix[0]; i++)
        holeMap[k*ns2+j*ix[0]+i] = 1;

  for(int k = 0; k < ix[2]; k++)
    for(int j = 0; j < ix[1]; j += (ix[1]-1))
      for(int i = 0; i < ix[0]; i++)
        holeMap[k*ns2+j*ix[0]+i] = 1;

  for(int k = 0; k < ix[2]; k++)
    for(int j = 0; j < ix[1]; j++)
      for(int i = 0; i < ix[0]; i += (ix[0]-1))
        holeMap[k*ns2+j*ix[0]+i] = 1;

  int npaint = ns2*ix[2];
  while (npaint > 0)
  {
    npaint = 0;
    for(int k = 1; k < ix[2]-1; k++)
      for(int j = 1; j < ix[1]-1; j++)
        for(int i = 1; i < ix[0]-1; i++)
        {
          int ipaint = 1;
          int nneig = 0;
          int m = (k*ix[1]+j)*ix[0]+i;
          if (holeMap[m] == 0)
          {
            ipaint = 0;
            if (isym == 1)
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
            for (int kk = 0; (kk < nneig && ipaint == 0); kk++)
            {
              ipaint = (ipaint || holeMap[mk[kk]] == 1);
            }
            if (ipaint > 0)
            {
              holeMap[m] = 1;
              npaint++;
            }
          }
        }
  }

  // For Direct Cut - want to keep info about boundary-containing cells
  for (int i = 0; i < ix[2]*ix[1]*ix[0]; i++)
  {
    if (holeMap[i] != 2)
      holeMap[i] = !holeMap[i];
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

/*
 * Create a unique hash for list of coordinates with duplicates in 
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes(double *x,double *rtag,int *itag,int nn)
{
  int NSUB = 101;
  int nnodes = nn;

  double xmax[3], xmin[3];
  for (int j = 0; j < 3; j++) xmax[j] = -1E15;
  for (int j = 0; j < 3; j++) xmin[j] =  1E15;

  for (int i = 0; i < nnodes; i++)
  {
    for (int j = 0; j < 3; j++) 
    {
      xmax[j] = std::max(xmax[j], x[3*i+j]);
      xmin[j] = std::min(xmin[j], x[3*i+j]);
    }
  }

  double ds = (xmax[0]-xmin[0]+xmax[1]-xmin[1]+xmax[2]-xmin[2])/3.0/NSUB;
  double dsi = 1.0/ds;

  for (int j = 0; j < 3; j++) xmax[j] += ds;
  for (int j = 0; j < 3; j++) xmin[j] -= ds;

  int jmax = std::min((int)round((xmax[0]-xmin[0])*dsi),NSUB);
  double dsx = (xmax[0]-xmin[0])/jmax;
  double dsxi = 1./dsx;    
  int kmax = std::min((int)round((xmax[1]-xmin[1])*dsi),NSUB);
  double dsy = (xmax[1]-xmin[1])/kmax;
  double dsyi = 1./dsy;
  int lmax = std::min((int)round((xmax[2]-xmin[2])*dsi),NSUB);
  double dsz = (xmax[2]-xmin[2])/lmax;
  double dszi = 1./dsz;
  int nsblks = jmax*kmax*lmax;
  int jkmax = jmax*kmax;

  std::vector<int> cft(nsblks+1);
  std::vector<int> numpts(nsblks);
  std::vector<int> ilist(nnodes);

  for (int i = 0; i < nsblks; i++) numpts[i]=0;
  for (int i = 0; i < nnodes; i++)
  {
    int i3 = 3*i;
    int jj = (int)((x[i3]   - xmin[0])*dsxi);
    int kk = (int)((x[i3+1] - xmin[1])*dsyi);
    int ll = (int)((x[i3+2] - xmin[2])*dszi);
    int indx = ll*jkmax + kk*jmax + jj;
    numpts[indx] = numpts[indx] + 1;
  }

  cft[0] = 0;
  for (int i = 0; i < nsblks; i++) 
    cft[i+1] = cft[i] + numpts[i];

  for (int i = 0; i < nnodes ;i++)
  {
    int i3 = 3*i;
    int jj = (int)((x[i3]-xmin[0])*dsxi);
    int kk = (int)((x[i3+1]-xmin[1])*dsyi);
    int ll = (int)((x[i3+2]-xmin[2])*dszi);
    int indx = ll*jkmax + kk*jmax + jj;
    ilist[cft[indx] + numpts[indx]-1] = i;
    numpts[indx]--;
    itag[i] = i;
  }

  for (int i = 0; i < nsblks; i++)
  {
    for (int j = cft[i]; j < cft[i+1]; j++)
    {
      int p1 = ilist[j];
      for (int k = j+1; k < cft[i+1]; k++)
      {
        int p2 = ilist[k];
        if ( fabs(x[3*p1  ]-x[3*p2  ]) 
           + fabs(x[3*p1+1]-x[3*p2+1]) 
           + fabs(x[3*p1+2]-x[3*p2+2]) < TOL )
        {
          if (p1 > p2) 
          {
            rtag[p2] = std::max(rtag[p1],rtag[p2]);
            itag[p1] = p2;
          }
          else 
          {
            rtag[p1] = std::max(rtag[p1],rtag[p2]);
            itag[p2] = p1;
          }
        }
      }
    }
  }
}
