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
#include "MeshBlock.h"
extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  void writebbox(OBB *obb,int bid);
  void writePoints(double *x,int nsearch,int bid);
  void uniquenodes(double *x,double *rtag,int *itag,int *nn);
}


void MeshBlock::search(void)
{
  int i,j,k,l,m,n,p,i3;
  int ndim;
  int iptr,isum,nvert;
  OBB *obq;
  int *icell;
  int *itag;
  int cell_count; 
  int cellindex;
  double xd[3];
  double dxc[3];
  double xmin[3];
  double xmax[3];
  //
  // form the bounding box of the 
  // query points
  //

  if (nsearch == 0) {
    donorCount=0;
    return;
  }

  obq=(OBB *) malloc(sizeof(OBB));
  
findOBB(xsearch,obq->xc,obq->dxc,obq->vec,nsearch);


  //writebbox(obq,4);
  //writePoints(xsearch,nsearch,4);
  //
  // find all the cells that may have intersections with
  // the OBB
  //
  icell=(int *)malloc(sizeof(int)*ncells);
  for(i=0;i<ncells;i++) icell[i]=-1;
  iptr=-1;
  cell_count=0;
  p=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  //
	  // find each cell that has
	  // overlap with the bounding box
	  //
	  xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
	  xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
	  for(m=0;m<nvert;m++)
	    {
	      i3=3*(vconn[n][nvert*i+m]-BASE);	      
	      for(j=0;j<3;j++)
		{
		  xd[j]=0;
		  for(k=0;k<3;k++)
		    xd[j]+=(x[i3+k]-obq->xc[k])*obq->vec[j][k];
		  xmin[j]=MIN(xmin[j],xd[j]);
		  xmax[j]=MAX(xmax[j],xd[j]);
		}
	      for(j=0;j<3;j++)
		{
		  xd[j]=(xmax[j]+xmin[j])*0.5;
		  dxc[j]=(xmax[j]-xmin[j])*0.5;
		}
	    }
	  if (fabs(xd[0]) <= (dxc[0]+obq->dxc[0]) &&
	      fabs(xd[1]) <= (dxc[1]+obq->dxc[1]) &&
	      fabs(xd[2]) <= (dxc[2]+obq->dxc[2])) 
	    {
	      //
	      // create a LIFO stack
	      // with all the cells that 
	      // have bounding box intersection with
	      // the QP bounding box
	      //
	      icell[p]=iptr;
	      iptr=p;
	      cell_count++;
	    }
	  p++;
	}
    }
  //
  // now find the axis aligned bounding box
  // of each cell in the LIFO stack to build the
  // ADT
  //

  if (elementBbox) free(elementBbox);
  if (elementList) free(elementList);
  elementBbox=(double *)malloc(sizeof(double)*cell_count*6);
  elementList=(int *)malloc(sizeof(int)*cell_count);
  //
  k=iptr;
  l=0;
  p=0;
  //for(k=0;k<ncells;k++)
  while(k!=-1)
    {
      cellindex=k;
      isum=0;
      for(n=0;n<ntypes;n++) 
	{
	  isum+=nc[n];
	  if (cellindex < isum)
	    {
	      i=cellindex-(isum-nc[n]);
	      break;
	    }
	}
      nvert=nv[n];
      xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
      xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
      for(m=0;m<nvert;m++)
	{
	  i3=3*(vconn[n][nvert*i+m]-BASE);
	  for(j=0;j<3;j++)
	    {
	      xmin[j]=MIN(xmin[j],x[i3+j]);
	      xmax[j]=MAX(xmax[j],x[i3+j]);
	    }
	}
      //
      elementBbox[l++]=xmin[0];
      elementBbox[l++]=xmin[1];
      elementBbox[l++]=xmin[2];
      elementBbox[l++]=xmax[0];
      elementBbox[l++]=xmax[1];
      elementBbox[l++]=xmax[2];
      //
      elementList[p++]=k;
      //
      k=icell[k];
    }
  //
  // build the ADT now
  //
  if (adt) 
   {
    adt->clearData();
   }
  else
   {
    adt=new ADT[1];
   }
  ndim=6;
  //
  adt->buildADT(ndim,cell_count,elementBbox);
  //
  if (donorId) free(donorId);
  donorId=(int*)malloc(sizeof(int)*nsearch);
  if (xtag) free(xtag);
  xtag=(int *)malloc(sizeof(int)*nsearch);
  if (res_search0) free(res_search0);
  res_search0=(double *)malloc(sizeof(double)*nsearch);
  //
  // create a unique hash
  //
  for(i=0;i<nsearch;i++) res_search0[i]=res_search[i];
  uniquenodes(xsearch,res_search,xtag,&nsearch);
  //
  donorCount=0;
  ipoint=0; 
  for(i=0;i<nsearch;i++)
    {
      if (xtag[i]==i) {
	adt->searchADT(this,&(donorId[i]),&(xsearch[3*i]));
      }
      else {
	donorId[i]=donorId[xtag[i]];
      }
      if (donorId[i] > -1) {
	  donorCount++;
	}
      ipoint+=3;
    }
  //
  free(icell);
  free(obq);
}
