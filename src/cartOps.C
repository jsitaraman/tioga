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

void MeshBlock::getUnresolvedMandatoryReceptors(void)
{
  int i,j,k,m,n,nvert,i3,fcount;
  int inode[8];
  int *iflag;

  iflag=(int *) malloc(sizeof(int) *ncells);
  if (pickedCart !=NULL) free(pickedCart);
  pickedCart=(int *) malloc(sizeof(int)*nnodes);

  for(i=0;i<ncells;i++) iflag[i]=0;
  for(i=0;i<nnodes;i++) pickedCart[i]=0;

  k=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  fcount=0;
	  for(m=0;m<nvert;m++)
	    {
	      inode[m]=vconn[n][nvert*i+m]-BASE;
	      if (nodeRes[inode[m]]==BIGVALUE) fcount++;
	    }
	  if (fcount==nvert && iblank_cell[k]==1) 
	    {
	      iflag[k]=1;
	      for(m=0;m<nvert;m++)
		pickedCart[inode[m]]=1;
	    }
	  k++;
	}
    }
  //
  if (ctag_cart!=NULL) free(ctag_cart);
  ctag_cart=(int *)malloc(sizeof(int)*ncells);
  nreceptorCellsCart=0;
  for(i=0;i<ncells;i++)
    if (iflag[i]==-1) ctag_cart[nreceptorCellsCart++]=i+1;
  //
  if (ihigh) 
    {
      if (pointsPerCell!=NULL) free(pointsPerCell);
      pointsPerCell=(int *)malloc(sizeof(int)*nreceptorCellsCart);
      //
      maxPointsPerCell=0;
      ntotalPointsCart=0;
      //
      for(i=0;i<nreceptorCellsCart;i++)
	{
	  get_nodes_per_cell(&(ctag_cart[i]),&(pointsPerCell[i]));
	  ntotalPointsCart+=pointsPerCell[i];
	  maxPointsPerCell=max(maxPointsPerCell,pointsPerCell[i]);
      }
      //
      if (rxyzCart !=NULL) free(rxyzCart);
      //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
      rxyzCart=(double *)malloc(sizeof(double)*ntotalPointsCart*3);
      //
      m=0;
      for(i=0;i<nreceptorCellsCart;i++)
	{
	  get_receptor_nodes(&(ctag_cart[i]),&(pointsPerCell[i]),&(rxyzCart[m]));
	  m+=(3*pointsPerCell[i]);
	}
    }
  else
    {
      ntotalPointsCart=0;      
      for(i=0;i<nnodes;i++) 
	{
	  if (pickedCart[i]==-1) 
	    {
	      pickedCart[i]=1;
	      ntotalPointsCart++;
	    }
	}
      if (rxyzCart !=NULL) free(rxyzCart);
      rxyzCart=(double *)malloc(sizeof(double)*ntotalPointsCart*3);
      m=0;
      for(i=0;i<nnodes;i++)
	if (pickedCart[i]) 
	  {
	    i3=3*i;
	    for(j=0;j<3;j++)
	      rxyzCart[m++]=x[i3+j];
	  }
    }
 free(iflag);
}
