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
//
// routine for extra query points
// have to unify both routines here 
// FIX later ...
//
#include "MeshBlock.h"
#include "parallelComm.h"
#include "CartGrid.h"
extern "C"{
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
                        double vB[3][3],double xB[3],double dxB[3]);

  void get_amr_index_xyz(int nq,int i,int j,int k,
			 int pBasis,
			 int nX,int nY,int nZ,
			 int nf,
			 double *xlo,double *dx,
			 double *qnodes,
			 int* index, double* xyz);
  void deallocateLinkList3(INTEGERLIST2 *);
}

void MeshBlock::getCartReceptors(CartGrid *cg,parallelComm *pc)
{
  int i,j,k,l,m,ploc,c,n,ntm;
  int i3;
  int iflag;
  int icount,dcount;
  int nsend,nrecv;
  int *pmap;
  int *sndMap,*rcvMap;
  OBB *obcart;
  INTEGERLIST2 *head;
  INTEGERLIST2 *dataPtr;
  int *itm;
  double *xtm;
  double xd[3];
  //
  // limit case we communicate to everybody
  //
  pmap=(int *)malloc(sizeof(int)*pc->numprocs);
  for(i=0;i<pc->numprocs;i++) pmap[i]=0;
  //
  obcart=(OBB *)malloc(sizeof(OBB));
  for(j=0;j<3;j++)
    for(k=0;k<3;k++)
      obcart->vec[j][k]=0;
  obcart->vec[0][0]=obcart->vec[1][1]=obcart->vec[2][2]=1.0;
  //
  head=(INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
  dataPtr=head;
  //
  nsearch=0;
  //
  for(c=0;c<cg->ngrids;c++)
    {
      for(n=0;n<3;n++)
	{
	  obcart->dxc[n]=cg->dx[3*c+n]*cg->dims[3*c+n]*0.5;
	  obcart->xc[n]=cg->xlo[3*c+n]+obcart->dxc[n];
	}
      if (obbIntersectCheck(obb->vec,obb->xc,obb->dxc,
			    obcart->vec,obcart->xc,obcart->dxc) ||
	  obbIntersectCheck(obcart->vec,obcart->xc,obcart->dxc,
			    obb->vec,obb->xc,obb->dxc))
	{
	  ntm=(cg->porder[c]+1)*(cg->porder[c]+1)*(cg->porder[c]+1);      
	  xtm=(double *)malloc(sizeof(double)*3*ntm);
	  itm=(int *) malloc(sizeof(int)*ntm);
	  ploc=(cg->porder[c])*(cg->porder[c]+1)/2;
	  for(j=0;j<cg->dims[3*c];j++)
	    for(k=0;k<cg->dims[3*c+1];k++)
	      for(l=0;l<cg->dims[3*c+2];l++)
		{
		  get_amr_index_xyz(cg->qstride,j,k,l,
				    cg->porder[c],cg->dims[3*c],cg->dims[3*c+1],cg->dims[3*c+2],
				    cg->nf,
				    &cg->xlo[3*c+n],
				    &cg->dx[3*c+n],
				    &cg->qnode[ploc],
				    itm,
				    xtm);
		  iflag=0;
		  for(n=0;n<ntm;n++)
		    {
		      i3=3*n;
		      for(j=0;j<3;j++) xd[j]=0;
		      for(j=0;j<3;j++)
			for(k=0;k<3;k++)
			  xd[j]+=(xtm[i3+k]-obb->xc[k])*obb->vec[j][k];
		      
		      if (fabs(xd[0]) <= obb->dxc[0] &&
			  fabs(xd[1]) <= obb->dxc[1] &&
			  fabs(xd[2]) <= obb->dxc[2])
			{
			  iflag++;
			}
		    }
		  
		  if (iflag==ntm) 
		    {
		      pmap[cg->proc_id[c]]=1;
		      dataPtr->next=(INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
		      dataPtr=dataPtr->next;
		      dataPtr->realData=(double *) malloc(sizeof(double)*3*ntm);
		      dataPtr->intData=(int *) malloc(sizeof(int)*ntm+2);
		      dataPtr->intData[0]=cg->proc_id[c];
		      dataPtr->intData[1]=cg->local_id[c];
		      dataPtr->intDataSize=ntm+2;
		      dataPtr->realDataSize=ntm;
		      nsearch+=ntm;
		      for(n=0;n<ntm;n++)
			{
			  for(k=0;k<3;k++)
			    dataPtr->realData[3*n+k]=xtm[3*n+k];
			  dataPtr->intData[n+2]=itm[n];
			}
		      dataPtr->next=NULL;
		    }
		}
	  free(xtm);
	  free(itm);
	}
    }
  //
  // create the communication map
  //
  nsend=0;
  for(i=0;i<pc->numprocs;i++) if (pmap[i]==1) nsend++;
  nrecv=nsend;
  sndMap=(int *)malloc(sizeof(int)*nsend);
  rcvMap=(int *)malloc(sizeof(int)*nrecv);
  m=0;
  for(i=0;i<pc->numprocs;i++)
    {
      if (pmap[i]==1) 
	{
	  sndMap[m]=rcvMap[m]=i;
	  m++;
	}
    }
  pc->setMap(nsend,nrecv,sndMap,rcvMap);
  //
  // if these were already allocated
  // get rid of them
  //
  if (xsearch) free(xsearch);
  if (isearch) free(isearch);
  if (donorId) free(donorId);
  if (rst) free(rst);
  //
  xsearch=(double *)malloc(sizeof(double)*3*nsearch);
  isearch=(int *)malloc(3*sizeof(int)*nsearch);
  donorId=(int *)malloc(sizeof(int)*nsearch);
  rst=(double *) malloc(sizeof(double)*3*nsearch);
  //
  dataPtr=head->next;
  m=n=0;
  while(dataPtr!=NULL)
    {
      for(j=2;j<dataPtr->intDataSize;j++)
	{
	  isearch[m++]=dataPtr->intData[0];
	  isearch[m++]=dataPtr->intData[1];	  
	  isearch[m++]=dataPtr->intData[j];
	}
      for(j=0;j<3*dataPtr->realDataSize;j++)
	xsearch[n++]=dataPtr->realData[j];
    }
  deallocateLinkList3(head);
  free(obcart);
  free(pmap);
  free(sndMap);
  free(rcvMap);
}

  
