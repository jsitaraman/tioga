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
#include "linklist.h"
#include "utils.h"

extern "C"{
  void get_amr_index_xyz(int nq,int i,int j,int k,
			 int pBasis,
			 int nX,int nY,int nZ,
			 int nf,
			 double *xlo,double *dx,
			 double *qnodes,
			 int* index, double* xyz);
}

void MeshBlock::getCartReceptors(CartGrid *cg,parallelComm *pc,int itype)
{
  int *sndMap,*rcvMap;
  double xd[3];
  int intersectCount=0;

  // limit case we communicate to everybody
  std::vector<int> pmap(pc->numprocs);
  
  OBB* obcart = (OBB *)malloc(sizeof(OBB));
  for (int j = 0;j<3;j++)
    for (int k = 0;k<3;k++)
      obcart->vec[j][k]=0;
  obcart->vec[0][0]=obcart->vec[1][1]=obcart->vec[2][2]=1.0;
  
  INTEGERLIST2* head = (INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
  head->intData = NULL;
  head->realData = NULL;

  INTEGERLIST2* dataPtr = head;
  dataPtr->next = NULL;
  
  nsearch=0;
  
  //writeOBB(myid);
  
  for (int c = 0; c < cg->ngrids; c++)
  {
    for (int n = 0; n < 3; n++)
    {
      obcart->dxc[n]=cg->dx[3*c+n]*(cg->dims[3*c+n])*0.5;
      obcart->xc[n]=cg->xlo[3*c+n]+obcart->dxc[n];
    }

    if (obbIntersectCheck(obb->vec,obb->xc,obb->dxc,
          obcart->vec,obcart->xc,obcart->dxc) ||
        obbIntersectCheck(obcart->vec,obcart->xc,obcart->dxc,
          obb->vec,obb->xc,obb->dxc))
    {
      intersectCount++;
      //if (myid==0 && intersectCount==0) writeOBB2(obcart,c);
      int ntm = 1;
      int iadd = 0;
      if (itype == 0) 
      {         
        ntm = (cg->porder[c]+1)*(cg->porder[c]+1)*(cg->porder[c]+1);      
        iadd = 0;
      }

      std::vector<double> xtm(3*ntm);
      std::vector<int> itm(ntm);

      int ploc = (cg->porder[c])*(cg->porder[c]+1)/2;
      for (int j = 0; j < cg->dims[3*c]+iadd; j++)
        for (int k = 0; k < cg->dims[3*c+1]+iadd; k++)
          for (int l = 0; l < cg->dims[3*c+2]+iadd; l++)
          {
            if (itype==0) 
            {
              get_amr_index_xyz(cg->qstride,j,k,l,
                  cg->porder[c],cg->dims[3*c],cg->dims[3*c+1],cg->dims[3*c+2],
                  cg->nf,
                  &cg->xlo[3*c],
                  &cg->dx[3*c],
                  &cg->qnode[ploc],
                  itm.data(),
                  xtm.data());
            }
            else
            {
              xtm[0]=cg->xlo[3*c  ]+cg->dx[3*c  ]*j;
              xtm[1]=cg->xlo[3*c+1]+cg->dx[3*c+1]*k;
              xtm[2]=cg->xlo[3*c+2]+cg->dx[3*c+2]*l;
              int nf=cg->nf;
              itm[0]=l*(cg->dims[3*c])*(cg->dims[3*c+1])+k*(cg->dims[3*c])+j;
              //itm[0]=(l+nf)*(cg->dims[3*c]+2*nf)*(cg->dims[3*c+1]+2*nf)+(k+nf)*(cg->dims[3*c]+2*nf)+j+nf;
            }   

            int iflag = 0;

            for (int n = 0; n < ntm; n++)
            {
              int i3 = 3*n;
              //if (intersectCount==1 && myid==0) fprintf(fp,"%lf %lf %lf\n",xtm[i3],xtm[i3+1],xtm[i3+2]);
              for (int jj = 0; jj < 3; jj++) xd[jj]=0;
              for (int jj = 0; jj < 3; jj++)
                for (int kk = 0; kk < 3; kk++)
                  xd[jj]+=(xtm[i3+kk]-obb->xc[kk])*obb->vec[jj][kk];

              if (fabs(xd[0]) <= obb->dxc[0] &&
                  fabs(xd[1]) <= obb->dxc[1] &&
                  fabs(xd[2]) <= obb->dxc[2])
              {
                iflag++;
              }
            }

            if (iflag > 0) 
            {
              pmap[cg->proc_id[c]] = 1;
              dataPtr->next = (INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
              dataPtr = dataPtr->next;
              dataPtr->realData = (double *) malloc(sizeof(double)*3*ntm); /// 3->4 : adding nodeRes
              dataPtr->intData  = (int *) malloc(sizeof(int)*(ntm+2));
              dataPtr->intData[0] = cg->proc_id[c];
              dataPtr->intData[1] = cg->local_id[c];
              dataPtr->intDataSize  = ntm+2;
              dataPtr->realDataSize = ntm;
              nsearch += ntm;
              for (int n = 0; n < ntm; n++)
              {
                for (int kk = 0; kk < 3; kk++)
                  dataPtr->realData[3*n+kk] = xtm[3*n+kk];
                //dataPtr->realData[4*n+3] = cg->dx[3*c]*cg->dx[3*c+1]*cg->dx[3*c+2]; // nodeRes / cell volume (dx^3)
                dataPtr->intData[n+2] = itm[n];
              }
              dataPtr->next=NULL;
            }
          }
    }
  }
  
  // create the communication map
  int nsend=0;
  for (int i = 0; i < pc->numprocs; i++) 
    if (pmap[i]==1) nsend++;
  int nrecv = nsend;

  sndMap = (int *)malloc(sizeof(int)*nsend);
  rcvMap = (int *)malloc(sizeof(int)*nrecv);

  int m = 0;
  for (int i = 0; i < pc->numprocs; i++)
  {
    if (pmap[i]==1) 
    {
      sndMap[m]=rcvMap[m]=i;
      m++;
    }
  }
  pc->setMap(nsend,nrecv,sndMap,rcvMap);

  xsearch.resize(3*nsearch);
  isearch.resize(3*nsearch);
  donorId.resize(nsearch);
  rst.resize(3*nsearch);

  dataPtr=head->next;
  int n = 0, p = 0;
  m = 0;
  while(dataPtr!=NULL)
  {
    for (int j = 2; j < dataPtr->intDataSize; j++)
    {
      isearch[m++] = dataPtr->intData[0];
      isearch[m++] = dataPtr->intData[1];	  
      isearch[m++] = dataPtr->intData[j];
    }

    for (int j = 0; j < dataPtr->realDataSize; j++)
    {
      // Point locations
      for (int k = 0; k < 3; k++)
        xsearch[n++] = dataPtr->realData[3*j+k]; //3->4
    }

    dataPtr = dataPtr->next;
  }

  deallocateLinkList3(head);
  free(obcart);
  free(sndMap);
  free(rcvMap);
}
