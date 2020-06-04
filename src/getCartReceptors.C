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
#include "codetypes.h"
#include "MeshBlock.h"
#include "parallelComm.h"
#include "CartGrid.h"
extern "C"{
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
                        double vB[3][3],double xB[3],double dxB[3]);
  void deallocateLinkList3(INTEGERLIST2 *);
}

void MeshBlock::getCartReceptors(CartGrid *cg,parallelComm *pc)
{
  int i,j,k,l,m,c,n,itm,jj,kk;
  int i3;
  int iflag;
  int icount,dcount,cell_count;
  int nsend,nrecv;
  int *pmap;
  int *sndMap,*rcvMap;
  OBB *obcart;
  INTEGERLIST2 *head;
  INTEGERLIST2 *dataPtr;
  double *xtm;
  double xd[3],vol;
  //char qstr[2];
  //char fname[80];
  //char intstring[7];
  //FILE *fp;
  int intersectCount=0;

  //sprintf(intstring,"%d",100000+myid);
  //sprintf(fname,"zsearch_%s.dat",&(intstring[1]));
  //fp=fopen(fname,"w");
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
  head->intData=NULL;
  head->realData=NULL;
  dataPtr=head;
  dataPtr->next=NULL;
  //
  nsearch=0;
  //
  //writeOBB(myid);
  //
  for(c=0;c<cg->ngrids;c++)
  {
    cell_count = (cg->dims[3*c]+2*cg->nf)
        * (cg->dims[3*c+1]+2*cg->nf)
        * (cg->dims[3*c+2]+2*cg->nf);

    vol = cg->dx[3*c]*cg->dx[3*c+1]*cg->dx[3*c+2];

    for(n=0;n<3;n++)
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
           
	  xtm=(double *)malloc(sizeof(double)*3);

	  // loop over all cells
	  for(j=0;j<cg->dims[3*c];j++)
	    for(k=0;k<cg->dims[3*c+1];k++)
	      for(l=0;l<cg->dims[3*c+2];l++)
		{
	    //Q[nq,nZ+2*nf,nY+2*nf,nX+2*nf]--> C++ Cell storage
      itm = (cg->dims[3*c+1]+2*cg->nf)*(cg->dims[3*c]+2*cg->nf)*(l+cg->nf)
          + (cg->dims[3*c]+2*cg->nf)*(k+cg->nf) + (j+cg->nf);

      xtm[0] = cg->xlo[3*c]   + (j+0.5)*cg->dx[3*c];
      xtm[1] = cg->xlo[3*c+1] + (k+0.5)*cg->dx[3*c+1];
      xtm[2] = cg->xlo[3*c+2] + (l+0.5)*cg->dx[3*c+2];

		  iflag=0;

		    {
                      //if (intersectCount==1 && myid==0) fprintf(fp,"%lf %lf %lf\n",xtm[0],xtm[1],xtm[2]);
		      for(jj=0;jj<3;jj++) xd[jj]=0;
		      for(jj=0;jj<3;jj++)
			for(kk=0;kk<3;kk++)
			  xd[jj]+=(xtm[kk]-obb->xc[kk])*obb->vec[jj][kk];
		      
		      if (fabs(xd[0]) <= obb->dxc[0] &&
			  fabs(xd[1]) <= obb->dxc[1] &&
			  fabs(xd[2]) <= obb->dxc[2])
			{
			  iflag++;
			}
		    }
		  
		  if (iflag > 0) 
		    {
		      pmap[cg->proc_id[c]]=1;
		      dataPtr->next=(INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
		      dataPtr=dataPtr->next;
          dataPtr->intDataSize=4;
          dataPtr->realDataSize=4;
          dataPtr->realData=(double *) malloc(sizeof(double)*dataPtr->realDataSize);
          dataPtr->intData=(int *) malloc(sizeof(int)*dataPtr->intDataSize);
		      dataPtr->intData[0]=cg->proc_id[c];
		      dataPtr->intData[1]=cg->local_id[c];
	        dataPtr->intData[2]=itm;
	        dataPtr->intData[3]=0;
		      nsearch+=1;

          for(kk=0;kk<dataPtr->realDataSize-1;kk++)
            dataPtr->realData[kk]=xtm[kk];

          dataPtr->realData[dataPtr->realDataSize-1]=vol;

		      dataPtr->next=NULL;
		    }
		}

    // loop over all nodes
    for(j=0;j<cg->dims[3*c]+1;j++)
      for(k=0;k<cg->dims[3*c+1]+1;k++)
        for(l=0;l<cg->dims[3*c+2]+1;l++)
        {
          // Q[nq,nZ+1+2*nf,nY+1+2*nf,nX+1+2*nf]--> C++ Node storage
          // indexing starts from number of cells onwards
          itm = (cg->dims[3*c+1]+1+2*cg->nf)*(cg->dims[3*c]+1+2*cg->nf)*(l+cg->nf)
              + (cg->dims[3*c]+1+2*cg->nf)*(k+cg->nf) + (j+cg->nf)
              + cell_count;

          xtm[0] = cg->xlo[3*c]   + j*cg->dx[3*c];
          xtm[1] = cg->xlo[3*c+1] + k*cg->dx[3*c+1];
          xtm[2] = cg->xlo[3*c+2] + l*cg->dx[3*c+2];

          iflag=0;

          {
            //if (intersectCount==1 && myid==0) fprintf(fp,"%lf %lf %lf\n",xtm[0],xtm[1],xtm[2]);
            for(jj=0;jj<3;jj++) xd[jj]=0;
            for(jj=0;jj<3;jj++)
              for(kk=0;kk<3;kk++)
                xd[jj]+=(xtm[kk]-obb->xc[kk])*obb->vec[jj][kk];

            if (fabs(xd[0]) <= obb->dxc[0] &&
                fabs(xd[1]) <= obb->dxc[1] &&
                fabs(xd[2]) <= obb->dxc[2])
            {
              iflag++;
            }
          }

          if (iflag > 0)
          {
            pmap[cg->proc_id[c]]=1;
            dataPtr->next=(INTEGERLIST2 *) malloc(sizeof(INTEGERLIST2));
            dataPtr=dataPtr->next;
            dataPtr->intDataSize=4;
            dataPtr->realDataSize=4;
            dataPtr->realData=(double *) malloc(sizeof(double)*dataPtr->realDataSize);
            dataPtr->intData=(int *) malloc(sizeof(int)*dataPtr->intDataSize);
            dataPtr->intData[0]=cg->proc_id[c];
            dataPtr->intData[1]=cg->local_id[c];
            dataPtr->intData[2]=itm;
            dataPtr->intData[3]=0;
            nsearch+=1;

            for(kk=0;kk<dataPtr->realDataSize-1;kk++)
              dataPtr->realData[kk]=xtm[kk];

            dataPtr->realData[dataPtr->realDataSize-1]=vol;

            dataPtr->next=NULL;
          }
        }

	  TIOGA_FREE(xtm);
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
  if (donorId) TIOGA_FREE(donorId);
  if (tagsearch) TIOGA_FREE(tagsearch);
  if (isearch) TIOGA_FREE(isearch);
  if (rst) TIOGA_FREE(rst);
  if (res_search) TIOGA_FREE(res_search);
  if (xsearch) TIOGA_FREE(xsearch);
  //
  donorId=(int *)malloc(sizeof(int)*nsearch);
  tagsearch=(int *)malloc(sizeof(int)*nsearch);
  isearch=(int *)malloc(3*sizeof(int)*nsearch);
  res_search=(double *)malloc(sizeof(double)*nsearch);
  xsearch=(double *)malloc(sizeof(double)*3*nsearch);
  rst=(double *) malloc(sizeof(double)*3*nsearch);
  //
  dataPtr=head->next;
  k=l=m=n=0;
  while(dataPtr!=NULL)
  {
    for(j=0;j<dataPtr->intDataSize-1;j++)
      isearch[m++]=dataPtr->intData[j];
    tagsearch[k++]=dataPtr->intData[dataPtr->intDataSize-1];

    for(j=0;j<dataPtr->realDataSize-1;j++)
      xsearch[n++]=dataPtr->realData[j];
    res_search[l++]=dataPtr->realData[dataPtr->realDataSize-1];

    dataPtr=dataPtr->next;
  }

  deallocateLinkList3(head);
  //fclose(fp);
  TIOGA_FREE(obcart);
  TIOGA_FREE(pmap);
  TIOGA_FREE(sndMap);
  TIOGA_FREE(rcvMap);
}

  
