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
extern "C" 
{
  void insertInList(DONORLIST **donorList,DONORLIST *temp1);
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);  
  int checkHoleMap(double *x,int *nx,int *sam,double *extents);
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}

void MeshBlock::getDonorPacket(PACKET *sndPack, int nsend)
{
  int i,k;
  int *icount;
  int *dcount;

  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nsend);
  //
  // count numbers to send first
  //
  for(i=0;i<nsearch;i++)
    {
      if (donorId[i] > -1)
        {
          k=isearch[2*i];
          sndPack[k].nints+=3;
          sndPack[k].nreals++;
        }
    }
  for(k=0;k<nsend;k++)
    {
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
    }

  for(i=0;i<nsend;i++) {icount[i]=dcount[i]=0;};
  //
  for(i=0;i<nsearch;i++)
    {
      if (donorId[i] > -1)
        {
          k=isearch[2*i];
          sndPack[k].intData[icount[k]++]=meshtag;               // mesh tag
          sndPack[k].intData[icount[k]++]=isearch[2*i+1];        // point id
          sndPack[k].intData[icount[k]++]=i;                         // point id on the donor side
          sndPack[k].realData[dcount[k]++]=cellRes[donorId[i]];  // donor resolution
        }
    }
  free(icount);
  free(dcount);
}
void MeshBlock::initializeDonorList(void)
{
  int i;
  if (donorList) 
    {
      for(i=0;i<donorListLength;i++){
	deallocateLinkList(donorList[i]);
	//printf("\n\t Deallocate (nnodes) i: %d %d ",nnodes,i);
       }
      free(donorList);
    }

  donorListLength = nnodes;
  donorList=(DONORLIST **)malloc(sizeof(DONORLIST *)*donorListLength);
  for(i=0;i<donorListLength;i++)
    donorList[i]=NULL;
}

void MeshBlock::insertAndSort(int pointid,int senderid,int meshtagdonor, int remoteid,
			      double donorRes)
{
  DONORLIST *temp1;
  temp1=(DONORLIST *)malloc(sizeof(DONORLIST));
  temp1->donorData[0]=senderid;
  temp1->donorData[1]=meshtagdonor;
  temp1->donorData[2]=remoteid;
  temp1->donorRes=donorRes;
  insertInList(&donorList[pointid],temp1);
}

void MeshBlock::processDonors(HOLEMAP *holemap, int nmesh, int **donorRecords,double **receptorResolution,
			      int *nrecords)
{
  //
  // first mark hole points
  //
  int *iflag = (int *)malloc(sizeof(int)*nmesh);

  for (int i = 0; i < nnodes; i++)
  {
    iblank[i] = NORMAL;
    if (donorList[i]==NULL)
    {
      for (int j = 0; j < nmesh; j++)
        if (j != (meshtag-BASE) && holemap[j].existWall)
        {
          if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents))
          {
            iblank[i] = HOLE;
            break;
          }
        }
    }
    else
    {
      DONORLIST *temp = donorList[i];

      for (int j = 0; j < nmesh; j++) iflag[j] = 0;

      while (temp != NULL)
      {
        int meshtagdonor = temp->donorData[1]-BASE;
        iflag[meshtagdonor] = 1;
        temp = temp->next;
      }
      for (int j = 0; j < nmesh; j++)
      {
        if (j != (meshtag-BASE) && holemap[j].existWall)
        {
          if (!iflag[j])
            if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents))
            {
              iblank[i] = HOLE;
              break;
            }
        }
      }
    }
  }

  for (int i = 0; i < nwbc; i++)
  {
   if (iblank[wbcnode[i]-BASE] == HOLE) {
     printf("--------------------------------------------------------------------\n");
     printf("Alarm from process %d : wall node is being tagged as a hole %d %p\n",myid,wbcnode[i]-BASE,
        donorList[wbcnode[i]-BASE]);
     int ii=wbcnode[i]-BASE;
     printf("xloc=%e %e %e\n",x[3*ii],x[3*ii+1],x[3*ii+2]);
     printf("Computations will continue, but may suffer from accuracy problems\n");
     printf("Please recheck positions of your grids\n");
     printf("--------------------------------------------------------------------\n");
    }
  }
  //
  // mark mandatory fringes as neighbors (up to nfringe depth)
  // of hole points 
  //
  int *mtag = (int *)malloc(sizeof(int)*nnodes);
  int *mtag1 = (int *)malloc(sizeof(int)*nnodes);
 
  for (int i = 0; i < nnodes; i++)
  {
    mtag[i] = mtag1[i] = HOLE;
    if (iblank[i] == HOLE)
      mtag[i] = mtag1[i] = NORMAL;
  }

  for (int iter = 0; iter < nfringe; iter++)
  {
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for (int i = 0; i < nc[n]; i++)
      {
        for (int m = 0; m < nvert; m++)
        {
          if (mtag[(vconn[n][nvert*i+m]-BASE)] == NORMAL)
          {
            for (int mm = 0; mm < nvert; mm++)
              if (m != mm && mtag[vconn[n][nvert*i+mm]-BASE] != NORMAL)
                mtag1[vconn[n][nvert*i+mm]-BASE] = NORMAL;
          }
        }
      }
    }

    for (int i = 0; i < nnodes; i++) mtag[i] = mtag1[i];
  }

  for (int i = 0; i < nnodes; i++)
    if (mtag1[i] && iblank[i])
      nodeRes[i]=BIGVALUE;

  free(mtag);
  free(mtag1);

  //
  // now find fringes
  //
  *nrecords=0;
  for (int i = 0; i < nnodes; i++)
  {
    if (donorList[i] != NULL && iblank[i] != HOLE)
    {
      DONORLIST *temp = donorList[i];
      while (temp != NULL)
      {
        if (temp->donorRes < nodeRes[i])
        {
          iblank[i] = FRINGE;
          (*nrecords)++;
          break;
        }
        temp = temp->next;
      }
    }
  }

  //
  // set the records to send back to the donor
  // process
  //
  (*donorRecords)=(int *)malloc(sizeof(int)*2*(*nrecords));
  (*receptorResolution)=(double *)malloc(sizeof(double)*(*nrecords));

  int m = 0;
  int k = 0;
  for (int i = 0; i < nnodes; i++)
  {
    if (iblank[i] == FRINGE)
    {
      DONORLIST *temp = donorList[i];
      (*receptorResolution)[k++] = nodeRes[i];
      (*donorRecords)[m++] = temp->donorData[0];
      (*donorRecords)[m++] = temp->donorData[2];
    }
  }

  //
  // release local memory
  //
  free(iflag);
}

void MeshBlock::initializeInterpList(int ninterp_input)
{
  int i;
  if (interpList) {
    //for(i=0;i<ninterp;i++)
    for(i=0;i<interpListSize;i++)
      {
	if (interpList[i].inode) free(interpList[i].inode);
	if (interpList[i].weights) free(interpList[i].weights);
      }
    free(interpList);
  }
  ninterp=ninterp_input;   
  interpListSize=ninterp_input;
  interpList=(INTERPLIST *)malloc(sizeof(INTERPLIST)*interpListSize);
  for(i=0;i<interpListSize;i++) {
   interpList[i].inode=NULL;
   interpList[i].weights=NULL;
  }
  if (cancelList) deallocateLinkList2(cancelList);
  cancelList=NULL;
  ncancel=0;
  if (interp2donor) free(interp2donor);
  interp2donor=(int *)malloc(sizeof(int)*nsearch);
  for(i=0;i<nsearch;i++) interp2donor[i]=-1;
    
}
		
void MeshBlock::findInterpData(int *recid,int irecord,double receptorRes)
{
  int i,j,i3,m,n;
  int nvert;
  int isum;
  int procid,pointid;
  double xv[8][3];
  double xp[3];
  double frac[8];
  int inode[8];
  int acceptFlag;
  INTEGERLIST *clist;
  //
  procid=isearch[2*irecord];
  pointid=isearch[2*irecord+1];
  i3=3*irecord;
  xp[0]=xsearch[i3];
  xp[1]=xsearch[i3+1];
  xp[2]=xsearch[i3+2]; 
  //
  isum=0;
  for(n=0;n<ntypes;n++)
    {
      isum+=nc[n];
      if (donorId[irecord] < isum) 
	{
	  i=donorId[irecord]-(isum-nc[n]);
          break;
	}
    }
  nvert=nv[n];
  acceptFlag=1;
  for(m=0;m<nvert;m++)
    {
      inode[m]=vconn[n][nvert*i+m]-BASE;
      i3=3*inode[m];
      if (iblank[inode[m]] <=0)
        {
         acceptFlag=0;
        }
      for(j=0;j<3;j++)
        xv[m][j]=x[i3+j];
    }
  //
  if (acceptFlag==0 && receptorRes!=BIGVALUE) return;
  if (receptorRes==BIGVALUE)
    {
      clist=cancelList;
      //
      // go to the end of the list 
      //
      if (clist !=NULL) while(clist->next !=NULL) clist=clist->next;
      //
      for(m=0;m<nvert;m++)
	{
          inode[m]=vconn[n][nvert*i+m]-BASE;
	  if (iblank[inode[m]]==-1 && nodeRes[inode[m]]!=BIGVALUE) 
	    {
	      iblank[inode[m]]=1;
	      if (clist == NULL) 
		{
		  clist=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
		  clist->inode=inode[m];
		  clist->next=NULL;
                  cancelList=clist;
		}
	      else
		{
		  clist->next=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
		  clist->next->inode=inode[m];
		  clist->next->next=NULL;
		  clist=clist->next;
		}
	      ncancel++;
	    }
	}
    }

  //  Find shape-function interpolation weights for query point
  // if (!ihighGlobal) // Unnecessary - will be done later if high-order grids present
  computeNodalWeights(xv,xp,frac,nvert);

  interp2donor[irecord]=*recid;
  interpList[*recid].cancel=0;
  interpList[*recid].nweights=nvert;
  interpList[*recid].receptorInfo[0]=procid;
  interpList[*recid].receptorInfo[1]=pointid;
  interpList[*recid].inode=(int *)malloc(sizeof(int)*nvert);
  interpList[*recid].weights=(double *)malloc(sizeof(double)*nvert);
  for(m=0;m<nvert;m++)
    {
      interpList[*recid].inode[m]=inode[m];
      interpList[*recid].weights[m]=frac[m];
    }
  (*recid)++;
}

void MeshBlock::set_ninterp(int ninterp_input)
{
  ninterp=ninterp_input;
}
  
void MeshBlock::getCancellationData(int *nrecords,int **intData)
{
  int i;
  int inode;
  INTEGERLIST *clist;   
  *nrecords=ncancel;
  if (ncancel > 0) 
    {
      (*intData)=(int *)malloc(sizeof(int)*(*nrecords)*2);
      i=0;
      for(clist=cancelList;clist!=NULL;clist=clist->next) 
	{
	  inode=clist->inode;
	  (*intData)[i++]=donorList[inode]->donorData[0];
	  (*intData)[i++]=donorList[inode]->donorData[2];
	}
    }
}

void MeshBlock::cancelDonor(int irecord)
{
  int iptr;
  iptr=interp2donor[irecord];
  if (iptr > -1) interpList[iptr].cancel=1;
}

void MeshBlock::getInterpData(int *nrecords, int **intData)
{
  int i,k;
  //
  *nrecords=0;
  for(i=0;i<ninterp;i++)
    if (!interpList[i].cancel) (*nrecords)++;
  //
  (*intData)=(int *)malloc(sizeof(int)*2*(*nrecords));
  for(i=0,k=0;i<ninterp;i++)
    if (!interpList[i].cancel) {
       (*intData)[k++]=interpList[i].receptorInfo[0];
       (*intData)[k++]=interpList[i].receptorInfo[1];
    }
}

void MeshBlock::clearIblanks(void)
{
  int i;
  for(i=0;i<nnodes;i++)
     if (iblank[i] < 0) iblank[i]=1;
}

void MeshBlock::setIblanks(int inode)
{
  iblank[inode] = FRINGE;
}
