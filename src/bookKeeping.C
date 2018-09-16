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
#include "codetypes.h"
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
  TIOGA_FREE(icount);
  TIOGA_FREE(dcount);
}

void MeshBlock::getMBDonorPktSizes
(
  std::vector<int>& nints,
  std::vector<int>& nreals
)
{
  for(int i=0; i< nsearch; i++) {
    if (donorId[i] > -1) {
      int ii = isearch[3*i];
      nints[ii] += 4;
      nreals[ii]+=2;
    }
  }
}

void MeshBlock::getMBDonorPackets
(
  std::vector<int>& ixOffset,
  std::vector<int>& rxOffset,
  PACKET* sndPack
)
{
  for(int i=0; i<nsearch; i++) {
    if (donorId[i] < 0) continue;

    int k = isearch[3*i];
    int& ix = ixOffset[k];
    int& rx = rxOffset[k];

    sndPack[k].intData[ix++] = meshtag;           // Unique mesh tag
    sndPack[k].intData[ix++] = isearch[3*i + 1];  // point ID
    sndPack[k].intData[ix++] = i;                 // point ID on donor side
    sndPack[k].intData[ix++] = isearch[3*i + 2];  // receptor block ID
    sndPack[k].realData[rx++]=cellRes[donorId[i]];
    sndPack[k].realData[rx++]=res_search[xtag[i]];
  }  
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
      TIOGA_FREE(donorList);
    }

  donorListLength = nnodes;
  donorList=(DONORLIST **)malloc(sizeof(DONORLIST *)*donorListLength);
  for(i=0;i<donorListLength;i++)
    donorList[i]=NULL;
}

void MeshBlock::insertAndSort(int pointid,int senderid,int meshtagdonor, int remoteid,
			      double donorRes,double receptorRes)
{
  DONORLIST *temp1;
  temp1=(DONORLIST *)malloc(sizeof(DONORLIST));
  temp1->donorData[0]=senderid;
  temp1->donorData[1]=meshtagdonor;
  temp1->donorData[2]=remoteid;
  temp1->donorRes=donorRes;
  temp1->receptorRes=receptorRes;
  insertInList(&donorList[pointid],temp1);
}

void MeshBlock::processDonors(HOLEMAP *holemap, int nmesh, int **donorRecords,double **receptorResolution,
			      int *nrecords)
{
  int i,j,k,m,n,mm,ii;
  int nvert;
  DONORLIST *temp;
  int *iflag; 
  int meshtagdonor;
  int *mtag,*mtag1;
  int iter;
  int verbose;
  //
  // first mark hole points
  //
  iflag=(int *)malloc(sizeof(int)*nmesh);
  //
  for(i=0;i<nnodes;i++)
    {
      iblank[i]=1;
      verbose=0;
      //if (meshtag==1 && myid==0 && i==76639) verbose=1;
      //if (meshtag==2 && i==240304 && myid==1) verbose=1;
      //if (meshtag==3 && i==241402 && myid==1) verbose=1;
      /*
      if (fabs(x[3*i]-1.68) < 1e-4 && 
          fabs(x[3*i+1]-0.88) < 1e-4 &&
          fabs(x[3*i+2]) < 1e-4 && meshtag==3 && myid==1) {
       verbose=1;
      }
      */

      if (verbose) TRACEI(i);
      if (verbose) TRACEI(iblank[i]);
      if (verbose) TRACED(nodeRes[i]);
      if (verbose) printf("%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);
      if (donorList[i]==NULL)
	{
          if (verbose) {
            printf("No donor found for %d\n",i);
          }
	  for(j=0;j<nmesh;j++)
           if (j!=(meshtag-BASE) && holemap[j].existWall) 
            {
	     if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents)) 
	      {
		iblank[i]=0;
		break;
	      }
           }
	}
      else
	{
	  temp=donorList[i];
	  for(j=0;j<nmesh;j++) iflag[j]=0;
	  while(temp!=NULL) 
	    {
	      meshtagdonor=temp->donorData[1]-BASE;
	      iflag[meshtagdonor]=1;
              if (verbose) {
               TRACEI(meshtagdonor);
	       TRACED(temp->donorRes);
               TRACEI(temp->donorData[0]);
               TRACEI(temp->donorData[1]);
               TRACEI(temp->donorData[2]);
              }
              nodeRes[i]=TIOGA_MAX(nodeRes[i],temp->receptorRes);
	      temp=temp->next;
	    }
	  for(j=0;j<nmesh;j++)
	    {
             if (j!=(meshtag-BASE) && holemap[j].existWall)
              {
	       if (!iflag[j])
		if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents))
		  {
		    iblank[i]=0;
		    break;
		  }
              }
	    }
	}
    }
  for(i=0;i<nwbc;i++)
   if (iblank[wbcnode[i]-BASE]==0) {
     printf("--------------------------------------------------------------------\n");
     printf("Alarm from process %d : wall node is being tagged as a hole %d %p\n",myid,wbcnode[i]-BASE,
        donorList[wbcnode[i]-BASE]);
     ii=wbcnode[i]-BASE;
     printf("xloc=%e %e %e\n",x[3*ii],x[3*ii+1],x[3*ii+2]);
     printf("Computations will continue, but may suffer from accuracy problems\n");
     printf("Please recheck positions of your grids\n");
     printf("--------------------------------------------------------------------\n");
    }
  //
  // mark mandatory fringes as neighbors (up to nfringe depth)
  // of hole points 
  //
  mtag=(int *)malloc(sizeof(int)*nnodes);
  mtag1=(int *)malloc(sizeof(int)*nnodes);
 
  for(i=0;i<nnodes;i++)
    {
     mtag[i]=mtag1[i]=0;
     if (iblank[i]==0) mtag[i]=mtag1[i]=1;
    }

 for(iter=0;iter<nfringe;iter++)
 {
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  for(m=0;m<nvert;m++)
	    {
	      if (mtag[(vconn[n][nvert*i+m]-BASE)]==1)
		{
		  for(mm=0;mm<nvert;mm++)
		    if (m!=mm && mtag[vconn[n][nvert*i+mm]-BASE] !=1) 
	                 mtag1[vconn[n][nvert*i+mm]-BASE] = 1;
		}
	    }
	}
    }
   for(i=0;i<nnodes;i++) mtag[i]=mtag1[i];
  }
  for(i=0;i<nnodes;i++)
     if (mtag1[i] && iblank[i]) nodeRes[i]=BIGVALUE;
  TIOGA_FREE(mtag);
  TIOGA_FREE(mtag1);
  //
  // now find fringes
  //
  *nrecords=0;
  for(i=0;i<nnodes;i++)
    {
      verbose=0;
      //if (meshtag==1 && myid==0 && i==76639) verbose=1;
      //if (meshtag==2 && i==240304 && myid==1) verbose=1;
      //if (meshtag==3 && i==241402 && myid==1) verbose=1;
      //if (meshtag==3 && i==34299) verbose=1;
      if (verbose) {
         TRACEI(i);
         TRACEI(iblank[i]);
         TRACED(nodeRes[i]);
      }
      if (donorList[i]!=NULL && iblank[i]!=0)
	{ 
	  temp=donorList[i];
	  while(temp!=NULL)
	    {
	      if (verbose) TRACED(temp->donorRes);
	      if (temp->donorRes < nodeRes[i])
		{
		  iblank[i]=-temp->donorData[1];
                  if (verbose) TRACEI(iblank[i]);
                  if (verbose) {
                  TRACEI(temp->donorData[0]);
                  TRACEI(temp->donorData[1]);
                  TRACEI(temp->donorData[2]);}
		  (*nrecords)++;
		  break;
		}
	      temp=temp->next;
	    }
	}
    }
  //
  // set the records to send back to the donor
  // process
  //
  (*donorRecords)=(int *)malloc(sizeof(int)*3*(*nrecords));
  (*receptorResolution)=(double *)malloc(sizeof(double)*(*nrecords));
  m=0;
  k=0;
  for(i=0;i<nnodes;i++)
    {
      verbose=0;
      //if (meshtag==1 && myid==0 && i==76639) verbose=1;
      //if (meshtag==3 && myid==1 && i==245609) verbose=1;
      //if (meshtag==2 && i==240304 && myid==1) verbose=1;
      //if (meshtag==3 && i==241402 && myid==1) verbose=1;
      //if (meshtag==3 && i==34299) verbose=1;
      if (iblank[i] < 0) 
	{
	  temp=donorList[i];
          while(temp!=NULL)
           { 
            if (temp->donorRes < nodeRes[i])
             {
              break;
             }
           }
          //if (temp->donorRes < 0) nodeRes[i]=BIGVALUE;
	  (*receptorResolution)[k++]=(resolutionScale > 1.0) ? -nodeRes[i]:nodeRes[i];
	  (*donorRecords)[m++]=temp->donorData[0];
	  (*donorRecords)[m++]=temp->donorData[2];
          (*donorRecords)[m++]=temp->donorData[1];
          if (verbose) {
            TRACEI(iblank[i]);           
            TRACEI(m);
            TRACEI((*donorRecords)[m-3]);
            TRACEI((*donorRecords)[m-1]);
            TRACEI((*donorRecords)[m-2]);
	    TRACED((*receptorResolution)[k-1]);
          }
	}
    }  	      
  //
  // release local memory
  //
  TIOGA_FREE(iflag);
}

void MeshBlock::initializeInterpList(int ninterp_input)
{
  int i;
  if (interpList) {
    //for(i=0;i<ninterp;i++)
    for(i=0;i<interpListSize;i++)
      {
	if (interpList[i].inode) TIOGA_FREE(interpList[i].inode);
	if (interpList[i].weights) TIOGA_FREE(interpList[i].weights);
      }
    TIOGA_FREE(interpList);
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
  if (interp2donor) TIOGA_FREE(interp2donor);
  interp2donor=(int *)malloc(sizeof(int)*nsearch);
  for(i=0;i<nsearch;i++) interp2donor[i]=-1;
    
}
		
void MeshBlock::findInterpData(int *recid,int irecord,double receptorRes2)
{
  int i,j,i3,m,n;
  int nvert;
  int isum;
  int procid,pointid,blockid;
  double xv[8][3];
  double xp[3];
  double frac[8];
  int inode[8];
  int acceptFlag; 
  double receptorRes;
  int verbose;
  int meshtagrecv;

  INTEGERLIST *clist;
  //
  verbose=0;
  //if (myid==3 && irecord==4878 && meshtag==2) verbose=1;
  //if (myid==63 && irecord==3224) verbose=1;
  //if (myid==1 && irecord==158192 && meshtag==1) verbose=1;
  receptorRes=fabs(receptorRes2);
  procid=isearch[3*irecord];
  pointid=isearch[3*irecord+1];
  blockid=isearch[3*irecord+2];
  meshtagrecv=tagsearch[irecord];
  if (verbose) {
      TRACEI(procid);
      TRACEI(pointid);
      TRACED(receptorRes);
      TRACEI(meshtagrecv);
  }
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
  if (verbose) TRACEI(donorId[irecord])
  if (verbose) TRACEI(n)
  if (verbose) TRACEI(nvert)
  if (verbose) TRACEI(i)
  if (verbose) TRACED(nodeRes[241291]);
  for(m=0;m<nvert;m++)
    {
      inode[m]=vconn[n][nvert*i+m]-BASE;
      if (verbose) TRACEI(inode[m]);
      if (verbose) TRACEI(iblank[inode[m]])
      if (verbose) TRACED(nodeRes[inode[m]])
      i3=3*inode[m];
      if (iblank[inode[m]] <=0 && receptorRes2 > 0.0)
        {
         if (nodeRes[inode[m]]==BIGVALUE) acceptFlag=0;
         if (abs(iblank[inode[m]])==meshtagrecv) acceptFlag=0;
        }
      for(j=0;j<3;j++)
        xv[m][j]=x[i3+j];
    }
  //
  if (verbose) TRACEI(acceptFlag);
  if (receptorRes==BIGVALUE && resolutionScale==1.0)
    {
      clist=cancelList;
      //
      // go to the end of the list 
      //
      if (clist !=NULL) while(clist->next !=NULL) clist=clist->next;
      //
      for(m=0;m<nvert;m++)
	{
          verbose=0;
          inode[m]=vconn[n][nvert*i+m]-BASE;
          //if (myid==763 && inode[m]==9515) verbose=1;
          //if (myid==1 && meshtag==2 && inode[m]==240304) verbose=1;
          if (verbose) TRACEI(inode[m]);
          if (verbose) TRACED(nodeRes[inode[m]]);
          if (verbose) {
              TRACEI(procid);
              TRACEI(pointid);
              TRACED(receptorRes);
              TRACEI(meshtagrecv);
              TRACEI(irecord);
              TRACEI(donorId[irecord]);
          }
	  if (iblank[inode[m]]<=0 && nodeRes[inode[m]]!=BIGVALUE) 
	    {
	      if (iblank[inode[m]] < 0) iblank[inode[m]]=1;
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
  //  
  computeNodalWeights(xv,xp,frac,nvert);
  //
  interp2donor[irecord]=*recid;
  interpList[*recid].cancel=0;
  interpList[*recid].nweights=nvert;
  interpList[*recid].receptorInfo[0]=procid;
  interpList[*recid].receptorInfo[1]=pointid;
  interpList[*recid].receptorInfo[2]=blockid;
  if (verbose) {
    TRACEI(*recid);
    TRACEI(interpList[*recid].receptorInfo[0]);
    TRACEI(interpList[*recid].receptorInfo[1]);
  }
  interpList[*recid].inode=(int *)malloc(sizeof(int)*(nvert+1));
  interpList[*recid].weights=(double *)malloc(sizeof(double)*(nvert+1));
  for(m=0;m<nvert;m++)
    {
      interpList[*recid].inode[m]=inode[m];
      interpList[*recid].weights[m]=frac[m];
      if ( frac[m] < -0.2 || frac[m] > 1.2) {
       TRACEI(myid);
       TRACEI(irecord);
       TRACEI(meshtag);
       TRACEI(donorId[irecord]);
       TRACED(frac[m]);
       int ierr;
       MPI_Abort(MPI_COMM_WORLD,ierr);
      }
    }
  interpList[*recid].inode[m]=donorId[irecord];
  interpList[*recid].weights[m]=0.0;
  if (acceptFlag==0 && receptorRes!=BIGVALUE) interpList[*recid].cancel=1;
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
      (*intData)=(int *)malloc(sizeof(int)*(*nrecords)*3);
      i=0;
      for(clist=cancelList;clist!=NULL;clist=clist->next) 
	{
	  inode=clist->inode;
          if (donorList[inode]!=NULL) {
      	    (*intData)[i++]=donorList[inode]->donorData[0];
	    (*intData)[i++]=donorList[inode]->donorData[2];
	    (*intData)[i++]=donorList[inode]->donorData[1];
         }
	}
      *nrecords=i/3;
    }
}

void MeshBlock::cancelDonor(int irecord)
{
  int iptr;
  iptr=interp2donor[irecord];
  if (iptr > -1) interpList[iptr].cancel=1;
}

void MeshBlock::resetCoincident(void)
{
  int i,iptr;
  int *ireset;

  ireset=(int *)malloc(sizeof(int)*nsearch);
  for(i=0;i<nsearch;i++) ireset[i]=1;
  
  for(i=0;i<nsearch;i++)
    {
      iptr=interp2donor[i];
      if (iptr > -1) {
        ireset[xtag[i]]=TIOGA_MIN(ireset[xtag[i]],interpList[iptr].cancel);
      }
    }	
  for(i=0;i<nsearch;i++)
    {
      iptr=interp2donor[i];
      if (iptr > -1) {
	if (interpList[iptr].cancel==1) {
	  interpList[iptr].cancel=ireset[xtag[i]];
	}
      }
    }
  TIOGA_FREE(ireset);
}
void MeshBlock::getInterpData(int *nrecords, int **intData)
{
  int i,k;
  //
  *nrecords=0;
  for(i=0;i<ninterp;i++)
    if (!interpList[i].cancel) (*nrecords)++;
  //
  (*intData)=(int *)malloc(sizeof(int)*3*(*nrecords));
  for(i=0,k=0;i<ninterp;i++)
    if (!interpList[i].cancel) {
       (*intData)[k++]=interpList[i].receptorInfo[0];
       (*intData)[k++]=interpList[i].receptorInfo[1];
       (*intData)[k++]=interpList[i].receptorInfo[2];
    }
}

void MeshBlock::clearIblanks(void)
{
  int i;
  for(i=0;i<nnodes;i++)
     if (iblank[i] < 0) iblank[i]=1;
  if (iblank_reduced) {
   for(i=0;i<nnodes;i++)
     if (iblank_reduced[i]==0) iblank[i]=0;
   TIOGA_FREE(iblank_reduced);
  }
}

void MeshBlock::getStats(int mstats[2])
{
  int i;
  mstats[0]=mstats[1]=0;
  for (i=0;i<nnodes;i++)
    {
      if (iblank[i]==0) mstats[0]++;
      if (iblank[i] < 0) mstats[1]++;
    }
}

void MeshBlock::setIblanks(int inode)
{
/*  if (fabs(nodeRes[inode]-BIGVALUE) < TOL)
    {
      iblank[inode]=-2;
    }
  else*/
//    {
   iblank[inode]=-1;
//    }
}

void MeshBlock::reduce_fringes(void)
{
  int *ibltmp;
  int m,n,nvert,i,flag,ncount,inode[8];
  int verbose,iter;
  INTEGERLIST *clist;
  //
  if (iblank_reduced) TIOGA_FREE(iblank_reduced);
  iblank_reduced=(int *)malloc(sizeof(int)*nnodes);
  for(int i=0;i< nnodes;i++) iblank_reduced[i]=iblank[i] > 0 ? iblank[i]:0;
  ibltmp=(int *)malloc(sizeof(int)*nnodes);
  //
  // make sure only partial iblank=1 + iblank=-1
  // cell nodes are tagged
  //
  //for(n=0;n<ntypes;n++)
  //  {
  //    nvert=nv[n];
  //    for(i=0;i<nc[n];i++)
  //	{
  //        verbose=0;
  //	  ncount=0;
  //	  for(m=0;m<nvert;m++)
  //	    {
  //	      inode[m]=vconn[n][nvert*i+m]-BASE;
  //	      ncount=ncount+(iblank[inode[m]] <=0);
  //	    }
  //	  if (ncount==nvert) {
  //	    for(m=0;m<nvert;m++) iblank_reduced[inode[m]]=0;
  //	  }
  //	}
  //  }
  //
  // now make sure neighbors of neighbors are tagged
  // to create a two point depth
  //
  for(i=0;i<nnodes;i++) ibltmp[i]=iblank_reduced[i];
  for(iter=0;iter < nfringe+1 ;iter++) 
  {
   for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
        {
          verbose=0;
          ncount=0;
          for(m=0;m<nvert;m++)
            {
              inode[m]=vconn[n][nvert*i+m]-BASE;
              ncount=ncount+(iblank_reduced[inode[m]] == 1 || iblank_reduced[inode[m]] < 0);
            }
	  if (ncount > 0 && ncount < nvert) {
	    for(m=0;m<nvert;m++)
	      if (ibltmp[inode[m]]==0 && iblank[inode[m]] < 0) ibltmp[inode[m]]=iblank[inode[m]];
	  }
        }
    }
    for(i=0;i<nnodes;i++) iblank_reduced[i]=ibltmp[i];
  }

  if (cancelList) deallocateLinkList2(cancelList);
  cancelList=NULL;
  ncancel=0;
  clist=cancelList;

  for(i=0;i<nnodes;i++) 
    {
      if (iblank[i] < 0 && iblank_reduced[i]==0) 
	{
	  if (clist == NULL) 
	    {
	      clist=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
	      clist->inode=i;
	      clist->next=NULL;
	      cancelList=clist;
	    }
	  else
	    {
	      clist->next=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
	      clist->next->inode=i;
	      clist->next->next=NULL;
	      clist=clist->next;
	    }
	  ncancel++;	
	}
    }

  // int norph=0;
  // for(n=0;n<ntypes;n++)
  //   {
  //     nvert=nv[n];
  //     for(i=0;i<nc[n];i++)
  //       {
  //         verbose=0;
  //         ncount=0;
  //         flag=0;
  //         for(m=0;m<nvert;m++)
  //           {
  //             inode[m]=vconn[n][nvert*i+m]-BASE;
  //             if (iblank_reduced[inode[m]]==1) flag=1;
  //             ncount=ncount+(iblank_reduced[inode[m]] == 1);
  //           }
  //         if (flag &&  ncount < nvert) {
  //           for(m=0;m<nvert;m++)
  //             if (iblank_reduced[inode[m]]==0) norph=norph+1;
  //         }
  //       }
  //   }
  // TRACEI(norph);
  TIOGA_FREE(ibltmp);
}





























