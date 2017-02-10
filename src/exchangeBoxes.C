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
#include "tioga.h"
extern "C"{
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
			double vB[3][3],double xB[3],double dxB[3]);			   
}
void tioga::exchangeBoxes(void)
{
  int i,j,k,m;
  int *alltags;
  int *sndMap;
  int *rcvMap;
  int *blockcount;
  int *displs;
  int nsend;
  int nrecv;
  int ntotalblks;
  int overlap_present;
  PACKET *sndPack,*rcvPack;
  //
  MPI_Allreduce(&nblocks,&ntotalblks,1,MPI_INT,MPI_SUM,scomm);
  //
  alltags=(int *)malloc(sizeof(int)*ntotalblks);
  pid=(int *)malloc(sizeof(int)*ntotalblks);
  blockcount=(int *)malloc(sizeof(int)*numprocs);
  displs=(int *)malloc(sizeof(int)*(numprocs+1));
  cflag=(int *)malloc(sizeof(int)*numprocs);
  //
  MPI_Allgather(&nblocks, 1, MPI_INT, blockcount,1,MPI_INT,scomm);
  m=0;
  for(i=0;i<numprocs;i++)
    {
      for(j=0;j<blockcount[i];j++)
	pid[m++]=i;
      cflag[i]=0;
    }
  //
  displs[0]=0;
  for(i=1;i<numprocs+1;i++)
    displs[i]=displs[i-1]+blockcount[i-1];
  //
  MPI_Allgatherv(&mytags, nblocks, MPI_INT, alltags,blockcount,displs,
		 MPI_INT,scomm);
  //
  // count number of other processors to communicate to
  // in overset grid scenario, usually you do not communicate
  // to mesh blocks that carry the same mesh tag (i.e. you don't
  // talk to your sister partitions)
  //
  nsend=nrecv=0;
  for(i=0;i<ntotalblks;i++) 
    {
      for(iblk=0;iblk<nblocks;iblk++)
	{
	  if (alltags[i] != mytags[iblk]) 
	    {
	      if (cflag[pid[i]]==0) nsend++;
	      cflag[pid[i]]=1;
	    }
	}
    }
  //
  // In general we communicate forward
  // and backward, separate lists are maintained for
  // flexibility
  //
  nrecv=nsend;
  sndMap=(int *)malloc(sizeof(int)*nsend);
  rcvMap=(int *)malloc(sizeof(int)*nrecv);
  //
  for(i=0,m=0;i<numprocs;i++)
    {
      if (cflag[i]==1)
	{
	  sndMap[m]=rcvMap[m]=i;
	  m++;
	}
    }
  //
  pc->setMap(nsend,nrecv,sndMap,rcvMap);
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  // 
  //  
  for(k=0;k<nsend;k++)
    {
      sndPack[k].nints=0;
      //
      for(iblk=0;iblk<nblocks;iblk++)
	for(d=displs[sndMap[k]];d<displs[sndMap[k]+1];d++)
	  {
	    if (alltags[d]!=mytags[iblk]) sndPack[k].nints++;
	  }
      //
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      m=0;
      sndPack[k].nreals=15*sndPack[k].nints;
      sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
      m=im=0;
      for(iblk=0;iblk<nblocks;iblk++)
	{
	  for(d=displs[sndMap[k]];d<displs[sndMap[k]+1];d++)
	    if (alltags[d]!=mytags[iblk]) 
	      {
		sndPack[k].intData[im++]=d-displs[sndMap[k]];
		for(i=0;i<3;i++)
		  for(j=0;j<3;j++)
		    sndPack[k].realData[m++]=mb[iblk]->obb->vec[i][j];
		for(i=0;i<3;i++)
		  sndPack[k].realData[m++]=mb[iblk]->obb->xc[i];
		for(i=0;i<3;i++)
		  sndPack[k].realData[m++]=mb[iblk]->obb->dxc[i];
	      }
	}
    }
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  if (obblist) free(obblist);
  nobb=0;
  for(k=0;k<nrecv;k++) nobb+=rcvPack[k].nints;
  //
  obblist=(OBB *) malloc(sizeof(OBB)*nobb);
  bid=(int *)malloc(sizeof(int)*nobb);
  //
  d=0;
  for(k=0;k<nrecv;k++)
    {
      m=im=0;
      for(n=0;n<rcvPack[k].nints;n++)
	{
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      obblist[d].vec[i][j]=rcvPack[k].realData[m++];
	  for(i=0;i<3;i++)
	    obblist[d].xc[i]=rcvPack[k].realData[m++];
	  for(i=0;i<3;i++)
	    obblist[d].dxc[i]=rcvPack[k].realData[m++];
	  bid[d++]=rcvPack[k].intData[im++];
	}
    }
  //
  m=0;
  for(k=0;k<nobb;k++)
    {
      iblk=bid[k];
      if ( obbIntersectCheck(mb[iblk]->obb->vec,mb[iblk]->obb->xc,mb[iblk]->obb->dxc,
			     obblist[k].vec,obblist[k].xc,obblist[k].dxc) ||
	   obbIntersectCheck(obblist[k].vec,obblist[k].xc,obblist[k].dxc,
			     mb[iblk]->obb->vec,mb[iblk]->obb->xc,mb[iblk]->obb->dxc)) 
	{
	  rcvMap[m]=sndMap[k];
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      obblist[m].vec[i][j]=obblist[k].vec[i][j];
	  for(i=0;i<3;i++)
	    obblist[m].xc[i]=obblist[k].xc[i];
	  for(i=0;i<3;i++)
	    obblist[m].dxc[i]=obblist[k].dxc[i];
	  bid[m]=bid[k];
	  m++;
	}
    }
  nsend=nrecv=m;
  for(i=0;i<nsend;i++) sndMap[i]=rcvMap[i];
  //printf("%d %d %d\n",myid,nsend,nrecv);
  //

  // clear packets before nsend and nrecv are modified in pc->setMap
  pc->clearPackets(sndPack,rcvPack);
  pc->setMap(nsend,nrecv,sndMap,rcvMap);  
  for(k=0;k<nsend;k++) {
   sndPack[k].nints=0;
   sndPack[k].nreals=6;
   sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
   mb->getReducedOBB(&obblist[k],sndPack[k].realData);
  }
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for(k=0;k<nrecv;k++) {
   for(j=0;j<3;j++) obblist[k].xc[j]=rcvPack[k].realData[j];
   for(j=0;j<3;j++) obblist[k].dxc[j]=rcvPack[k].realData[j+3];
   }  
  pc->clearPackets(sndPack,rcvPack);
  //
  // Free local memory
  //
  free(alltags);
  free(sndMap);
  free(rcvMap);
  free(sndPack);
  free(rcvPack);
}
