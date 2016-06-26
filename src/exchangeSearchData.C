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

void tioga::exchangeSearchData(void)
{
  // Get the processor map for sending and receiving
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  // Create packets to send and receive points, and initialize them to zero
  PACKET *sndPack = (PACKET *)malloc(sizeof(PACKET)*nsend);
  PACKET *rcvPack = (PACKET *)malloc(sizeof(PACKET)*nrecv);

  for (int i = 0; i < nsend; i++)
  {
    sndPack[i].nints = 0;
    sndPack[i].nreals = 0;
    sndPack[i].intData = NULL;
    sndPack[i].realData = NULL;
  }

  for (int i = 0; i < nrecv; i++)
  {
    rcvPack[i].nints = 0;
    rcvPack[i].nreals = 0;
    rcvPack[i].intData = NULL;
    rcvPack[i].realData = NULL;
  }

  // Find all mesh nodes on this grid which lie within the bounding box of each
  // other grid
  for (int k = 0; k < nsend; k++)
    mb->getQueryPoints(&obblist[k],
		       &sndPack[k].nints,&sndPack[k].intData,
		       &sndPack[k].nreals,&sndPack[k].realData);

  // Exchange the data
  pc->sendRecvPackets(sndPack,rcvPack);

  // now assort the data into the search list arrays
  mb->nsearch = 0;
  for (int k = 0; k < nrecv; k++)
    mb->nsearch += rcvPack[k].nints;

  // If these were already allocated, free & re-allocate
  free(mb->xsearch);
  free(mb->isearch);
  free(mb->donorId);

  mb->xsearch = (double *)malloc(sizeof(double)*3*mb->nsearch);
  mb->isearch = (int *)malloc(2*sizeof(int)*mb->nsearch);
  mb->donorId = (int *)malloc(sizeof(int)*mb->nsearch);

  // now fill the query point arrays
  int icount = 0;
  int dcount = 0;
  for (int k = 0; k < nrecv; k++)
  {
    int l = 0;
    for (int j = 0; j < rcvPack[k].nints; j++)
    {
      mb->isearch[icount++] = k;
      mb->isearch[icount++] = rcvPack[k].intData[j];
      mb->xsearch[dcount++] = rcvPack[k].realData[l++];
      mb->xsearch[dcount++] = rcvPack[k].realData[l++];
      mb->xsearch[dcount++] = rcvPack[k].realData[l++];
    }
  }

  for (int i = 0; i < nsend; i++)
  {
    free(sndPack[i].intData);
    free(sndPack[i].realData);
  }

  for (int i = 0; i < nrecv; i++)
  {
    free(rcvPack[i].intData);
    free(rcvPack[i].realData);
  }

  free(sndPack);
  free(rcvPack);
}
  
//
// routine for extra query points
// have to unify both routines here 
// FIX later ...
//
void tioga::exchangePointSearchData(void)
{
  int i,j,k,l,m;
  int icount,dcount;
  int nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int *sndMap;
  int *rcvMap;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  for(i=0;i<nsend;i++)     
    {
      sndPack[i].nints=sndPack[i].nreals=0;
      sndPack[i].intData=NULL;
      sndPack[i].realData=NULL;
    }
  //
  for(i=0;i<nrecv;i++)     
    {
      rcvPack[i].nints=rcvPack[i].nreals=0;
      rcvPack[i].intData=NULL;
      rcvPack[i].realData=NULL;
    }
  //
  // now get data for each packet
  //
  for(k=0;k<nsend;k++)
    mb->getExtraQueryPoints(&obblist[k],
			    &sndPack[k].nints,&sndPack[k].intData,
			    &sndPack[k].nreals,&sndPack[k].realData);
  MPI_Barrier(scomm);
  //if (myid==0) printf("AAAAA\n");
  //
  // exchange the data
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // now assort the data into the search
  // list arrays
  //
  mb->nsearch=0;
  for(k=0;k<nrecv;k++)
    mb->nsearch+=rcvPack[k].nints;
  //
  // if these were already allocated
  // get rid of them
  //
  if (mb->xsearch) free(mb->xsearch);
  if (mb->isearch) free(mb->isearch);
  if (mb->donorId) free(mb->donorId);
  if (mb->rst) free(mb->rst);
  //
  // allocate query point storage
  //
  mb->xsearch=(double *)malloc(sizeof(double)*3*mb->nsearch);
  mb->isearch=(int *)malloc(2*sizeof(int)*mb->nsearch);
  mb->donorId=(int *)malloc(sizeof(int)*mb->nsearch);
  mb->rst=(double *) malloc(sizeof(double)*3*mb->nsearch);
  //
  // now fill the query point arrays
  //
  icount=0;
  dcount=0;
  for(k=0;k<nrecv;k++)
    {
      l=0;
      for(j=0;j<rcvPack[k].nints;j++)
	{
	  mb->isearch[icount++]=k;
	  mb->isearch[icount++]=rcvPack[k].intData[j];	  
	  mb->xsearch[dcount++]=rcvPack[k].realData[l++];
	  mb->xsearch[dcount++]=rcvPack[k].realData[l++];
	  mb->xsearch[dcount++]=rcvPack[k].realData[l++];
	}
    }
  for(i=0;i<nsend;i++)
    {
      if (sndPack[i].nints > 0) free(sndPack[i].intData);
      if (sndPack[i].nreals >0) free(sndPack[i].realData);
    }
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nints > 0) free(rcvPack[i].intData);
      if (rcvPack[i].nreals >0) free(rcvPack[i].realData);
    }
  free(sndPack);
  free(rcvPack);
  //printf("%d %d\n",myid,mb->nsearch);
}
  
