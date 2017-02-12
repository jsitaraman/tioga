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

#include <algorithm>
#include <vector>

#include "tioga.h"

void tioga::exchangeSearchData(void)
{
  int i;
  int nsend, nrecv;
  PACKET *sndPack, *rcvPack;
  int* sndMap;
  int* rcvMap;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack = (PACKET*)malloc(sizeof(PACKET) * nsend);
  rcvPack = (PACKET*)malloc(sizeof(PACKET) * nrecv);
  //
  for (i = 0; i < nsend; i++) {
    sndPack[i].nints = sndPack[i].nreals = 0;
    sndPack[i].intData = NULL;
    sndPack[i].realData = NULL;
  }
  //
  for (i = 0; i < nrecv; i++) {
    rcvPack[i].nints = rcvPack[i].nreals = 0;
    rcvPack[i].intData = NULL;
    rcvPack[i].realData = NULL;
  }

  // Process each intersection pair and determine the total data that needs to
  // be sent
  int nobb = obblist.size();
  std::vector<int> nintsSend(nobb);
  std::vector<int> nrealsSend(nobb);
  int** int_data = (int**)malloc(sizeof(int*) * nobb);
  double** real_data = (double**)malloc(sizeof(double*) * nobb);

  // Determine the ints and reals that need to be exchanged for each intersected
  // block
  for (int i = 0; i < nobb; i++) {
    int ib = obblist[i].iblk_local;
    int k = obblist[i].comm_idx;

    auto& mb = mblocks[ib];
    mb->getQueryPoints(
      &obblist[i], &nintsSend[i], &int_data[i], &nrealsSend[i], &real_data[i]);
    sndPack[k].nints += 2;
  }

  // Prepare data package to exchage array sizes
  for (int k = 0; k < nsend; k++) {
    sndPack[k].nreals = 0;
    sndPack[k].realData = NULL;
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
  }

  // Populate data package with array dimensions to be exchanged
  std::vector<int> ixOffset(nsend, 0);
  for (int i = 0; i < nobb; i++) {
    int k = obblist[i].comm_idx;
    int ioff = ixOffset[k];

    sndPack[k].intData[ioff] = nintsSend[i];
    sndPack[k].intData[ioff + 1] = nrealsSend[i];
    ixOffset[k] += 2;
  }
  pc->sendRecvPackets(sndPack, rcvPack);

  std::fill(ixOffset.begin(), ixOffset.end(), 0);
  std::vector<int> nintsRecv(nobb);
  std::vector<int> nrealsRecv(nobb);

  // Save off nints and nreals from the other proc in recv arrays
  for (int i = 0; i < nobb; i++) {
    int k = obblist[i].comm_idx;
    int ioff = ixOffset[k];

    nintsRecv[i] = rcvPack[k].intData[ioff];
    nrealsRecv[i] = rcvPack[k].intData[ioff + 1];
    ixOffset[k] += 2;
  }
  pc->clearPackets(sndPack, rcvPack);

  // Prepare packets to send the actual node indices and coordinate data
  //
  // Estimate array sizes
  for (int i = 0; i < nobb; i++) {
    int k = obblist[i].comm_idx;

    sndPack[k].nints += nintsSend[i];
    sndPack[k].nreals += nrealsSend[i];
  }

  // Allocate memory
  for (int k = 0; k < nsend; k++) {
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
    sndPack[k].realData = (double*)malloc(sizeof(double) * sndPack[k].nreals);
  }

  // Index of the next data entry point in the packet arrays
  std::vector<int> rxOffset(nsend, 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);

  // Populate packets with data
  for (int i = 0; i < nobb; i++) {
    int k = obblist[i].comm_idx;
    int ioff = ixOffset[k];
    int roff = rxOffset[k];

    for (int j = 0; j < nintsSend[i]; j++)
      sndPack[k].intData[ioff + j] = int_data[i][j];

    for (int j = 0; j < nrealsSend[j]; j++)
      sndPack[k].intData[roff + j] = real_data[i][j];

    ixOffset[k] += nintsSend[i];
    rxOffset[k] += nrealsSend[i];
  }

  pc->sendRecvPackets(sndPack, rcvPack);

  // Reset MeshBlock data structures
  for (auto& mb : mblocks) {
    mb->nsearch = 0;

    if (mb->xsearch)
      free(mb->xsearch);
    if (mb->isearch)
      free(mb->isearch);
    if (mb->donorId)
      free(mb->donorId);
  }

  // Calculate query point sizes for each mesh block
  for (int i = 0; i < nobb; i++) {
    int ib = obblist[i].iblk_local;
    auto& mb = mblocks[ib];

    mb->nsearch += nintsRecv[i];
  }

  for (auto& mb : mblocks) {
    mb->xsearch = (double*)malloc(sizeof(double) * 3 * mb->nsearch);
    mb->isearch = (int*)malloc(2 * sizeof(int) * mb->nsearch);
    mb->donorId = (int*)malloc(sizeof(int) * mb->nsearch);
  }

  // Finally populate data into individual mesh block arrays
  std::fill(ixOffset.begin(), ixOffset.end(), 0);
  std::vector<int> ibOffset(nblocks, 0);

  for (int i = 0; i < nobb; i++) {
    int ib = obblist[i].iblk_local;
    int k = obblist[i].comm_idx;
    int ioff = ixOffset[k];
    int iboff = ibOffset[ib];

    auto& mb = mblocks[ib];

    for (int j = 0; j < nintsRecv[i]; j++) {
      // Change first index from storing sndMap index to obblist index
      mb->isearch[2 * iboff + j] = i;
      mb->isearch[2 * iboff + j + 1] = rcvPack[k].intData[ioff + j];
    }

    for (int j = 0; j < nrealsRecv[i]; j++) {
      mb->xsearch[3 * iboff + j] = rcvPack[k].realData[3 * ioff + j];
    }

    ixOffset[k] += nintsRecv[i];
    ibOffset[ib] += 2 * nintsRecv[i];
  }
  pc->clearPackets(sndPack, rcvPack);

  free(sndPack);
  free(rcvPack);
  // printf("%d %d\n",myid,mb->nsearch);

  if (int_data) {
    for (int i=0; i<nobb; i++) {
      if (int_data[i]) free(int_data[i]);
    }
    free(int_data);
  }
  if (real_data) {
    for (int i=0; i<nobb; i++) {
      if (real_data[i]) free(real_data[i]);
    }
    free(real_data);
  }
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
  
