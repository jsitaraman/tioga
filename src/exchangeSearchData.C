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

  for (int ii=0; ii < nobb; ii++) {
    int ib = obblist[ii].iblk_local;
    auto& mb = mblocks[ib];
    mb->getQueryPoints(
      &obblist[ii], &nintsSend[ii], &int_data[ii], &nrealsSend[ii],
      &real_data[ii]);
  }

  // Populate send packets and exchange data with other processors
  for (int k=0; k<nsend; k++) {
    sndPack[k].nints = 3 * ibsPerProc[k];
    sndPack[k].nreals = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].nints += nintsSend[ii];
      sndPack[k].nreals += nrealsSend[ii];
    }
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
    sndPack[k].realData = (double*)malloc(sizeof(double) * sndPack[k].nreals);

    int n = 0;
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].intData[n++] = obblist[ii].send_tag;
      sndPack[k].intData[n++] = nintsSend[ii];
      sndPack[k].intData[n++] = nrealsSend[ii];

      for (int j=0; j < nintsSend[ii]; j++)
        sndPack[k].intData[n++] = int_data[ii][j];

      for (int j=0; j < nrealsSend[ii]; j++)
        sndPack[k].realData[m++] = real_data[ii][j];
    }
  }
  pc->sendRecvPackets(sndPack, rcvPack);

  // Reset MeshBlock data structures
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb = mblocks[ib];
    mb->nsearch = 0;

    if (mb->xsearch)
      free(mb->xsearch);
    if (mb->isearch)
      free(mb->isearch);
    if (mb->tagsearch)
      free(mb->tagsearch);
    if (mb->donorId)
      free(mb->donorId);
    if (mb->res_search)
      free(mb->res_search);
  }

  // Loop through recv packets and estimate search array sizes in MeshBlock data
  // structures
  std::vector<int> nintsRecv(nobb);
  std::vector<int> nrealsRecv(nobb);
  for (int k=0; k<nrecv; k++) {
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      nintsRecv[ii] = rcvPack[k].intData[m++];
      nrealsRecv[ii] = rcvPack[k].intData[m++];

      mb->nsearch += nintsRecv[ii];
      // Skip the node indices
      m += nintsRecv[ii];
    }
  }

  // Resize MeshBlock array sizes
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb = mblocks[ib];
    if (mb->nsearch < 1) continue;
    mb->xsearch = (double*)malloc(sizeof(double) * 3 * mb->nsearch);
    mb->res_search = (double*)malloc(sizeof(double) * mb->nsearch);
    mb->isearch = (int*)malloc(3 * sizeof(int) * mb->nsearch);
    mb->tagsearch = (int*)malloc(sizeof(int) * mb->nsearch);
    mb->donorId = (int*)malloc(sizeof(int) * mb->nsearch);
  }
  //
  // Update search arrays in mesh blocks from recv packets
  //
  std::vector<int> icOffset(nblocks,0); // Index of isearch arrays where next fill happens
  std::vector<int> dcOffset(nblocks, 0); // Index of xsearch arrays where next fill happens
  std::vector<int> rcOffset(nblocks, 0); // Index of res_search arrays where next fill happens
  for (int k=0; k < nrecv; k++) {
    int l = 0;
    int m = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      int ioff = icOffset[ib];
      int doff = dcOffset[ib];
      int roff = rcOffset[ib];

      m += 2; // Skip nints and nreals information
      for (int j=0; j < nintsRecv[ii]; j++) {
        mb->isearch[ioff++] = k;
        mb->isearch[ioff++] = rcvPack[k].intData[m++];
	mb->isearch[ioff++] = obblist[ii].iblk_remote;
        mb->tagsearch[ioff/3-1]=obblist[ii].tag_remote;
      }

      for (int j=0; j < nrealsRecv[ii]/4; j++) {
        for(int mm=0;mm<3;mm++)
          mb->xsearch[doff++] = rcvPack[k].realData[l++];
        mb->res_search[roff++]= rcvPack[k].realData[l++];
      }

      icOffset[ib] = ioff;
      dcOffset[ib] = doff;
      rcOffset[ib] = roff;
    }
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
  
