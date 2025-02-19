// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include <algorithm>
#include <vector>
#include <cstring>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
void tioga::exchangeSearchData(int at_points)
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
    if (at_points==0) 
    {
      mb->getQueryPoints2(
      &obblist[ii], &nintsSend[ii], &int_data[ii], &nrealsSend[ii],
      &real_data[ii]);
    } 
    else
    {
     mb->getExtraQueryPoints(&obblist[ii],
                            &nintsSend[ii],&int_data[ii],
                            &nrealsSend[ii],&real_data[ii]);
    }
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

    if (mb->xsearch) {
      TIOGA_FREE(mb->xsearch);
      mb->xsearch=NULL;
    }
    if (mb->isearch) {
      TIOGA_FREE(mb->isearch);
      mb->isearch=NULL;
    }
    if (mb->tagsearch) {
      TIOGA_FREE(mb->tagsearch);
      mb->tagsearch=NULL;
    }
    if (mb->donorId) {
      TIOGA_FREE(mb->donorId);
      mb->donorId=NULL;
    }
    if (mb->res_search) {
      TIOGA_FREE(mb->res_search);
      mb->res_search=NULL;
    }
   if (at_points==1) {
     if (mb->rst) {
       TIOGA_FREE(mb->rst);
       mb->rst=NULL;
     }
    }
#ifdef TIOGA_HAS_NODEGID
   mb->gid_search.clear();
#endif
  }

#ifdef TIOGA_HAS_NODEGID
  int nintsPerNode = 3;
#else
  int nintsPerNode = 1;
#endif

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

      // Skip the nodal data (indices and global IDs)
      m += nintsRecv[ii];
      // Adjust for node GIDs if available
      nintsRecv[ii] /= nintsPerNode;
      mb->nsearch += nintsRecv[ii];
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
#ifdef TIOGA_HAS_NODEGID
    mb->gid_search.resize(mb->nsearch);
#endif
    if (at_points==1) mb->rst = (double*)malloc(sizeof(double) * 3 * mb->nsearch);
  }

  //
  // Update search arrays in mesh blocks from recv packets
  //
  std::vector<int> icOffset(nblocks,0); // Index of isearch arrays where next fill happens
  std::vector<int> dcOffset(nblocks, 0); // Index of xsearch arrays where next fill happens
  std::vector<int> rcOffset(nblocks, 0); // Index of res_search arrays where next fill happens
  std::vector<int> igOffset(nblocks, 0); // Index of gid_search where next fill happens
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
      int goff = igOffset[ib];

      m += 2; // Skip nints and nreals information
      for (int j=0; j < nintsRecv[ii]; j++) {
        mb->isearch[ioff++] = k;
        mb->isearch[ioff++] = rcvPack[k].intData[m++];
        mb->isearch[ioff++] = obblist[ii].iblk_remote;
        mb->tagsearch[ioff/3-1]=obblist[ii].tag_remote;

#ifdef TIOGA_HAS_NODEGID
        std::memcpy(&(mb->gid_search[goff]), &(rcvPack[k].intData[m]), sizeof(uint64_t));
        m += 2;
        goff++;
#endif
      }

      for (int j = 0; j < nrealsRecv[ii] / 4; j++) {
        for (int mm = 0; mm < 3; mm++)
          mb->xsearch[doff++] = rcvPack[k].realData[l++];
        mb->res_search[roff++] = rcvPack[k].realData[l++];
      }

      icOffset[ib] = ioff;
      dcOffset[ib] = doff;
      rcOffset[ib] = roff;
      igOffset[ib] = goff;
    }
  }

  pc->clearPackets(sndPack, rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  // printf("%d %d\n",myid,mb->nsearch);

  if (int_data) {
    for (int i=0; i<nobb; i++) {
      if (int_data[i]) TIOGA_FREE(int_data[i]);
    }
    TIOGA_FREE(int_data);
  }
  if (real_data) {
    for (int i=0; i<nobb; i++) {
      if (real_data[i]) TIOGA_FREE(real_data[i]);
    }
    TIOGA_FREE(real_data);
  }
}
