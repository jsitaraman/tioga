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
#include "tioga.h"

void tioga::exchangeDonors(void)
{
  int nsend,nrecv;
  int *sndMap;
  int *rcvMap;
  PACKET *sndPack,*rcvPack;

  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) return;  
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);

  int nobb = obblist.size();
  std::vector<int> nintsSend(nobb,0), nrealsSend(nobb,0);

  // First pass through all mesh blocks and determine data sizes for each
  // intersection pair
  for (auto& mb: mblocks) {
    mb->getMBDonorPktSizes(nintsSend, nrealsSend);
  }

  // Setup sndPack 
  std::vector<int> ixOffset(nobb), rxOffset(nobb);
  for (int k=0; k<nsend; k++) {
    // (SEND_TAG, nints, nreals)
    sndPack[k].nints = 3 * ibsPerProc[k];

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

      ixOffset[ii] = n;
      n += nintsSend[ii];
      rxOffset[ii] = m;
      m += nrealsSend[ii];
    }
  }

  // Populate send packets with data from each mesh block in this partition
  for (auto& mb: mblocks) {
    mb->getMBDonorPackets(obblist, ixOffset, rxOffset, sndPack);
  }
  pc->sendRecvPackets(sndPack,rcvPack);

  // Initialize linked lists and populate donor data from rcvPack
  for (auto& mb: mblocks) {
    mb->initializeDonorList();
  }

  for (int k=0; k<nrecv; k++) {
    int m = 0;
    int l = 0;

    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      m++; // int nints = rcvPack[k].intData[m++];
      int nreals = rcvPack[k].intData[m++];

      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      for (int j=0; j < nreals; j++) {
        int meshtag = rcvPack[k].intData[m++];
        int pointid = rcvPack[k].intData[m++];
        int remoteid = rcvPack[k].intData[m++];
        double donorRes = rcvPack[k].realData[l++];

        mb->insertAndSort(pointid, ii, meshtag, remoteid, donorRes);
      }
    }
  }

  // Figure out the state of each point (i.e., if it is a hole, fringe, or a
  // field point)
  std::vector<int> nrecords(nblocks,0);
  int** donorRecords = (int**)malloc(sizeof(int*)*nblocks);
  double** receptorResolution = (double**)malloc(sizeof(double*)*nblocks);

  for (int i=0; i<nblocks; i++) {
    auto& mb = mblocks[i];
    mb->processDonors(holeMap, nmesh, &(donorRecords[i]),
                      &(receptorResolution[i]),&(nrecords[i]));
  }

  // Reset all send/recv data structures
  pc->clearPackets(sndPack, rcvPack);
  for (int i=0; i<nobb; i++) {
    nintsSend[i] = 0;
    nrealsSend[i] = 0;
    ixOffset[i] = 0;
    rxOffset[i] = 0;
  }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      nintsSend[ii]++;
      nrealsSend[ii]++;
    }
  }

  for(int k=0; k<nsend; k++) {
    sndPack[k].nints = 3 * ibsPerProc[k];

    for (int i=0; i< ibsPerProc[k]; i++) {
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

      ixOffset[ii] = n;
      n += nintsSend[ii];
      rxOffset[ii] = m;
      m += nrealsSend[ii];
    }
  }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i < nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      int k = obblist[ii].comm_idx;

      sndPack[k].intData[ixOffset[ii]++] = donorRecords[n][2*i+1];
      sndPack[k].realData[rxOffset[ii]++] = receptorResolution[n][i];
    }
  }

  pc->sendRecvPackets(sndPack,rcvPack);

  std::vector<int> ninterp(nblocks,0);

  for (int k=0; k<nrecv; k++) {
    int m = 0;
    for (int i=0; i < ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;

      int nints =  rcvPack[k].intData[m++];
      ninterp[ib] += nints;
      m += (1 + nints); // Skip nreals and actual nints info
    }
  }

  for (int i=0; i<nblocks; i++) {
    mblocks[i]->initializeInterpList(ninterp[i]);
  }
  std::fill(ninterp.begin(), ninterp.end(), 0);

  for (int k=0; k<nrecv; k++) {
    int l = 0;
    int m = 0;

    for (int i=0; i<ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      m++; // int nints = rcvPack[k].intData[m++];
      int nreals = rcvPack[k].intData[m++];

      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      for (int j=0; j<nreals; j++) {
        mb->findInterpData(&(ninterp[ib]),rcvPack[k].intData[m++],
                           rcvPack[k].realData[l++]);
      }
    }
  }

  for (int i=0; i<nblocks; i++) {
    mblocks[i]->set_ninterp(ninterp[i]);
  }

  pc->clearPackets(sndPack, rcvPack);

  // Find cancellation data
  for (int i=0; i<nblocks; i++) {
    if (donorRecords[i]) {
      free(donorRecords[i]);
      donorRecords[i] = NULL;
    }
    nrecords[i] = 0;

    mblocks[i]->getCancellationData(&(nrecords[i]),
                                    &(donorRecords[i]));
  }
  std::fill(nintsSend.begin(), nintsSend.end(), 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);

  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      nintsSend[ii]++;
    }
  }

  for (int k=0; k<nsend; k++) {
    sndPack[k].nints = 2 * ibsPerProc[k];
    sndPack[k].nreals = 0;
    sndPack[k].realData = NULL;

    for (int i=0; i<ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];
      sndPack[k].nints += nintsSend[ii];
    }

    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);

    int n = 0;
    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].intData[n++] = obblist[ii].send_tag;
      sndPack[k].intData[n++] = nintsSend[ii];

      ixOffset[ii] = n;
      n += nintsSend[ii];
    }
  }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i < nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      int k = obblist[ii].comm_idx;

      sndPack[k].intData[ixOffset[ii]++] = donorRecords[n][2*i+1];
    }
  }

  pc->sendRecvPackets(sndPack,rcvPack);

  for (int k=0; k<nrecv; k++) {
    int m = 0;

    for (int i=0; i<ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int nints = rcvPack[k].intData[m++];

      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      for (int j=0; j<nints; j++) {
        mb->cancelDonor(rcvPack[k].intData[m++]);
      }
    }
  }

  pc->clearPackets(sndPack, rcvPack);

  // Find cancellation data
  for (int i=0; i<nblocks; i++) {
    if (donorRecords[i]) {
      free(donorRecords[i]);
      donorRecords[i] = NULL;
    }
    nrecords[i] = 0;

    mblocks[i]->getInterpData(&(nrecords[i]),
                              &(donorRecords[i]));
  }

  std::fill(nintsSend.begin(), nintsSend.end(), 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);

  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      nintsSend[ii]++;
    }
  }

  for (int k=0; k<nsend; k++) {
    sndPack[k].nints = 2 * ibsPerProc[k];
    sndPack[k].nreals = 0;
    sndPack[k].realData = NULL;

    for (int i=0; i<ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];
      sndPack[k].nints += nintsSend[ii];
    }

    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);

    int n = 0;
    for (int i=0; i < ibsPerProc[k]; i++) {
      int ii = ibProcMap[k][i];

      sndPack[k].intData[n++] = obblist[ii].send_tag;
      sndPack[k].intData[n++] = nintsSend[ii];

      ixOffset[ii] = n;
      n += nintsSend[ii];
    }
  }

  for (int n=0; n<nblocks; n++) {
    for (int i=0; i < nrecords[n]; i++) {
      int ii = donorRecords[n][2*i];
      int k = obblist[ii].comm_idx;

      sndPack[k].intData[ixOffset[ii]++] = donorRecords[n][2*i+1];
    }
  }

  pc->sendRecvPackets(sndPack,rcvPack);

  for (auto& mb: mblocks)
    mb->clearIblanks();

  for (int k=0; k<nrecv; k++) {
    int m = 0;

    for (int i=0; i<ibsPerProc[k]; i++) {
      int key = rcvPack[k].intData[m++];
      int nints = rcvPack[k].intData[m++];

      int ii = intBoxMap[key];
      int ib = obblist[ii].iblk_local;
      auto& mb = mblocks[ib];

      for (int j=0; j<nints; j++) {
        mb->setIblanks(rcvPack[k].intData[m++]);
      }
    }
  }

  // for (int i=0; i<nobb; i++) {
  //   std::cout << myid << " " << i << " " << ixOffset[i] <<  " "
  //             << nintsSend[i] << std::endl;
  //   std::cout << myid << " " << i << " " << ixOffset[i] << " " << rxOffset[i] << " "
  //             << nintsSend[i] << " " << nrealsSend[i] << std::endl;
  // }

  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);

  if (donorRecords) {
    for (int i=0; i<nblocks; i++) {
      if (donorRecords[i]) free(donorRecords[i]);
    }
    free(donorRecords);
  }
  if (receptorResolution) {
    for (int i=0; i<nblocks; i++) {
      if (receptorResolution[i]) free(receptorResolution[i]);
    }
    free(receptorResolution);
  }
}

void tioga::outputStatistics(void)
{
  int mstats[2], mstats_sum[2], mstats_global[2];
  mstats_sum[0] = 0;
  mstats_sum[1] = 0;
  for (auto& mb: mblocks) {
    mb->getStats(mstats);
    mstats_sum[0] += mstats[0];
    mstats_sum[1] += mstats[1];
  }
  MPI_Reduce(mstats,mstats_global,2,MPI_INT,MPI_SUM,0,scomm);
  if (myid==0) {
    printf("#tioga -----------------------------------------\n");
    printf("#tioga : total receptors:\t%d\n",mstats_global[1]);
    printf("#tioga : total holes    :\t%d\n",mstats_global[0]);
    printf("#tioga -----------------------------------------\n");
  }
}
