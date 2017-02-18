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

void tioga::exchangeDonors(void)
{
  // get the processor map for sending and receiving
  int nsend, nrecv;
  int *sndMap, *rcvMap;
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  if (nsend == 0) return;

  int* icount = (int *)malloc(sizeof(int)*nsend);
  int* dcount = (int *)malloc(sizeof(int)*nsend);

  // create packets to send and receive
  // and initialize them to zero
  PACKET* sndPack = (PACKET *)malloc(sizeof(PACKET)*nsend);
  PACKET* rcvPack = (PACKET *)malloc(sizeof(PACKET)*nrecv);

  pc->initPackets(sndPack,rcvPack);

  // Get a list of donor-cell resolutions for all fringe points found within this grid
  mb->getDonorPacket(sndPack,nsend);

  // Exchange the data (comm1)
  pc->sendRecvPackets(sndPack,rcvPack);

  // Packet received has all the donor data for this grid's potential fringe points
  // Use this to populate link lists per point
  mb->initializeDonorList();

  for (int k = 0; k < nrecv; k++)
  {
    int m = 0;
    for (int i = 0; i < rcvPack[k].nints/3; i++)
    {
      int meshtag = rcvPack[k].intData[m++];
      int pointid = rcvPack[k].intData[m++];
      int remoteid = rcvPack[k].intData[m++];
      double donorRes = rcvPack[k].realData[i];
      mb->insertAndSort(pointid,k,meshtag,remoteid,donorRes);
    }
  }

  // figure out the state of each point now (i.e. if its a hole or fringe or field)
  double* receptorResolution = NULL;
  int* donorRecords = NULL;
  int nrecords;
  mb->processDonors(holeMap,nmesh,&donorRecords,&receptorResolution,&nrecords);
  /// ^ acutally processing fringe/receptor nodes based partly on their available donors

//  if (iartbnd)
//  {
//    free(donorRecords);
//    free(receptorResolution);
//    free(icount);
//    free(dcount);

//    // free Send and recv data
//    pc->clearPackets(sndPack,rcvPack);
//    free(sndPack);
//    free(rcvPack);

//    return;
//  }

  // free Send and recv data
  pc->clearPackets(sndPack,rcvPack);

  // count number of records in each packet
  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].nints++;
    sndPack[k].nreals++;
  }

  // allocate the data containers
  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double *)malloc(sizeof(double)*sndPack[k].nreals);
    icount[k] = dcount[k] = 0;
  }

  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].intData[icount[k]++] = donorRecords[2*i+1];
    sndPack[k].realData[dcount[k]++] = receptorResolution[i];
  }

  // now communicate the data (comm2)
  pc->sendRecvPackets(sndPack,rcvPack);

  // create the interpolation list now
  int ninterp=0;
  for (int k = 0; k < nrecv; k++)
    ninterp += rcvPack[k].nints;

  mb->initializeInterpList(ninterp);

  // Figure out which points from other grids we should be interpolating to
  int m = 0;
  for (int k = 0; k < nrecv; k++)
    for (int i = 0; i < rcvPack[k].nints; i++)
      mb->findInterpData(&m,rcvPack[k].intData[i],rcvPack[k].realData[i]);

  mb->set_ninterp(m);

  //printf("process %d has (%d,%d) points to interpolate out %d donors\n",myid,ninterp,m,mb->donorCount);

  pc->clearPackets(sndPack,rcvPack);
  free(donorRecords);
  donorRecords=NULL;

  // cancel donors that have conflict 
  mb->getCancellationData(nrecords, donorRecords);
  //printf("process %d has %d cancelled receptors\n",myid,nrecords);

  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].nints++;
  }

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    icount[k] = 0;
  }

  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].intData[icount[k]++] = donorRecords[2*i+1];
  }

  // communicate the cancellation data (comm 3)
  pc->sendRecvPackets(sndPack,rcvPack);

  for (int k = 0; k < nrecv; k++)
  {
    for (int i = 0; i < rcvPack[k].nints; i++)
      mb->cancelDonor(rcvPack[k].intData[i]);
  }

  free(donorRecords);
  donorRecords = NULL;
  pc->clearPackets(sndPack,rcvPack);

  mb->getInterpData(nrecords, donorRecords);

  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].nints++;
  }

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    icount[k] = 0;
  }

  for (int i = 0; i < nrecords; i++)
  {
    int k = donorRecords[2*i];
    sndPack[k].intData[icount[k]++] = donorRecords[2*i+1];
  }

  // communicate the final receptor data (comm 4)
  pc->sendRecvPackets(sndPack,rcvPack);
  mb->clearIblanks();
  for (int k = 0; k < nrecv; k++)
    for (int i = 0; i < rcvPack[k].nints; i++)
      mb->setIblanks(rcvPack[k].intData[i]);

  // finished the communication, free all memory now
  free(donorRecords);
  free(receptorResolution);
  free(icount);
  free(dcount);

  // free Send and recv data
  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
}
