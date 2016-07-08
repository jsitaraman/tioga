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
  int i,j,k,l,m;
  int i3;
  int nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int meshtag,procid,pointid,remoteid;
  double donorRes;
  double *receptorResolution;
  int *donorRecords;
  int ninterp;
  int nrecords;
  int *sndMap;
  int *rcvMap;
  int *icount;
  int *dcount;
  //
  donorRecords=NULL;
  icount=NULL;
  dcount=NULL;
  receptorResolution=NULL;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) return;  
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nsend);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);
  //
  // get the data to send now
  //
  mb->getDonorPacket(sndPack,nsend);
  //
  // exchange the data (comm1)
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // packet received has all the donor data
  // use this to populate link lists per point
  //
  mb->initializeDonorList();
  //
  for(k=0;k<nrecv;k++)
    {
      m=0;
      for(i=0;i<rcvPack[k].nints/3;i++)
	{
	  meshtag=rcvPack[k].intData[m++];
	  pointid=rcvPack[k].intData[m++];
	  remoteid=rcvPack[k].intData[m++];
	  donorRes=rcvPack[k].realData[i];
	  mb->insertAndSort(pointid,k,meshtag,remoteid,donorRes);
	}
    }
  //
  // figure out the state of each point now (i.e. if its a hole or fringe or field)
  //
  mb->processDonors(holeMap,nmesh,&donorRecords,&receptorResolution,&nrecords);
  //
  // free Send and recv data
  //
  pc->clearPackets(sndPack,rcvPack);
  //
  // count number of records in each
  // packet
  //
  for(i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].nints++;
      sndPack[k].nreals++;
    }
  //
  // allocate the data containers
  //
  for(k=0;k<nsend;k++)
    {
     sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
     sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
     icount[k]=dcount[k]=0;
    }
  //
  for (i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
      sndPack[k].realData[dcount[k]++]=receptorResolution[i];
    }
  //
  // now communicate the data (comm2)
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // create the interpolation list now
  //
  ninterp=0;
  for(k=0;k<nrecv;k++)
    ninterp+=rcvPack[k].nints;
  //
  mb->initializeInterpList(ninterp);
  //
  m=0;
  for(k=0;k<nrecv;k++)
    for(i=0;i<rcvPack[k].nints;i++)
      mb->findInterpData(&m,rcvPack[k].intData[i],rcvPack[k].realData[i]);
  mb->set_ninterp(m);
  //
  //printf("process %d has (%d,%d) points to interpolate out %d donors\n",myid,ninterp,m,mb->donorCount);
  //
  pc->clearPackets(sndPack,rcvPack);
  free(donorRecords);
  donorRecords=NULL;
  //
  // cancel donors that have conflict 
  //
  mb->getCancellationData(&nrecords,&donorRecords);
  //printf("process %d has %d cancelled receptors\n",myid,nrecords);
  //
  
  for(i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].nints++;
   }
  for(k=0;k<nsend;k++)
    {
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      icount[k]=0;
    }  
  for(i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
    }
  //
  // communicate the cancellation data (comm 3)
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for(k=0;k<nrecv;k++)
    {
      for(i=0;i<rcvPack[k].nints;i++)
	mb->cancelDonor(rcvPack[k].intData[i]);
    }
  //
  if (donorRecords) free(donorRecords);  
  donorRecords=NULL;
  pc->clearPackets(sndPack,rcvPack);
  //
  mb->getInterpData(&nrecords,&donorRecords);
  //
  for(i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].nints++;
    }
  for(k=0;k<nsend;k++)
    {
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      icount[k]=0;
    }  
  for(i=0;i<nrecords;i++)
    {
      k=donorRecords[2*i];
      sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
    }
  //
  // communicate the final receptor data (comm 4)
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  mb->clearIblanks();
  for(k=0;k<nrecv;k++)
    for(i=0;i<rcvPack[k].nints;i++)
      mb->setIblanks(rcvPack[k].intData[i]);
  //
  // finished the communication, free all 
  // memory now
  //
  if (donorRecords) free(donorRecords);
  if (receptorResolution) free(receptorResolution);
  if (icount) free(icount);
  if (dcount) free(dcount);
  //
  // free Send and recv data
  //
  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
}
