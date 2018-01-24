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
#include <assert.h>
void tioga::exchangeAMRDonors(int itype)
{
  //FILE *fp;
  //char fname[80];
  //char qstr[2];
  //char intstring[7];

  //sprintf(intstring,"%d",100000+myid);
  //sprintf(fname,"ysearch_%s.dat",&(intstring[1]));  
  //fp=fopen(fname,"w");
  //
  // setup communicator for all to all now
  // since the receiver side is unknown
  // FIXME:
  // add sophisticated code later to fix the all_to_all
  // using MPI-2 standard

  int nsend_sav,nrecv_sav;
  int *sndMap,*rcvMap;

  pc_cart->getMap(&nsend_sav,&nrecv_sav,&sndMap,&rcvMap);

  int* sndMapAll = (int *)malloc(sizeof(int)*pc_cart->numprocs);
  int* rcvMapAll = (int *)malloc(sizeof(int)*pc_cart->numprocs);
  int nsend = pc_cart->numprocs;
  int nrecv = pc_cart->numprocs;
  int* imap = (int *)malloc(sizeof(int)*pc_cart->numprocs);
  int* icount = (int *)malloc(sizeof(int)*pc_cart->numprocs);

  for (int i = 0; i < pc_cart->numprocs; i++)
  {
    sndMapAll[i] = rcvMapAll[i] = imap[i] = i;
    icount[i] = 0;
  }

  pc_cart->setMap(nsend,nrecv,sndMapAll,rcvMapAll);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  PACKET* sndPack = (PACKET *)malloc(sizeof(PACKET)*nsend);
  PACKET* rcvPack = (PACKET *)malloc(sizeof(PACKET)*nrecv);
  int* obdonors = (int *)malloc(sizeof(int)*nsend);
  int* obreceptors = (int *)malloc(sizeof(int)*nsend);

  pc_cart->initPackets(sndPack,rcvPack);

  for (int i = 0; i < nsend; i++) obdonors[i] = obreceptors[i] = 0;

  if (nblocks > 0)
  {
    //
    // count data to send
    //
    for (int i = 0; i < mb->ntotalPointsCart; i++)
    {
      if (mb->donorIdCart[i] != -1)
      {
        int procid = cg->proc_id[mb->donorIdCart[i]];
        assert((procid < nsend && procid >= 0));
        obdonors[imap[procid]]++;
      }
    }

    for (int i = 0; i < mb->nsearch; i++)
    {
      if (mb->donorId[i] != -1)
      {
        assert((mb->isearch[3*i] < nsend && mb->isearch[3*i]>=0));
        obreceptors[imap[mb->isearch[3*i]]]++;
      }
    }
  }

  //
  // allocate data packets
  //
  for (int i = 0; i < nsend; i++)
  {
    sndPack[i].nints = obdonors[i]*2+obreceptors[i]*4+2;
    sndPack[i].nreals = obdonors[i]*3+obreceptors[i];
    sndPack[i].intData = (int *)malloc(sizeof(int)*sndPack[i].nints);
    sndPack[i].realData = (double *)malloc(sizeof(double)*sndPack[i].nreals);
    sndPack[i].intData[0] = obdonors[i];
    sndPack[i].intData[1] = obreceptors[i];
  }

  //
  // pack the data
  //
  int* intcount = (int *)malloc(sizeof(int)*nsend);
  int* realcount = (int *)malloc(sizeof(int)*nsend);

  for (int i = 0; i < nsend; i++) { intcount[i]=2; realcount[i]=0; }

  if (nblocks > 0)
  {
    for (int i = 0; i < mb->ntotalPointsCart; i++)
    {
      if (mb->donorIdCart[i] != -1)
      {
        int gid = mb->donorIdCart[i];
        int procid = imap[cg->proc_id[gid]];
        int localid = cg->local_id[gid];
        sndPack[procid].intData[intcount[procid]++] = localid;
        sndPack[procid].intData[intcount[procid]++] = i;
        sndPack[procid].realData[realcount[procid]++] = mb->rxyzCart[3*i];
        sndPack[procid].realData[realcount[procid]++] = mb->rxyzCart[3*i+1];
        sndPack[procid].realData[realcount[procid]++] = mb->rxyzCart[3*i+2];
      }
    }

    for (int i = 0; i < mb->nsearch; i++)
    {
      if (mb->donorId[i] != -1)
      {
        int procid = mb->isearch[3*i];
        sndPack[procid].intData[intcount[procid]++] = mb->isearch[3*i+1];
        sndPack[procid].intData[intcount[procid]++] = mb->isearch[3*i+2];
        sndPack[procid].intData[intcount[procid]++] = mb->meshtag;
        sndPack[procid].intData[intcount[procid]++] = i;
        sndPack[procid].realData[realcount[procid]++] = mb->cellRes[mb->donorId[i]];
        //      fprintf(fp,"%lf %lf %lf\n",mb->xsearch[3*i],mb->xsearch[3*i+1],mb->xsearch[3*i+2]);
      }
    }
  }
  //if (myid==0) {
  //for(i=0;i<nsend;i++)
  //  printf("intcount/intcount=%d %d\n",sndPack[i].nints,intcount[i]);
  //for(i=0;i<nsend;i++)
  //  printf("intcount/intcount=%d %d\n",sndPack[i].nreals,realcount[i]);
  //}
  free(realcount);
  //
  // exchange data
  //
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the data now
  //
  int* bcount = (int *)malloc(sizeof(int)*ncart);
  for (int i = 0; i < ncart; i++)
  {
    cb[i].clearLists(); 
    cb[i].initializeLists();
    bcount[i] = 0;
  }

  double xtmp[3] = {};
  for (int i = 0; i < nrecv; i++)
  {
    if (rcvPack[i].nreals > 0)
    {
      int m = 2;
      int n = 0;
      int interpCount = rcvPack[i].intData[0];
      int donorCount = rcvPack[i].intData[1];
      icount[i] = interpCount+donorCount;

      // Get interpolation weights from AMR grid for this grid's receptors
      for (int j = 0; j < interpCount; j++)
      {
        int localid=rcvPack[i].intData[m++];
        int remoteid=rcvPack[i].intData[m++];
        xtmp[0]=rcvPack[i].realData[n++];
        xtmp[1]=rcvPack[i].realData[n++];
        xtmp[2]=rcvPack[i].realData[n++];
        cb[localid].insertInInterpList(i,remoteid,xtmp,itype);
        bcount[localid]++;
      }

      for (int j = 0; j < donorCount; j++)
      {
        int localid = rcvPack[i].intData[m++];
        int index = rcvPack[i].intData[m++];
        int meshtag = rcvPack[i].intData[m++];
        int remoteid = rcvPack[i].intData[m++];
        int cellRes = rcvPack[i].realData[n++];
        cb[localid].insertInDonorList(i,index,meshtag,remoteid,cellRes,itype);
        bcount[localid]++;
      }
    }
  }

  for (int i = 0; i < ncart; i++) cb[i].processDonors(holeMap,nmesh,itype);

  pc_cart->clearPackets2(sndPack,rcvPack);  

  for (int i = 0; i < nsend; i++)
  {
    if (icount[i] > 0)
    {
      sndPack[i].nints = 2*icount[i];
      sndPack[i].nreals = 0;
      sndPack[i].intData = (int *) malloc(sizeof(int)*sndPack[i].nints);
    }
    intcount[i] = 0;
  }

  int* cancelledData = NULL;

  for (int i = 0; i < ncart; i++)
  {
    //tracei(i);
    free(cancelledData);
    cancelledData = NULL;
    int ncancel = bcount[i];
    if (ncancel > 0)
    {
      cancelledData = (int *)malloc(sizeof(int)*3*ncancel);
      cb[i].getCancellationData(cancelledData,&ncancel,1);

      for (int j = 0; j < ncancel; j++)
      {
        int procid = cancelledData[3*j];
        int ctype = cancelledData[3*j+1];
        int remoteid = cancelledData[3*j+2];
        sndPack[procid].intData[intcount[procid]++]=ctype;
        sndPack[procid].intData[intcount[procid]++]=remoteid;
        //tracei(intcount[procid]);
        //tracei(sndPack[procid].nints);
      }
    }
  }

  for (int i = 0; i < nsend; i++) sndPack[i].nints = intcount[i];

  pc_cart->sendRecvPackets(sndPack,rcvPack);

  for (int i = 0; i < nrecv; i++)
  {
    if (rcvPack[i].nints > 0)
    {
      int m = 0;
      for (int j = 0; j < rcvPack[i].nints/2; j++)
      {
        int ctype = rcvPack[i].intData[m++];
        int id = rcvPack[i].intData[m++];
        if (ctype == 0)
        {
          mb->donorIdCart[id]=-1;
        }
        else
        {
          mb->donorId[id]=-1;
        }
      }
    }
  }

  if (itype == 0) mb->setCartIblanks();

  pc_cart->clearPackets2(sndPack,rcvPack);

  mb->findInterpListCart();

  if (cancelledData) free(cancelledData);
  //fclose(fp);
  free(bcount);
  free(intcount);
  free(sndMapAll);
  free(rcvMapAll);
  free(imap);
  free(icount);
  free(sndPack);
  free(rcvPack);
}
void tioga::checkComm(void)
{
  int i;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;

  nsend=nrecv=pc_cart->numprocs;
  sndMap=(int *)malloc(sizeof(int)*nsend);
  rcvMap=(int *)malloc(sizeof(int)*nrecv);
  for(i=0;i<pc_cart->numprocs;i++) {sndMap[i]=rcvMap[i]=i;}
  pc_cart->setMap(nsend,nrecv,sndMap,rcvMap);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc_cart->initPackets(sndPack,rcvPack);
  //
  // allocate data packets
  //
  for(i=0;i<nsend;i++)
    {
      sndPack[i].nints=2;
      sndPack[i].nreals=2;//obdonors[i]*3+obreceptors[i]+1;
      sndPack[i].intData=(int *)malloc(sizeof(int)*sndPack[i].nints);
      sndPack[i].realData=(double *)malloc(sizeof(double)*sndPack[i].nreals);
      sndPack[i].intData[0]=0;
      sndPack[i].intData[1]=1;
    }
  //
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  pc_cart->clearPackets2(sndPack,rcvPack);
  free(sndMap);
  free(rcvMap);
  free(sndPack);
  free(rcvPack);
  printf("checkComm complete in %d\n",myid);
}
