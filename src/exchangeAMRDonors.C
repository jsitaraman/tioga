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

void tioga::exchangeAMRDonors(void)
{
  int i,j,k,l,m,n,i3;
  int nsend_sav,nrecv_sav,nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int *sndMap,*rcvMap,*sndMapAll,*rcvMapAll;
  int *imap,*icount;
  int *obdonors,*obreceptors;
  int *cancelledData;
  int gid,procid,localid,meshtag,remoteid,id;
  int interpCount,donorCount,index,ncancel,ctype;
  double cellRes;
  double xtmp[3];
  //
  // setup communicator for all to all now
  // since the receiver side is unknown
  // FIXME:
  // add sophisticated code later to fix the all_to_all
  // using MPI-2 standard
  //
  pc_cart->getMap(&nsend_sav,&nrecv_sav,&sndMap,&rcvMap);
  sndMapAll=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  rcvMapAll=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  nsend=nrecv=pc_cart->numprocs;
  nsend=pc_cart->numprocs;
  imap=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  icount=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  for(i=0;i<pc_cart->numprocs;i++) {sndMapAll[i]=rcvMapAll[i]=imap[i]=i; icount[i]=0;}
  pc_cart->setMap(nsend,nrecv,sndMapAll,rcvMapAll);
  //  for(i=0;i<nsend;i++) imap[sndMap[i]]=i;
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  obdonors=(int *)malloc(sizeof(int)*nsend);
  obreceptors=(int *)malloc(sizeof(int)*nsend);
  //
  pc_cart->initPackets(sndPack,rcvPack);
  //
  for(i=0;i<nsend;i++) obdonors[i]=obreceptors[i]=0;
  //
  if (nblocks > 0) 
    {
      //
      // count data to send
      //
      for(i=0;i<mb->ntotalPointsCart;i++) 
	{ 
	  if (mb->donorIdCart[i]!=-1)
	    {
	      gid=mb->donorIdCart[i];
	      obdonors[imap[cg->proc_id[gid]]]++;
	    }
	}
      for(i=0;i<mb->nsearch;i++) 
	{
	  if (mb->donorId[i]!=-1) 
	    {
	  obreceptors[imap[mb->isearch[3*i]]]++;
	    }
	}
    }
  //
  // allocate data packets
  //
  for(i=0;i<nsend;i++)
    {
      sndPack[i].nints=obdonors[i]*2+obreceptors[i]*4+2;
      sndPack[i].nreals=obdonors[i]*3+obreceptors[i];
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
      sndPack[k].intData[0]=obdonors[i];
      sndPack[k].intData[1]=obreceptors[i];
    }
  //
  // pack the data
  //
  m=2;
  n=0;
  if (nblocks > 0) 
    {
      for(i=0;i<mb->ntotalPointsCart;i++) 
	{ 
	  if (mb->donorIdCart[i]!=-1)
	    {
	      gid=mb->donorIdCart[i];
	      procid=imap[cg->proc_id[gid]];
	      localid=cg->local_id[gid];	  
	      sndPack[procid].intData[m++]=localid;
	      sndPack[procid].intData[m++]=i;
	      sndPack[procid].realData[n++]=mb->rxyzCart[3*i];
	      sndPack[procid].realData[n++]=mb->rxyzCart[3*i+1];
	      sndPack[procid].realData[n++]=mb->rxyzCart[3*i+2];	  
	    }
	}
      for(i=0;i<mb->nsearch;i++) 
	{
	  if (mb->donorId[i]!=-1) 
	    {
	      procid=imap[mb->isearch[3*i]];
	      sndPack[procid].intData[m++]=mb->isearch[3*i+1];
	      sndPack[procid].intData[m++]=mb->isearch[3*i+2];
	      sndPack[procid].intData[m++]=mb->meshtag;
	      sndPack[procid].intData[m++]=i;
	      sndPack[procid].realData[n++]=mb->cellRes[mb->donorId[i]];
	    }
	}
    }
  //
  // exchange data
  //
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the data now
  //
  for(i=0;i<ncart;i++) cb[i].initializeLists();
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nreals > 0) 
	{
	  m=0;
	  n=0;
	  interpCount=rcvPack[i].intData[0];
	  donorCount=rcvPack[i].intData[1];
	  icount[i]=interpCount+donorCount;
	  for(j=0;j<interpCount;j++)
	    {
	      localid=rcvPack[i].intData[m++];
	      remoteid=rcvPack[i].intData[m++];
	      xtmp[0]=sndPack[i].realData[n++];
	      xtmp[1]=sndPack[i].realData[n++];
	      xtmp[2]=sndPack[i].realData[n++];	      
	      cb[localid].insertInInterpList(i,remoteid,xtmp);
	    }
	  for(j=0;j<donorCount;j++)
	    {
	      localid=rcvPack[i].intData[m++];
	      index=rcvPack[i].intData[m++];
	      meshtag=rcvPack[i].intData[m++];
	      remoteid=rcvPack[i].intData[m++];
	      cellRes=rcvPack[i].realData[n++];
	      cb[localid].insertInDonorList(i,index,meshtag,remoteid,cellRes);
	    }
	}
    }
  for(i=0;i<ncart;i++) cb[i].processDonors(holeMap,nmesh);
  pc_cart->clearPackets2(sndPack,rcvPack);  
  for(i=0;i<nsend;i++)
    {
      if (icount[i] > 0) 
	{
	  sndPack[i].nints=2*icount[i];
	  sndPack[i].nreals=0;
	  sndPack[i].intData=(int *) malloc(sizeof(int)*sndPack[i].nints);
	  icount[i]=0;
	}
    }
  for(i=0;i<ncart;i++) 
    {
      if (cancelledData) free(cancelledData);
      cancelledData=NULL;
      ncancel=icount[i];
      if (ncancel > 0) {
	cancelledData=(int *)malloc(sizeof(int)*3*icount[i]);
	cb[i].getCancellationData(cancelledData,&ncancel);
	for(j=0;j<ncancel;j++)
	  {
	    procid=cancelledData[3*j];
	    ctype=cancelledData[3*j+1];
	    remoteid=cancelledData[3*j+2];
	    sndPack[procid].intData[icount[procid]++]=ctype;
	    sndPack[procid].intData[icount[procid]++]=remoteid;
	  }
      }
    }
  for(i=0;i<nsend;i++) sndPack[i].nints=2*icount[i];
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nints > 0) 
	{
	  m=0;
	  for(i=0;i<rcvPack[i].nints/2;i++)
	    {
	      ctype=rcvPack[i].intData[2*i];
	      id=rcvPack[i].intData[2*i+1];
	      if (ctype==0) 
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
  pc_cart->clearPackets2(sndPack,rcvPack);
  free(sndMap);
  free(rcvMap);
  free(imap);
  free(icount);
  free(sndPack);
  free(rcvPack);
}






































































































































































































































































