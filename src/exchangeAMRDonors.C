// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"
#include "tioga.h"
#include <assert.h>
using namespace TIOGA;
void tioga::exchangeAMRDonors(void)
{
  int i,j,k,l,m,n,i3;
  int nsend_sav,nrecv_sav,nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int *sndMap,*rcvMap,*sndMapAll,*rcvMapAll;
  int *imap,*icount,*intcount,*realcount;
  int *obdonors,*obreceptors;
  int *cancelledData;
  int *bcount;
  int gid,procid,localid,meshtag,remoteid,remoteblockid,id;
  int interpCount,donorCount,index,ncancel,ctype;
  double cellRes;
  double xtmp[3];
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
  //
  pc_cart->getMap(&nsend_sav,&nrecv_sav,&sndMap,&rcvMap);
  sndMapAll=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  rcvMapAll=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  nsend=nrecv=pc_cart->numprocs;
  imap=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  icount=(int *)malloc(sizeof(int)*pc_cart->numprocs);
  for(i=0;i<pc_cart->numprocs;i++) {sndMapAll[i]=rcvMapAll[i]=imap[i]=i; icount[i]=0;}
  pc_cart->setMap(nsend,nrecv,sndMapAll,rcvMapAll);
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
    for(int ib=0;ib<nblocks;ib++)
      {
      auto& mb = mblocks[ib];
      for(i=0;i<mb->ntotalPointsCart;i++) 
	{ 
	  if (mb->donorIdCart[i]!=-1)
	    {
	      gid=mb->donorIdCart[i];
              assert((cg->proc_id[gid] < nsend && cg->proc_id[gid]>=0));
	      obdonors[imap[cg->proc_id[gid]]]++;
	    }
	}
       for(i=0;i<mb->nsearch;i++) 
	{
	  if (mb->donorId[i]!=-1) 
	    {
              assert((mb->isearch[3*i] < nsend && mb->isearch[3*i]>=0));
	      obreceptors[imap[mb->isearch[3*i]]]++;
	    }
	}
      } 
    }
  //
  // allocate data packets
  //
  for(i=0;i<nsend;i++)
    {
      sndPack[i].nints=obdonors[i]*3+obreceptors[i]*5+2;
      sndPack[i].nreals=obdonors[i]*3+obreceptors[i];
      sndPack[i].intData=(int *)malloc(sizeof(int)*sndPack[i].nints);
      sndPack[i].realData=(double *)malloc(sizeof(double)*sndPack[i].nreals);
      sndPack[i].intData[0]=obdonors[i];
      sndPack[i].intData[1]=obreceptors[i];
    }
  //
  // pack the data
  //
  intcount=(int *)malloc(sizeof(int)*nsend);
  realcount=(int *)malloc(sizeof(int)*nsend);
  for(i=0;i<nsend;i++) { intcount[i]=2; realcount[i]=0;}
  if (nblocks > 0) 
    {
     for(int ib=0;ib<nblocks;ib++)
     { 
      auto & mb = mblocks[ib];
      for(i=0;i<mb->ntotalPointsCart;i++) 
	{ 
	  if (mb->donorIdCart[i]!=-1)
	    {
	      gid=mb->donorIdCart[i];
	      procid=imap[cg->proc_id[gid]];
	      localid=cg->local_id[gid];	  
	      sndPack[procid].intData[intcount[procid]++]=localid;
	      sndPack[procid].intData[intcount[procid]++]=i;
	      sndPack[procid].intData[intcount[procid]++]=ib;
	      sndPack[procid].realData[realcount[procid]++]=mb->rxyzCart[3*i];
	      sndPack[procid].realData[realcount[procid]++]=mb->rxyzCart[3*i+1];
	      sndPack[procid].realData[realcount[procid]++]=mb->rxyzCart[3*i+2];	  
	    }
	}
      for(i=0;i<mb->nsearch;i++) 
	{
	  if (mb->donorId[i]!=-1) 
	    {
	      procid=imap[mb->isearch[3*i]];
	      sndPack[procid].intData[intcount[procid]++]=mb->isearch[3*i+1];
	      sndPack[procid].intData[intcount[procid]++]=mb->isearch[3*i+2];
	      sndPack[procid].intData[intcount[procid]++]=mb->meshtag;
	      sndPack[procid].intData[intcount[procid]++]=i;
	      sndPack[procid].intData[intcount[procid]++]=ib;
	      sndPack[procid].realData[realcount[procid]++]=mb->cellRes[mb->donorId[i]];
	      //      fprintf(fp,"%lf %lf %lf\n",mb->xsearch[3*i],mb->xsearch[3*i+1],mb->xsearch[3*i+2]);
	    }
	}
      }
    }
  //if (myid==0) {
  //for(i=0;i<nsend;i++)
  //  printf("intcount/intcount=%d %d\n",sndPack[i].nints,intcount[i]);
  //for(i=0;i<nsend;i++)
  //  printf("intcount/intcount=%d %d\n",sndPack[i].nreals,realcount[i]);
  //}
  TIOGA_FREE(realcount);
  //
  // exchange data
  //
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the data now
  //
  bcount=(int *)malloc(sizeof(int)*ncart);
  for(i=0;i<ncart;i++) { 
    cb[i].clearLists(); 
    cb[i].initializeLists();
    bcount[i]=0;
  }
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nreals > 0) 
	{
     	  m=2;
	  n=0;
	  interpCount=rcvPack[i].intData[0];
	  donorCount=rcvPack[i].intData[1];
	  icount[i]=interpCount+donorCount;
	  for(j=0;j<interpCount;j++)
	    {
     	      localid=rcvPack[i].intData[m++];
	      remoteid=rcvPack[i].intData[m++];
	      remoteblockid=rcvPack[i].intData[m++];
	      xtmp[0]=rcvPack[i].realData[n++];
	      xtmp[1]=rcvPack[i].realData[n++];
	      xtmp[2]=rcvPack[i].realData[n++];	      
	      cb[localid].insertInInterpList(i,remoteid,remoteblockid,xtmp);
	      bcount[localid]++;
	    }
	  for(j=0;j<donorCount;j++)
	    {
     	      localid=rcvPack[i].intData[m++];
	      index=rcvPack[i].intData[m++];
	      meshtag=rcvPack[i].intData[m++];
	      remoteid=rcvPack[i].intData[m++];
              remoteblockid=rcvPack[i].intData[m++];
	      cellRes=rcvPack[i].realData[n++];
 	      cb[localid].insertInDonorList(i,index,meshtag,remoteid,remoteblockid,cellRes);
	      bcount[localid]++;
	    }
	}
    }
  for(i=0;i<ncart;i++) cb[i].processDonors(holeMap,nmesh);
  pc_cart->clearPackets2(sndPack,rcvPack);  
  for(i=0;i<nsend;i++)
    {
      if (icount[i] > 0) 
	{
	  sndPack[i].nints=3*icount[i];
	  sndPack[i].nreals=0;
	  sndPack[i].intData=(int *) malloc(sizeof(int)*sndPack[i].nints);
	}
        intcount[i]=0;
    }
  cancelledData=NULL;
  //if (myid==30) {
  for(i=0;i<ncart;i++) 
    {
      //TRACEI(i);
      if (cancelledData) TIOGA_FREE(cancelledData);
      cancelledData=NULL;
      ncancel=bcount[i];
      if (ncancel > 0) {
	cancelledData=(int *)malloc(sizeof(int)*4*ncancel);
	cb[i].getCancellationData(cancelledData,&ncancel);
	for(j=0;j<ncancel;j++)
	  {
	    procid=cancelledData[4*j];
	    ctype=cancelledData[4*j+1];
	    remoteid=cancelledData[4*j+2];
            remoteblockid=cancelledData[4*j+3];
	    sndPack[procid].intData[intcount[procid]++]=ctype;
	    sndPack[procid].intData[intcount[procid]++]=remoteid;
	    sndPack[procid].intData[intcount[procid]++]=remoteblockid;
	    //TRACEI(intcount[procid]);
            //TRACEI(sndPack[procid].nints);
	  }
      }
    }
  for(i=0;i<nsend;i++) sndPack[i].nints=intcount[i];
  //}

  pc_cart->sendRecvPackets(sndPack,rcvPack);
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nints > 0) 
	{
	  m=0;
	  for(j=0;j<rcvPack[i].nints/3;j++)
	    {
	      ctype=rcvPack[i].intData[m++];
	      id=rcvPack[i].intData[m++];
              int ib = rcvPack[i].intData[m++];
	      if (ctype==0) 
		{
		  mblocks[ib]->donorIdCart[id]=-1;
		}
	      else
		{
		  mblocks[ib]->donorId[id]=-1;
		}
	    }
	}
    }

  for(int ib=0;ib<nblocks;ib++)
    mblocks[ib]->setCartIblanks();
  pc_cart->clearPackets2(sndPack,rcvPack);
  //
  for(int ib=0;ib<nblocks;ib++)
   mblocks[ib]->findInterpListCart();
  if (cancelledData) TIOGA_FREE(cancelledData);
  //fclose(fp);
  TIOGA_FREE(bcount);
  TIOGA_FREE(intcount);
  TIOGA_FREE(sndMapAll);
  TIOGA_FREE(rcvMapAll);
  TIOGA_FREE(imap);
  TIOGA_FREE(icount);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
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
  TIOGA_FREE(sndMap);
  TIOGA_FREE(rcvMap);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  printf("checkComm complete in %d\n",myid);
}
