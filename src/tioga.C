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
#include "codetypes.h"
#include "tioga.h"
#include <assert.h>
using namespace TIOGA;
/**
 * set communicator
 * and initialize a few variables
 */
extern "C"{
//void writeqnode_(int *myid,double *qnodein,int *qnodesize);
};
void tioga::setCommunicator(MPI_Comm communicator, int id_proc, int nprocs)
{
  scomm=communicator;
  myid=id_proc;
  numprocs=nprocs;
  sendCount=(int *) malloc(sizeof(int)*numprocs);
  recvCount=(int *) malloc(sizeof(int)*numprocs);
  //
  // only one mesh block per process for now
  // this can be changed at a later date
  // but will be a fairly invasive change
  //
  // nblocks=0;
  // mb=new MeshBlock[1];
  //
  // instantiate the parallel communication class
  //
  pc=new parallelComm[1];
  pc->myid=myid;
  pc->scomm=scomm;
  pc->numprocs=numprocs;  
 
  // instantiate the parallel communication class
  //   
  pc_cart=new parallelComm[1];
  pc_cart->myid=myid;
  pc_cart->scomm=scomm;
  pc_cart->numprocs=numprocs;
  //
}
/**
 * register grid data for each mesh block
 */
void tioga::registerGridData(int btag,int nnodes,double *xyz,int *ibl, int nwbc,int nobc,
                             int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn,
                             uint64_t* cell_gid, uint64_t* node_gid)
{
  int iblk;

  auto idxit = tag_iblk_map.find(btag);
  if (idxit == tag_iblk_map.end()) {
    mtags.push_back(btag);
    mytag.push_back(btag);
    mblocks.push_back(std::unique_ptr<MeshBlock>(new MeshBlock));
    nblocks = mblocks.size();
    iblk = nblocks - 1;
    tag_iblk_map[btag] = iblk;
  } else {
    iblk = idxit->second;
  }

  auto& mb = mblocks[iblk];
  mb->setData(btag, nnodes, xyz, ibl, nwbc, nobc, wbcnode, obcnode, ntypes, nv,
              nc, vconn, cell_gid, node_gid);
  mb->myid = myid;
}

void tioga::registerSolution(int btag,double *q)
{
  auto idxit=tag_iblk_map.find(btag);
  int iblk=idxit->second;
  qblock[iblk]=q;
}

void tioga::profile(void)
{
  this->myTimer("tioga::profile",0);
  for(int ib=0;ib<nblocks;ib++)
   {
    auto& mb = mblocks[ib];
    mb->mexclude=mexclude;
    mb->nfringe=nfringe;
    mb->preprocess();
    //mb->writeGridFile(myid);
   }
  //mb->writeOBB(myid);
  //if (myid==4) mb->writeOutput(myid);
  //if (myid==4) mb->writeOBB(myid);
  this->myTimer("tioga::profile",1);
}

void tioga::performConnectivity(void)
{
  this->myTimer("tioga::performConnectivity",0);
  this->myTimer("tioga::getHoleMap",0);
  getHoleMap();
  this->myTimer("tioga::getHoleMap",1);
  this->myTimer("tioga::exchangeBoxes",0);
  exchangeBoxes();
  this->myTimer("tioga::exchangeBoxes",1);
  this->myTimer("tioga::exchangeSearchData",0);
  exchangeSearchData();
  this->myTimer("tioga::exchangeSearchData",1);
  this->myTimer("tioga::search",0);
  for(int ib=0;ib < nblocks;ib++)
  {
   auto& mb = mblocks[ib];
   mb->ihigh=0;
   mb->resetInterpData();
   mb->search();
  }
  this->myTimer("tioga::search",1);
  this->myTimer("tioga::exchangeDonors",0);
  exchangeDonors();
  this->myTimer("tioga::exchangeDonors",1);
  //this->reduce_fringes();
  outputStatistics();
  MPI_Allreduce(&ihigh,&ihighGlobal,1,MPI_INT,MPI_MAX,scomm);
  //if (ihighGlobal) {
  this->myTimer("tioga::getCellIblanks",0);
  for (int ib=0;ib<nblocks;ib++) {
    auto& mb = mblocks[ib];
    if (ihighGlobal) {
      mb->getCellIblanks2();
    }
    else {
      mb->getCellIblanks();
    }
    //mb->writeGridFile(100*myid+mtags[ib]);
  }
  this->myTimer("tioga::getCellIblanks",1);
  if (qblock) TIOGA_FREE(qblock);
  qblock=(double **)malloc(sizeof(double *)*nblocks);
  for(int ib=0;ib<nblocks;ib++)
    qblock[ib]=NULL;
  //}
  //mb->writeOutput(myid);
  //TRACEI(myid);
  //this->myTimer("tioga::performConnectivity",1);
}

void tioga::performConnectivityHighOrder(void)
{
 for(int ib=0;ib<nblocks;ib++)
 {
  auto& mb = mblocks[ib];
  mb->ihigh=ihigh;
  mb->getInternalNodes();
 }
 exchangeSearchData(1);
 for(int ib=0;ib<nblocks;ib++)
  { 
   auto& mb = mblocks[ib];
   mb->search();
   mb->processPointDonors();
  }
  iorphanPrint=1;
}  

void tioga::performConnectivityAMR(void)
{
  int i,ierr;
  int iamr;

  iamr=(ncart >0)?1:0;
  MPI_Allreduce(&iamr,&iamrGlobal,1,MPI_INT,MPI_MAX,scomm);
  cg->preprocess();
  for(i=0;i<ncart;i++) cb[i].preprocess(cg);
  
  if (nblocks > 0) 
    {
     for(int ib=0;ib<nblocks;ib++)
     {
      auto& mb = mblocks[ib];
      mb->getCartReceptors(cg,pc_cart);
      mb->ihigh=ihigh;
      mb->search();
      mb->getUnresolvedMandatoryReceptors();
      cg->search(mb->rxyzCart,mb->donorIdCart,mb->ntotalPointsCart);
     }
    }    
  //checkComm();
  exchangeAMRDonors();
  for(int ib=0;ib<nblocks;ib++)
   {
    auto &mb = mblocks[ib];
    mb->getCellIblanks();
   }
  //  mb->writeCellFile(myid);
  //for(i=0;i<ncart;i++)
	//cb[i].writeCellFile(i);
  MPI_Barrier(scomm);
  //printf("Finished performConnectivityAMR in %d\n",myid);
  //ierr=0;
  MPI_Abort(scomm,ierr);
}

void tioga::dataUpdate_AMR(int nvar,int interptype)
{
  int i,j,k,m;
  int nints;
  int nreals;
  int *integerRecords;
  double *realRecords;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;
  int *icount,*dcount;
  int bid;
  //
  // initialize send and recv packets
  //
  for(int ib=0;ib<nblocks;ib++)
    if (qblock[ib]==NULL) {
     printf("Solution data not set, cannot update \n");
     return;
    }   
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;
  //
  pc_cart->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);
  //
  pc_cart->initPackets(sndPack,rcvPack);  
  //
  // get the interpolated solution now
  //
  integerRecords=NULL;
  realRecords=NULL;
  //
  // TODO : verify for nblocks > 1
  //
  nints=nreals=0;
  for(int ib=0;ib<nblocks;ib++) {
   auto & mb = mblocks[ib];
   mb->getInterpolatedSolutionAMR(&nints,&nreals,&integerRecords,&realRecords,qblock[ib],nvar,interptype);
  }
  for(i=0;i<ncart;i++)
    cb[i].getInterpolatedData(&nints,&nreals,&integerRecords,&realRecords,nvar);
  //
  // populate the packets
  //
  for(i=0;i<nints;i++)
    {
      k=integerRecords[3*i];
      if (k <0 || k > nsend) {
	TRACEI(nsend);
	TRACEI(i);
	TRACEI(nints);
	TRACEI(k);
      }
      assert(k < nsend);      
      sndPack[k].nints+=2;
      sndPack[k].nreals+=nvar;
    }

  for(k=0;k<nsend;k++)
    {
     sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
     sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
     icount[k]=dcount[k]=0;
    }  

  m=0;
  for(i=0;i<nints;i++)
    {
      k=integerRecords[3*i];
      sndPack[k].intData[icount[k]++]=integerRecords[3*i+1];
      sndPack[k].intData[icount[k]++]=integerRecords[3*i+2];
      for(j=0;j<nvar;j++)
	sndPack[k].realData[dcount[k]++]=realRecords[m++];
    }
  //
  // communicate the data across
  //
  pc_cart->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the packets and update the data
  //
  for(k=0;k<nrecv;k++)
    {
      m=0;
      for(i=0;i<rcvPack[k].nints/2;i++)
	{
	  bid=rcvPack[k].intData[2*i];
	  if (bid < 0) 
	    {
              int tmp1=rcvPack[k].intData[2*i+1];
              int inode=mblocks[-(bid+1)]->receptorIdCart[tmp1];
	      mblocks[-(bid+1)]->updateSolnData(inode,&rcvPack[k].realData[m],qblock[-(bid+1)],nvar,interptype);
	    }
	  else
	    {
	      cb[bid-1].update(&rcvPack[k].realData[m],rcvPack[k].intData[2*i+1],nvar);
	    }
	    m+=nvar;
	}
    }
  //
  // release all memory
  //
  pc_cart->clearPackets2(sndPack,rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  if (integerRecords) TIOGA_FREE(integerRecords);
  if (realRecords) TIOGA_FREE(realRecords);
  if (icount) TIOGA_FREE(icount);
  if (dcount) TIOGA_FREE(dcount);
}

void tioga::dataUpdate(int nvar,int interptype, int at_points)
{
  int **integerRecords;
  double **realRecords;
  double **qtmp;
  int **itmp;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;
  char ofname[100];
  FILE *fp;
  //
  for(int ib=0;ib<nblocks;ib++)
    if (qblock[ib]==NULL) {
     printf("Solution data not set, cannot update \n");
     return;
    } 
  //
  // initialize send and recv packets
  //
  integerRecords=NULL;
  realRecords=NULL;
  itmp=NULL;
  qtmp=NULL;
  fp=NULL;
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);  
  //
  // get the interpolated solution now
  //
  integerRecords=(int **)malloc(sizeof(int*)*nblocks);
  realRecords=(double **)malloc(sizeof(double*)*nblocks);
  for(int ib=0;ib<nblocks;ib++)
    {
     integerRecords[ib]=NULL;
     realRecords[ib]=NULL;
    }
   
  std::vector<int> nints(nblocks,0), nreals(nblocks,0);
  std::vector<int> icount(nsend,0),dcount(nsend,0);

  for(int ib=0;ib<nblocks;ib++)
    {
      auto &mb =mblocks[ib];
      double *q  =qblock[ib];
      if (at_points==0) {
        mb->getInterpolatedSolution(&(nints[ib]),&(nreals[ib]),&(integerRecords[ib]),&(realRecords[ib]),
	     			  q,nvar,interptype);
      }
      else {
        mb->getInterpolatedSolutionAtPoints(&(nints[ib]),&(nreals[ib]),&(integerRecords[ib]),
                                      &(realRecords[ib]),q,nvar,interptype);
      }
      //
      // populate the packets
      //
      for(int i=0;i<nints[ib];i++)
	{
	  int k=integerRecords[ib][3*i];
	  sndPack[k].nints+=2;
	  sndPack[k].nreals+=nvar;
	}
    }
  //
  for(int k=0;k<nsend;k++)
    {
     sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
     sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
     icount[k]=dcount[k]=0;
    }  

  for(int ib=0;ib<nblocks;ib++)
    {
      int m=0;
      for(int i=0;i<nints[ib];i++)
	{
	  int k=integerRecords[ib][3*i];
	  sndPack[k].intData[icount[k]++]=integerRecords[ib][3*i+1];
	  sndPack[k].intData[icount[k]++]=integerRecords[ib][3*i+2];
          for(int j=0;j<nvar;j++)
	    sndPack[k].realData[dcount[k]++]=realRecords[ib][m++];
        }
    }
  //
  // communicate the data across
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the packets and update the data
  //
  if (at_points) {
   qtmp=(double **)malloc(sizeof(double*)*nblocks);
   itmp=(int **)malloc(sizeof(double*)*nblocks);
   for(int ib=0;ib<nblocks;ib++) 
    {
     qtmp[ib]=(double *)malloc(sizeof(double)*mblocks[ib]->ntotalPoints*nvar);
     itmp[ib]=(int *)malloc(sizeof(int)*mblocks[ib]->ntotalPoints);
     for(int i=0;i < mblocks[ib]->ntotalPoints;i++) itmp[ib][i]=0;
    }
  }
  //
  for(int k=0;k<nrecv;k++)
    {
      int l=0;
      int m=0;
      for(int i=0;i<rcvPack[k].nints/2;i++)
	{
	  int pointid=rcvPack[k].intData[l++];
	  int ib=rcvPack[k].intData[l++];
	  auto &mb = mblocks[ib];
          if (at_points==0) {
           double *q  = qblock[ib];
	   mb->updateSolnData(pointid,&rcvPack[k].realData[m],q,nvar,interptype);
          }
          else {
            for (int j=0;j<nvar;j++)
               qtmp[ib][pointid*nvar+j]=rcvPack[k].realData[m+j];
            itmp[ib][pointid]=1;
          }  
	  m+=nvar;
	}
    }
   int norphanPoint=0;
   int ntotalPoints=0;
   for(int ib=0;ib < nblocks;ib++)
    {
     auto &mb = mblocks[ib];
     for(int i=0;i<mb->ntotalPoints;i++)
      {
       if (itmp[ib][i]==0 && iorphanPrint) {
        if (fp==NULL) 
          {
            sprintf(ofname,"orphan%d.%d.dat",myid,ib);
            fp=fopen(ofname,"w");
          }
        mb->outputOrphan(fp,i);
        norphanPoint++;
       }
      } 
     ntotalPoints+=mb->ntotalPoints;
    }
  if (fp!=NULL) fclose(fp);
  if (norphanPoint > 0 && iorphanPrint) {
   printf("Warning::number of orphans in %d = %d of %d\n",myid,norphanPoint,
        ntotalPoints);
    iorphanPrint=0;
   }
  //
  // change the state of cells/nodes who are orphans
  //
  if (at_points) {
  for (int ib=0;ib<nblocks;ib++)
   { 
    auto &mb = mblocks[ib];
    mb->clearOrphans(holeMap,nmesh,itmp[ib]);
    mb->updatePointData(qblock[ib],qtmp[ib],nvar,interptype);
   }
  }
  //
  // release all memory
  //
  pc->clearPackets(sndPack,rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  if (integerRecords) {
    for(int ib=0;ib<nblocks;ib++) if (integerRecords[ib]) TIOGA_FREE(integerRecords[ib]);
    TIOGA_FREE(integerRecords);
  }
  if (realRecords) {
    for(int ib=0;ib<nblocks;ib++)
      if (realRecords[ib]) TIOGA_FREE(realRecords[ib]);
    TIOGA_FREE(realRecords);
  }
  if (qtmp) {
    for (int ib =0;ib<nblocks;ib++)
      if (qtmp[ib]) TIOGA_FREE(qtmp[ib]);
      TIOGA_FREE(qtmp);
  }
  if (itmp) {
    for(int ib=0;ib<nblocks;ib++)
     if (itmp[ib]) TIOGA_FREE(itmp[ib]);
    TIOGA_FREE(itmp);
  }
}

void tioga::writeData(int nvar,int interptype)
{
  //mb->writeGridFile(myid);
  for(int ib=0;ib<nblocks;ib++) {
     mblocks[ib]->writeFlowFile(100*myid+ib,qblock[ib],nvar,interptype);
     mblocks[ib]->writeCellFile(100*myid+ib);
  }
}

void tioga::getDonorCount(int btag, int *dcount,int *fcount)
{
  int nsend, nrecv;
  int *sndMap, *rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) {
      *dcount = 0;
      *fcount = 0;
      return;
  }

  auto idxit = tag_iblk_map.find(btag);
  int iblk=idxit->second;
  auto &mb = mblocks[iblk];
  mb->getDonorCount(dcount,fcount);
}

void tioga::getDonorInfo(int btag,int *receptors,int *indices,double *frac,int *dcount)
{
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  int i;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) return;

  auto idxit=tag_iblk_map.find(btag);
  int iblk=idxit->second;
  auto &mb = mblocks[iblk];
  mb->getDonorInfo(receptors,indices,frac);
  //
  // change to actual processor id here
  //
  for(i=0;i<4*(*dcount);i+=4)
    receptors[i]=sndMap[receptors[i]];
      
}

void tioga::getReceptorInfo(std::vector<int>& receptors)
{
  int nsend, nrecv;
  int *sndMap, *rcvMap;
  PACKET *sndPack, *rcvPack;
  pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);

  if (nsend == 0) {
      receptors.clear();
      return;
  }

  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack = (PACKET*)malloc(sizeof(PACKET) * nsend);
  rcvPack = (PACKET*)malloc(sizeof(PACKET) * nrecv);

  for (int i = 0; i < nsend; i++) {
    sndPack[i].nints = sndPack[i].nreals = 0;
    sndPack[i].intData = NULL;
    sndPack[i].realData = NULL;
  }
  //
  for (int i = 0; i < nrecv; i++) {
    rcvPack[i].nints = rcvPack[i].nreals = 0;
    rcvPack[i].intData = NULL;
    rcvPack[i].realData = NULL;
  }

  std::vector<int> dcount(nblocks), fcount(nblocks);
  std::vector<std::vector<int> > fringeSend(nblocks);

  for (int ib=0; ib < nblocks; ib++) {
    mblocks[ib]->getDonorCount(&dcount[ib], &fcount[ib]);

    // For each fringe the following data is returned by mesh block
    //
    // - MPI rank of the fringe point
    // - local node ID for the fringe point (on that proc)
    // - local mesh block index containing the fringe point (on that proc)
    // - Global ID of the donor cell
    //
    // The donor cell GID is uint64_t and is packed using 2*sizeof(int) bytes
    // through the int array. This makes total number of ints per fringe data
    // be (3 + 2) = 5
    //
    fringeSend[ib].resize(5 * dcount[ib]);
    mblocks[ib]->getReceptorInfo(fringeSend[ib].data());

    std::vector<int>& fringeData = fringeSend[ib];
    for (int i=0; i<(5*dcount[ib]); i+=5) {
      int k = fringeData[i];
      sndPack[k].nints += 4;
    }
  }

  // Allocate buffers for data send
  for (int k=0; k<nsend; k++) {
    sndPack[k].intData = (int*) malloc(sizeof(int) * sndPack[k].nints);
  }

  std::vector<int> ix(nsend, 0);
  for (int ib=0; ib<nblocks; ib++) {
    std::vector<int>& fringeData = fringeSend[ib];

    for(size_t i=0; i<fringeData.size(); i+=5) {
      int k = fringeData[i];
      sndPack[k].intData[ix[k]++] = fringeData[i+1];  // nodeID
      sndPack[k].intData[ix[k]++] = fringeData[i+2];  // local block index at receiver

      // The uint64_t donor cell ID data (transferred as 2 4-byte entries)
      sndPack[k].intData[ix[k]++] = fringeData[i+3];
      sndPack[k].intData[ix[k]++] = fringeData[i+4];
    }
  }

  pc->sendRecvPackets(sndPack, rcvPack);

  int rsize=0;
  for (int k=0; k<nrecv; k++) {
    rsize += rcvPack[k].nints;
  }
  receptors.resize(rsize);

  size_t idx=0;
  for (int k=0; k<nrecv; k++) {
    for (int j=0; j<rcvPack[k].nints; j+=4) {
      receptors[idx++] = rcvPack[k].intData[j];
      receptors[idx++] = mtags[rcvPack[k].intData[j+1]];
      receptors[idx++] = rcvPack[k].intData[j+2];
      receptors[idx++] = rcvPack[k].intData[j+3];
    }
  }

  pc->clearPackets(sndPack, rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
}

tioga::~tioga()
{      
  int i;
  if (holeMap)
    {
      for(i=0;i<nmesh;i++)
	if (holeMap[i].existWall) TIOGA_FREE(holeMap[i].sam);
      delete [] holeMap;
    }
  if (pc) delete[] pc;
  if (pc_cart) delete[] pc_cart;
  if (sendCount) TIOGA_FREE(sendCount);
  if (recvCount) TIOGA_FREE(recvCount);
  if (cb) delete [] cb;
  if (cg) delete [] cg;
  if (qblock) TIOGA_FREE(qblock);
  if (myid==0) printf("#tioga :successfully cleared all the memory accessed\n");
};

void tioga::register_amr_global_data(int nf,int qstride,double *qnodein,int *idata,
				     double *rdata,int ngridsin,
				     int qnodesize)
{
  if (cg) delete [] cg;
  cg=new CartGrid[1];
  cg->myid=myid;
  cg->registerData(nf,qstride,qnodein,idata,rdata,ngridsin,qnodesize);
  //writeqnode_(&myid,qnodein,&qnodesize);
}

void tioga::set_amr_patch_count(int npatchesin)
{
  ncart=npatchesin;
  if (cb) delete [] cb;
  cb=new CartBlock[ncart];
}

void tioga::register_amr_local_data(int ipatch,int global_id,int *iblank,double *q)
{
  cb[ipatch].registerData(ipatch,global_id,iblank,q);
}

#ifdef TIOGA_ENABLE_TIMERS
void tioga::myTimer(char const *info,int type)
#else
void tioga::myTimer(char const*, int)
#endif
{
#ifdef TIOGA_ENABLE_TIMERS
  static double t_start;
  double t_end;

  MPI_Barrier(scomm);
  if (type==0) {
    t_start=MPI_Wtime();
    if (myid==0) printf("Begin %s\n",info);
  }
  if (type==1) {
   t_end=MPI_Wtime();
   if (myid==0) printf("End %s, time taken=%lf\n",info,t_end-t_start);
  }
#endif
}

void tioga::reduce_fringes(void)
{
  //
  int nsend,nrecv;
  int *sndMap;
  int *rcvMap;
  PACKET *sndPack,*rcvPack;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  //if (nsend == 0) return;  

  for(int ib=0;ib<nblocks;ib++)
   {
    auto& mb = mblocks[ib];
    mb->reduce_fringes();
   }  
  //
  // now redo the process in exchangeDonors to reduce
  // the data
  //
  // TODO make the stuff below into a subroutine to be more modular
  //
  // Find cancellation data (based on donor quality)
  //
  std::vector<int> nrecords(nblocks,0);
  int** donorRecords = (int**)malloc(sizeof(int*)*nblocks);
  //
  for (int i=0; i<nblocks; i++) {
    donorRecords[i]=NULL;
    nrecords[i] = 0;
    mblocks[i]->getCancellationData(&(nrecords[i]),
                                    &(donorRecords[i]));
  }
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);  
  //
  std::vector<int> nintsSend(nsend,0);
  std::vector<int> ixOffset(nsend,0);  
  std::fill(nintsSend.begin(), nintsSend.end(), 0);
  std::fill(ixOffset.begin(), ixOffset.end(), 0);
  //
  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].nints+=2;
    }
  }
  //
  for(int k=0;k<nsend;k++)
    sndPack[k].intData = (int*)malloc(sizeof(int) * sndPack[k].nints);
  //
  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k= donorRecords[n][3*i];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+1];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+2];
    }
  }
  //
  // communciate cancellation data comm 3
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for (int k=0; k<nrecv; k++) {
    int m = 0;
    for (int j=0;j<rcvPack[k].nints/2;j++) {
      int recid=rcvPack[k].intData[m++];
      int ib = tag_iblk_map[rcvPack[k].intData[m++]];
      mblocks[ib]->cancelDonor(recid);
    }
  }
  //
  for (int ib=0;ib<nblocks;ib++) {
    auto &mb =mblocks[ib];
    mb->resetCoincident();
  }
  //
  pc->clearPackets(sndPack, rcvPack);
  //
  // Find final interpolation data
  //
  for (int i=0; i<nblocks; i++) {
    if (donorRecords[i]) {
      TIOGA_FREE(donorRecords[i]);
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
      int k = donorRecords[n][3*i];
      sndPack[k].nints+=2;
    }
  }
  for(int k=0;k<nsend;k++)
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
  for (int n=0; n < nblocks; n++) {
    for (int i=0; i<nrecords[n]; i++) {
      int k = donorRecords[n][3*i];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+1];
      sndPack[k].intData[ixOffset[k]++]=donorRecords[n][3*i+2];
    }
  }
  //
  // comm 4
  // final receptor data to set iblanks
  //     
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for(int ib=0;ib<nblocks;ib++)
    mblocks[ib]->clearIblanks();
  
  for (int k=0; k<nrecv; k++) {
    int m = 0;
    for(int j=0;j< rcvPack[k].nints/2;j++)
      {
	int pointid=rcvPack[k].intData[m++];
	int ib=rcvPack[k].intData[m++];
	mblocks[ib]->setIblanks(pointid);
      }
  }
  pc->clearPackets(sndPack,rcvPack);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
  
  if (donorRecords) {
    for (int i=0; i<nblocks; i++) {
      if (donorRecords[i]) TIOGA_FREE(donorRecords[i]);
    }
    TIOGA_FREE(donorRecords);
  }
  outputStatistics();
  //mb->writeOBB(myid);
  //if (myid==4) mb->writeOutput(myid);
  //if (myid==4) mb->writeOBB(myid);
}
