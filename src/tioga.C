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
#ifdef USE_CUDA
#include "cuda_functions.h"
#endif
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
  nblocks=0;
  mb=new MeshBlock[1];
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
			     int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn)
{
  if (nnodes > 0) nblocks=1;

  mb->resolutionScale=1.0 - btag*1.0e-10;
  // mb->resolutionScale=1000.0;

  mb->setData(btag,nnodes,xyz,ibl,nwbc,nobc,wbcnode,obcnode,ntypes,nv,nc,vconn);
  mb->myid=myid;
  mytag=btag;
}

void tioga::profile(void)
{

  mb->preprocess();
  //mb->writeGridFile(myid);
  //mb->writeOBB(myid);
  //if (myid==4) mb->writeOutput(myid);
  //if (myid==4) mb->writeOBB(myid);
}

void tioga::performConnectivity(void)
{
  getHoleMap();
  exchangeBoxes();
  exchangeSearchData();
  mb->ihigh=0;
  mb->search();
  exchangeDonors();
  outputStatistics();
  MPI_Allreduce(&ihigh,&ihighGlobal,1,MPI_INT,MPI_MAX,scomm);
  //if (ihighGlobal) {
  mb->getCellIblanks();
  // mb->writeCellFile(myid);
  //}
  //mb->writeOutput(myid);
  //tracei(myid);
}

void tioga::performConnectivityHighOrder(void)
{
  mb->ihigh=ihigh;
  mb->getInternalNodes();
  exchangePointSearchData();
  mb->search();
  mb->processPointDonors();
  iorphanPrint=0;
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
      mb->getCartReceptors(cg,pc_cart);
      mb->ihigh=ihigh;
      mb->search();
      mb->getUnresolvedMandatoryReceptors();
      cg->search(mb->rxyzCart,mb->donorIdCart,mb->ntotalPointsCart);
    }    
  //checkComm();
  exchangeAMRDonors();
  mb->getCellIblanks();
  //mb->writeCellFile(myid);
  //for(i=0;i<ncart;i++)
	//cb[i].writeCellFile(i);
  //MPI_Barrier(scomm);
  //printf("Finished performConnectivityAMR in %d\n",myid);
  //ierr=0;
  //MPI_Abort(scomm,ierr);
}

void tioga::dataUpdate_AMR(int nvar,double *q,int interptype)
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
  mb->getInterpolatedSolutionAMR(&nints,&nreals,&integerRecords,&realRecords,q,nvar,interptype);
  for(i=0;i<ncart;i++)
    cb[i].getInterpolatedData(&nints,&nreals,&integerRecords,&realRecords,nvar);
  //
  // populate the packets
  //
  for(i=0;i<nints;i++)
    {
      k=integerRecords[3*i];
      if (k <0 || k > nsend) {
	tracei(nsend);
	tracei(i);
	tracei(nints);
	tracei(k);
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
	  if (bid< 0) 
	    {
	      mb->updateSolnData(rcvPack[k].intData[2*i+1],&rcvPack[k].realData[m],q,nvar,interptype);
	    }
	  else
	    {

	      cb[bid].update(&rcvPack[k].realData[m],rcvPack[k].intData[2*i+1],nvar);
	      m+=nvar;
	    }
	}
    }
  //
  // release all memory
  //
  pc_cart->clearPackets2(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);
}


void tioga::dataUpdate(int nvar,double *q,int interptype)
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
  //
  // initialize send and recv packets
  //
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);  
  //
  // get the interpolated solution now
  //
  integerRecords=NULL;
  realRecords=NULL;
  mb->getInterpolatedSolution(&nints,&nreals,&integerRecords,&realRecords,q,nvar,interptype);
  //
  // populate the packets
  //
  for(i=0;i<nints;i++)
    {
      k=integerRecords[2*i];
      sndPack[k].nints++;
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
      k=integerRecords[2*i];
      sndPack[k].intData[icount[k]++]=integerRecords[2*i+1];
      for(j=0;j<nvar;j++)
	sndPack[k].realData[dcount[k]++]=realRecords[m++];
    }
  //
  // communicate the data across
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the packets and update the data
  //
  for(k=0;k<nrecv;k++)
    {
      m=0;
      for(i=0;i<rcvPack[k].nints;i++)
	{
	  mb->updateSolnData(rcvPack[k].intData[i],&rcvPack[k].realData[m],q,nvar,interptype);
	  m+=nvar;
	}
    }
  //
  // release all memory
  //
  pc->clearPackets2(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);
}

#ifdef USE_CUDA
void tioga::dataUpdate(GPUvec<double> *vec)
{
  int i,j,k,m;
  int nints;
  int nreals;
  int *integerRecords = NULL;
  double *realRecords = NULL;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;
  int *icount,*dcount;

  int interptype = 0;
  int nvar  = vec->nvar;
  // double* q = new double[vec->nvar*vec->pts];
  // vec->to_cpu(q);

  // printf("CALLING GPU VERSION OF DATAUPDATE\n");

  //
  // initialize send and recv packets
  //
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);  
  //
  // get the interpolated solution now
  //
  integerRecords=NULL;
  realRecords=NULL;

  // mb->getInterpolatedSolution(&nints,&nreals,&integerRecords,&realRecords,q,nvar,interptype);
  mb->getInterpolatedSolution(&nints,&nreals,&integerRecords,&realRecords,vec);

  //
  // populate the packets
  //
  for(i=0;i<nints;i++)
    {
      k=integerRecords[2*i];
      sndPack[k].nints++;
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
      k=integerRecords[2*i];
      sndPack[k].intData[icount[k]++]=integerRecords[2*i+1];
      for(j=0;j<nvar;j++)
	sndPack[k].realData[dcount[k]++]=realRecords[m++];
    }

  //
  // communicate the data across
  //
  pc->sendRecvPackets(sndPack,rcvPack);

  //
  // decode the packets and update the data
  //
  // for(k=0;k<nrecv;k++)
  //   {
  //     m=0;
  //     for(i=0;i<rcvPack[k].nints;i++)
  // 	{
  // 	  printf("_c_ q[%d] = q[%d]\n", rcvPack[k].intData[i]*nvar+k, m);
  // 	  mb->updateSolnData(rcvPack[k].intData[i],&rcvPack[k].realData[m],q,nvar,interptype);
  // 	  m+=nvar;
  // 	}
  //   }
  updateSolnGPU(nrecv, rcvPack, vec);

  //
  // release all memory
  //
  pc->clearPackets2(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);

  // vec->to_gpu(q);
  // delete q;
}
#endif

void tioga::writeData(int nvar,double *q,int interptype)
{
  //mb->writeGridFile(myid);
  mb->writeFlowFile(myid,q,nvar,interptype);
}

void tioga::getDonorCount(int *dcount,int *fcount)
{
  mb->getDonorCount(dcount,fcount);
}

void tioga::getDonorInfo(int *receptors,int *indices,double *frac,int *dcount)
{
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  int i;
  mb->getDonorInfo(receptors,indices,frac);
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap); 
  //
  // change to actual processor id here
  //
  for(i=0;i<3*(*dcount);i+=3)
    receptors[i]=sndMap[receptors[i]];
      
}

tioga::~tioga()
{      
  int i;
  if (mb) delete[] mb;
  if (holeMap) 
    {
      for(i=0;i<nmesh;i++)
	if (holeMap[i].existWall) free(holeMap[i].sam);
      delete [] holeMap;
    }
  if (pc) delete[] pc;
  if (sendCount) free(sendCount);
  if (recvCount) free(recvCount);
  if (obblist) free(obblist);
  if (myid==0) printf("#tioga :successfully cleared all the memory accessed\n");
};

void tioga::dataUpdate_highorder(int nvar,double *q,int interptype)
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
  double *qtmp;
  int norphanPoint;
  int *itmp;
  FILE *fp;
  char ofname[100];
  //
  // initialize send and recv packets
  //
  fp=NULL;
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);
  //
  pc->initPackets(sndPack,rcvPack);  
  //
  // get the interpolated solution now
  //
  integerRecords=NULL;
  realRecords=NULL;
  mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
				      &realRecords,q,nvar,interptype);
  //
  // populate the packets
  //
  for(i=0;i<nints;i++)
    {
      k=integerRecords[2*i];
      sndPack[k].nints++;
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
      k=integerRecords[2*i];
      sndPack[k].intData[icount[k]++]=integerRecords[2*i+1];
      for(j=0;j<nvar;j++)
	sndPack[k].realData[dcount[k]++]=realRecords[m++];
    }
  //
  // communicate the data across
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // decode the packets and update the data
  //
  qtmp=(double *)malloc(sizeof(double)*nvar*mb->ntotalPoints);
  itmp=(int *) malloc(sizeof(int)*mb->ntotalPoints);
  for(i=0;i<mb->ntotalPoints;i++) itmp[i]=0;
  //
  for(k=0;k<nrecv;k++)
    {
      m=0;
      for(i=0;i<rcvPack[k].nints;i++)
	{
	  for(j=0;j<nvar;j++)
	    {
	      itmp[rcvPack[k].intData[i]]=1;
	      qtmp[rcvPack[k].intData[i]*nvar+j]=rcvPack[k].realData[m];
	      m++;
	    }
         //if (myid==12) printf("rcvPack[k].intData[i]=%d %d\n",rcvPack[k].intData[i],mb->ntotalPoints);
	}
    }
  norphanPoint=0;
  for(i=0;i<mb->ntotalPoints;i++)
    {
      if (itmp[i]==0 && iorphanPrint) {
	if (fp==NULL) 
	  {
	    sprintf(ofname,"orphan%d.dat",myid);
	    fp=fopen(ofname,"w");
	  }
	mb->outputOrphan(fp,i);
	norphanPoint++;
      }
    }
  if (fp!=NULL) fclose(fp);
  if (norphanPoint > 0 && iorphanPrint) {
   printf("Warning::number of orphans in %d = %d of %d\n",myid,norphanPoint,
	mb->ntotalPoints);
    iorphanPrint=0;
   }
  // change the state of cells/nodes who are orphans
  mb->clearOrphans(holeMap,nmesh,itmp);
  //
  mb->updatePointData(q,qtmp,nvar,interptype);
  //
  // release all memory
  //
  pc->clearPackets2(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  free(qtmp);
  free(itmp);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);
}

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
