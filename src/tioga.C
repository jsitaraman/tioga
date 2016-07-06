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
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

extern "C"{
//void writeqnode_(int *myid,double *qnodein,int *qnodesize);
}

/**
 * Set MPI communicator and initialize a few variables
 */
void tioga::setCommunicator(MPI_Comm communicator, int id_proc, int nprocs)
{
  scomm=communicator;
  myid=id_proc;
  numprocs=nprocs;
  sendCount=(int *) malloc(sizeof(int)*numprocs); /// Never used?
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
  mb->setData(btag,nnodes,xyz,ibl,nwbc,nobc,wbcnode,obcnode,ntypes,nv,nc,vconn);
  mb->myid=myid;
  mytag=btag;
}

void tioga::registerFaceConnectivity(int nftype, int *nf, int *nfv, int **fconn, int *f2c, int *c2f)
{
  mb->setFaceData(nftype, nf, nfv, fconn, f2c, c2f);
}

void tioga::profile(void)
{
  mb->preprocess();

  MPI_Allreduce(&iartbnd, &iabGlobal, 1, MPI_INT, MPI_MAX, scomm);
  MPI_Allreduce(&ihigh, &ihighGlobal, 1, MPI_INT, MPI_MAX, scomm);

  if (iabGlobal)
  {
    /// TODO: find a better location for this?
    // Split the global communicator into a mesh-local communicator
    MPI_Comm_split(scomm, mytag, myid, &meshcomm);
    MPI_Comm_rank(meshcomm, &meshRank);
    MPI_Comm_size(meshcomm, &nprocMesh);
  }
}

void tioga::performConnectivity(void)
{
  // Generate structured map of solid boundary (hole) locations
  getHoleMap();

  // Send/Recv the hole maps to/from all necessary ranks
  exchangeBoxes();

  // Find a list of all potential receptor points and send to all possible
  // donor ranks
  exchangeSearchData();

  // Find donors for all search points (other grids' possible receptor points)
  mb->search();

  // Exchange found donor data and do final iblank setting
  exchangeDonors();

  if (ihighGlobal)
  {
    // Calculate cell iblank values from nodal iblank values
    mb->getCellIblanks();

    if (iartbnd)  // Only done by ranks with high-order Artificial Boundary grids
    {
      // Find artificial boundary faces
      mb->calcFaceIblanks(meshcomm);

      // Get all AB face point locations
      mb->getBoundaryNodes();
    }
    else
    {
      // If this rank isn't high order, this will still load up fringe nodes
      mb->getInternalNodes();
    }

    // Exchange new list of points, including high-order Artificial Boundary
    // face points or internal points (or fringe nodes for non-high order)
    exchangePointSearchData();

    mb->search();

    // Setup interpolation weights and such for final interp-point list
    mb->processPointDonors();
  }
}

void tioga::performConnectivityHighOrder(void)
{
  mb->getInternalNodes();
  exchangePointSearchData();
  mb->search();
  mb->processPointDonors();
  iorphanPrint=1;
}

// !! PLACEHOLDERS - TO BE COMPLETED !!
void tioga::performConnectivityArtificialBoundary(void)
{
  mb->getBoundaryNodes();  //! Get all AB face point locations
  exchangePointSearchData();
  mb->search();
  mb->processPointDonors();
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
      mb->getCartReceptors(cg,pc_cart);
      mb->ihigh=ihigh;
      mb->search();
      mb->getUnresolvedMandatoryReceptors();
      cg->search(mb->rxyzCart,mb->donorIdCart,mb->ntotalPointsCart);
    }    
  //checkComm();
  exchangeAMRDonors();
  mb->getCellIblanks();
  mb->writeCellFile(myid);
  for(i=0;i<ncart;i++)
	cb[i].writeCellFile(i);
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

int tioga::findPointDonor(double *x_pt)
{
  return mb->findPointDonor(x_pt);
}

std::unordered_set<int> tioga::findCellDonors(double *bbox)
{
  return mb->findCellDonors(bbox);
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
  char ofname[10];
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
      if (itmp[i]==0) {
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
  mb->clearOrphans(itmp);
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

void tioga::dataUpdate_artBnd(int nvar, double *q_spts, double* q_fpts, int interpType)
{
  // initialize send and recv packets
  int nsend,nrecv;
  int *sndMap,*rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  if (nsend+nrecv == 0) return;

  PACKET *sndPack = new PACKET[nsend];
  PACKET *rcvPack = new PACKET[nrecv];

  std::vector<int> icount(nsend), dcount(nsend);

  pc->initPackets(sndPack,rcvPack);

  // get the interpolated solution now
  int *integerRecords = NULL;  // procID & pointID
  double *realRecords = NULL;  // Interpolated solution

  /// TODO: Replace 'interptype' with array of strides (simpler & more general)
  int nints, nreals;
  if (iartbnd)
  {
    interpType = 0;
    mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
                                        &realRecords,q_spts,nvar,interpType);
  }
  else if (ihigh && (ncart == 0))
  {
    mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
                                        &realRecords,q_spts,nvar,interpType);
  }
  else if (ncart > 0)
  {
    mb->getInterpolatedSolutionAMR(&nints,&nreals,&integerRecords,&realRecords,
                                   q_spts,nvar,interpType);
    for (int i = 0; i < ncart; i++)
      cb[i].getInterpolatedData(&nints,&nreals,&integerRecords,&realRecords,nvar);
  }
  else
  {
    // Note: same as original func, but using interpList2 instead
    mb->getInterpolatedSolution2(nints,nreals,integerRecords,realRecords,q_spts,
                                 nvar,interpType);
  }

  // Populate the packets [organize interp data by rank to send to]
  for (int i = 0; i < nints; i++)
  {
    int k = integerRecords[2*i]; // rank that interp point belongs to
    sndPack[k].nints++;
    sndPack[k].nreals += nvar;
  }

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double *)malloc(sizeof(double)*sndPack[k].nreals);
    icount[k] = dcount[k] = 0;
  }

  int m = 0;
  for (int i = 0; i < nints; i++)
  {
    int k = integerRecords[2*i];
    sndPack[k].intData[icount[k]++] = integerRecords[2*i+1];
    for (int j = 0; j < nvar; j++)
      sndPack[k].realData[dcount[k]++] = realRecords[m++];
  }

  // communicate the data across all partitions
  pc->sendRecvPackets(sndPack,rcvPack);

  // Decode the packets and update the values in the solver's data array
  if (ihigh)
  {
    std::vector<double> qtmp(nvar*mb->ntotalPoints);
    std::vector<int> itmp(mb->ntotalPoints);

    for (int k = 0; k < nrecv;k++)
    {
      int m = 0;
      for (int i = 0; i < rcvPack[k].nints; i++)
      {
        for (int j = 0; j < nvar; j++)
        {
          itmp[rcvPack[k].intData[i]] = 1;
          qtmp[rcvPack[k].intData[i]*nvar+j] = rcvPack[k].realData[m];
          m++;
        }
      }
    }

    // Print any 'orphan' points which may exist
    std::ofstream fp;
    int norphanPoint = 0;
    for (int i = 0; i < mb->ntotalPoints; i++)
    {
      if (itmp[i] == 0) {
        if (!fp.is_open())
        {
          std::stringstream ss;
          ss << "orphan" << myid << ".dat";
          fp.open(ss.str().c_str(),std::ofstream::out);
        }
        mb->outputOrphan(fp,i);
        norphanPoint++;
      }
    }
    fp.close();

    if (norphanPoint > 0 && iorphanPrint) {
      printf("Warning::number of orphans in %d = %d of %d\n",myid,norphanPoint,mb->ntotalPoints);
      iorphanPrint = 0;
    }

    // change the state of cells/nodes who are orphans
    mb->clearOrphans(itmp.data());

    if (iartbnd)
      mb->updateFluxPointData(q_fpts,qtmp.data(),nvar);
    else
      mb->updatePointData(q_spts,qtmp.data(),nvar,interpType);
  }
  else
  {
    for (int k = 0; k < nrecv; k++)
    {
      for (int i = 0; i < rcvPack[k].nints; i++)
      {
        mb->updateSolnData(rcvPack[k].intData[i],&rcvPack[k].realData[nvar*i],q_spts,nvar,interpType);
      }
    }
  }

  // release all memory
  pc->clearPackets2(sndPack,rcvPack);
  delete[] sndPack;
  delete[] rcvPack;
  delete[] integerRecords;
  delete[] realRecords;
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
