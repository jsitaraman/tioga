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

#ifdef _GPU
#include "cuda_funcs.h"
#endif

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

void tioga::registerFaceConnectivity(int nftype, int *nf, int *nfv, int **fconn,
                                     int *f2c, int *c2f, int *iblank_face,
                                     int nOverFaces, int nMpiFaces,
                                     int *overFaces, int *mpiFaces,
                                     int *mpiProcR, int *mpiFidR)
{
  mb->setFaceData(nftype, nf, nfv, fconn, f2c, c2f, iblank_face, nOverFaces,
                  nMpiFaces, overFaces, mpiFaces, mpiProcR, mpiFidR);
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
//  // Generate structured map of solid boundary (hole) locations
//  getHoleMap();

  // Send/Recv the hole maps to/from all necessary ranks
  exchangeBoxes();

//  // Find a list of all potential receptor points and send to all possible
//  // donor ranks
//  exchangeSearchData();

//  // Find donors for all search points (other grids' possible receptor points)
//  mb->search();

//  // Exchange found donor data and do final iblank setting
//  exchangeDonors();

  if (ihighGlobal)
  {
    // Calculate cell iblank values from nodal iblank values
    mb->getCellIblanks(meshcomm);

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

#ifdef _GPU
    setupCommBuffersGPU();
#endif
  }
}

#ifdef _GPU
void tioga::setupCommBuffersGPU(void)
{
  // Setup additional index arrars for GPU to automatically pack MPI buffers
  int nsend, nrecv, *sndMap, *rcvMap;
  pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);

  sndVPack.resize(nsend);
  rcvVPack.resize(nrecv);

  pc->initPacketsV(sndVPack,rcvVPack);

  mb->setupBuffersGPU(nsend, intData, sndVPack);

  ninterp = mb->ninterp2;
}
#endif

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
//  /* --- Direct Cut Method [copied from highOrder.C] --- */
//  int nDims = 3; /// TODO: add to MB.h ...

//  /// TODO: pre-process groups on this rank vs. groups on ALL ranks/grids
//  int maxID = 0;
//  for (int g = 0; g < nGroups; g++)
//    maxID = max(maxID, groupIDs[g]);

//  int nGroups_glob = maxID;
//  MPI_Allreduce(&maxId, &nGroups_glob, 1, MPI_INT, MPI_MAX, scomm);

//  // Full list of group types, and list of group IDs for each grid
//  std::vector<int> cutType_glob(nGroups_glob, -1);
//  std::vector<int> gridGroups(nGrids*nGroups_glob);
//  std::unordered_set<int> myGroups;

//  for (int g = 0; g < nGroups; g++)
//  {
//    int G = groupIDs[g];
//    cutType_glob[G] = cutType[g];
//    gridGroups[gridID*nGroups_glob+G] = 1;
//    myGroups.insert(G);
//  }

//  MPI_Allreduce(&cutType_glob[0], MPI_IN_PLACE, nGroups_glob, MPI_INT, MPI_MAX, scomm);
//  MPI_Allreduce(&gridGroups[0], MPI_IN_PLACE, nGrids*nGroups_glob, MPI_INT, MPI_MAX, scomm);


//  // Get cutting-group bounding-box data for this rank
//  std::vector<double> cutBox(6*nGroups_glob);
//  std::vector<std::vector<double>> faceBox(nGroups);
//  std::vector<std::vector<double>> faceNodes(nGroups);

//  mb->getCutGroupBoxes(cutBox, faceBox, nGroups_glob);
//  mb->getCutGroupFaces(faceNodes, nGroups_glob);

//  // Get the global bounding box info across all the partitions for all meshes
//  std::vector<double> cutBox_global(6*nGroups_glob);
//  for (int G = 0; G < nGroups_glob; G++)
//  {
//    MPI_Allreduce(&cutBox[6*G],  &cutBox_global[6*G],  3,MPI_DOUBLE,MPI_MIN,scomm);
//    MPI_Allreduce(&cutBox[6*G+3],&cutBox_global[6*G+3],3,MPI_DOUBLE,MPI_MAX,scomm);
//  }

//  // Figure out which cutting groups overlap with this rank
//  // [use rank-local OBB?]
//  // send/recv facebox data to/from ranks that need it
//  std::vector<std::unordered_set<int>> cellList(nGroups_glob);
//  mb->getDirectCutCells(cellList, cutBox_global, nGroups_global);

//  std::vector<int> nGFaces(nGroups_glob);
//  for (int g = 0; g < nGroups; g++)
//  {
//    nGFaces[groupIDs[g]] = nGf[g];
//  }

//  // Send all face box  data to all ranks?
//  std::vector<std::vector<int>> nFaces_glob(nGroups_glob);
//  std::vector<std::vector<double>> faceBox_glob(nGroups_glob);
//  for (int G = 0; G < nGroups_glob; G++)
//  {
//    nFaces_glob[G].resize(nproc);
//    MPI_Allgather(&nGFaces[G], 1, MPI_INT, &nFaces_glob[G][0], 1, MPI_INT, scomm);

//    std::vector<int> recvCnts(nproc);
//    std::vector<int> recvDisp(nproc);
//    int nface_glob = 0;
//    for (int p = 0; p < nproc; p++)
//    {
//      recvCnts[p] = nFaces_glob[G][p]*6;
//      recvDisp[p] += (p>0) ? nFaces_glob[G][p-1]*6 : 0;
//      nface_glob += nFaces_glob[G][p];
//    }

//    faceBox_glob[G].resize(nface_glob*6);
//    MPI_Allgatherv(&faceBox[G][0], nGFces[G], MPI_DOUBLE, &faceBox_glob[G][0], recvCnts.data(), recvDisp.data(), MPI_DOUBLE, scomm);
//  }

//  /// first try
////  // use 'PC' object to send/recv faceBox data to correct ranks
////  // Ranks for which cellList[G] != 0 must send data
////  int nsend = 0;
////  int nrecv = 0;
////  std::vector<int> sendMap, sendInts, recvMap, recvGroups;
////  for (int p = 0; p < numprocs; p++)
////  {
////    if (p == rank) continue;

////    int grid = gridIds[p];
////    for (int G = 0; G < nGroups_glob; G++)
////    {
////      if (gridGroups[grid*nGroups_glob+G]) // && cellList[G].size() > 0)
////      {
////        nrecv++;
////        nsend++;
////        recvMap.push_back(p);
////        recvGroups.push_back(G);
////        sendMap.push_back(p);
////        sendInts.push_back(cellList[G].size());
////      }
////    }
////  }
  
//  // Do the cutting
//  for (int G = 0; G < nGroups_glob; G++)
//  {
//    mb->directCut(faceBox_glob[G]);

//  }


//  mb->getBoundaryNodes();  //! Get all AB face point locations
//  exchangePointSearchData();
//  mb->search();
//  mb->processPointDonors();
//  iorphanPrint=1;
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
  mb->getCellIblanks(meshcomm);
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
  waitTime.showTime();
  interpTime.showTime();
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
   printf("Warning::number of orphans in rank %d = %d of %d\n",myid,norphanPoint,
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

#ifdef _GPU
void tioga::dataUpdate_artBnd(int nvar, double *q_spts, int dataFlag)
{
  // initialize send and recv packets
  int nsend,nrecv;
  int *sndMap,*rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  if (nsend+nrecv == 0) return;

  // get the interpolated solution now
  int stride = nvar;
  if (iartbnd && dataFlag == 1) stride = 3*nvar;

  /// TODO: consider placement of this
  // Allocate space for the interpolated fringe-point data if needed
  if (mb->ninterp2 > ninterp_d)
  {
    if (ninterp_d > 0)
      cuda_free(ubuf_d);

    resizeFlag = 1;
    ninterp_d = mb->ninterp2;
    cuda_malloc(ubuf_d, ninterp_d*stride);
  }
  else if (dataFlag && resizeFlag)
  {
    if (gradbuf_d)
      cuda_free(gradbuf_d);

    cuda_malloc(gradbuf_d, ninterp_d*stride);
    resizeFlag = 0;
  }

  int nints, nreals;
  interpTime.startTimer();

  /// TODO: make everything asyc but sync here to wait on corrected gradient
  nints = ninterp_d;
  nreals = stride*ninterp_d;
  if (dataFlag == 0)
    mb->interpSolution_gpu(ubuf_d, nvar);
  else
    mb->interpGradient_gpu(gradbuf_d, nvar);

  for (int k = 0; k < nsend; k++)
  {
    sndVPack[k].nreals = sndVPack[k].nints * stride;
    sndVPack[k].realData.resize(sndVPack[k].nreals);
  }

  dblData.resize(ninterp_d*stride);

  if (dataFlag == 0)
    cuda_copy_d2h(ubuf_d, dblData.data(), ninterp_d*stride);
  else
    cuda_copy_d2h(gradbuf_d, dblData.data(), ninterp_d*stride);

  interpTime.stopTimer();

  /// TODO: *** Split function here (like 'send_u_data', 'recv_u_data') ***

  PUSH_NVTX_RANGE("tg_packBuffers", 3);
  // Populate the packets [organize interp data by rank to send to]
  for (int p = 0; p < nsend; p++)
  {
    double *ptr = dblData.data()+mb->buf_disp[p]*stride;
    sndVPack[p].realData.assign(ptr,ptr+sndVPack[p].nints*stride);
  }
  POP_NVTX_RANGE;

  // communicate the data across all partitions
  PUSH_NVTX_RANGE("tg_pc_sendRecv", 0);
  MPI_Pcontrol(1, "tioga_dataUpdate_pc");
  pc->sendRecvPacketsV(sndVPack,rcvVPack);
  MPI_Pcontrol(-1, "tioga_dataUpdate_pc");
  POP_NVTX_RANGE;

  // Decode the packets and update the values in the solver's data array
  PUSH_NVTX_RANGE("tg_unpack_data", 1);
  if (ihigh)
  {
    MPI_Pcontrol(1, "tioga_unpack_buffers");
    std::vector<double> qtmp(stride*mb->ntotalPoints);
    std::vector<int> itmp(mb->ntotalPoints);

    for (int k = 0; k < nrecv;k++)
    {
      int m = 0;
      for (int i = 0; i < rcvVPack[k].nints; i++)
      {
        for (int j = 0; j < stride; j++)
        {
          itmp[rcvVPack[k].intData[i]] = 1;
          qtmp[rcvVPack[k].intData[i]*stride+j] = rcvVPack[k].realData[m];
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
      printf("Warning::number of orphans in rank %d = %d of %d\n",myid,norphanPoint,mb->ntotalPoints);
      iorphanPrint = 0;
    }

    // change the state of cells/nodes who are orphans
    if (!iartbnd)
      mb->clearOrphans(itmp.data());

    MPI_Pcontrol(-1, "tioga_unpack_buffers");

    if (iartbnd)
    {
      interpTime.startTimer();
      if (dataFlag == 0)
        mb->updateFluxPointData(qtmp.data(),nvar);
      else
        mb->updateFluxPointGradient(qtmp.data(),nvar);
      interpTime.startTimer();
    }
    else
      mb->updatePointData(q_spts,qtmp.data(),nvar,dataFlag);
  }
  POP_NVTX_RANGE;

  /// TODO: modify for moving grids
  // release all memory
  //pc->clearPacketsV(sndVPack, rcvVPack);
}
#endif

#ifndef _GPU
void tioga::dataUpdate_artBnd(int nvar, double *q_spts, int dataFlag)
{
  if (iartbnd && gpu)
    mb->getDonorDataGPU(dataFlag);

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

  int stride = nvar;
  if (iartbnd && dataFlag == 1) stride = 3*nvar;

  int nints, nreals;
  if (iartbnd)
  {
    if (dataFlag == 0)
      mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
                                          &realRecords,q_spts,nvar,dataFlag);
    else
      mb->getInterpolatedGradientAtPoints(nints,nreals,integerRecords,
                                          realRecords,q_spts,nvar);
  }
  else if (ihigh && (ncart == 0))
  {
    mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
                                        &realRecords,q_spts,nvar,dataFlag);
  }
  else if (ncart > 0)
  {
    mb->getInterpolatedSolutionAMR(&nints,&nreals,&integerRecords,&realRecords,
                                   q_spts,nvar,dataFlag);
    for (int i = 0; i < ncart; i++)
      cb[i].getInterpolatedData(&nints,&nreals,&integerRecords,&realRecords,nvar);
  }
  else
  {
    // Note: same as original func, but using interpList2 instead
    mb->getInterpolatedSolution2(nints,nreals,integerRecords,realRecords,q_spts,
                                 nvar,dataFlag);
  }

  // Populate the packets [organize interp data by rank to send to]
  for (int i = 0; i < nints; i++)
  {
    int k = integerRecords[2*i]; // rank that interp point belongs to
    sndPack[k].nints++;
    sndPack[k].nreals += stride;
  }

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double *)malloc(sizeof(double)*sndPack[k].nreals);
    icount[k] = dcount[k] = 0;
  }

  for (int i = 0; i < nints; i++)
  {
    int k = integerRecords[2*i];
    sndPack[k].intData[icount[k]++] = integerRecords[2*i+1];
    for (int j = 0; j < stride; j++)
      sndPack[k].realData[dcount[k]++] = realRecords[stride*i+j];
  }

  // communicate the data across all partitions
  pc->sendRecvPackets(sndPack,rcvPack);

  // Decode the packets and update the values in the solver's data array
  if (ihigh)
  {
    MPI_Pcontrol(1, "tioga_unpack_buffers");
    std::vector<double> qtmp(stride*mb->ntotalPoints);
    std::vector<int> itmp(mb->ntotalPoints);

    for (int k = 0; k < nrecv;k++)
    {
      int m = 0;
      for (int i = 0; i < rcvPack[k].nints; i++)
      {
        for (int j = 0; j < stride; j++)
        {
          itmp[rcvPack[k].intData[i]] = 1;
          qtmp[rcvPack[k].intData[i]*stride+j] = rcvPack[k].realData[m];
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
      printf("Warning::number of orphans in rank %d = %d of %d\n",myid,norphanPoint,mb->ntotalPoints);
      iorphanPrint = 0;
    }

    // change the state of cells/nodes who are orphans
    if (!iartbnd)
      mb->clearOrphans(itmp.data());

    if (iartbnd)
    {
      interpTime.startTimer();
      if (dataFlag == 0)
        mb->updateFluxPointData(qtmp.data(),nvar);
      else
        mb->updateFluxPointGradient(qtmp.data(),nvar);
      interpTime.startTimer();
    }
    else
      mb->updatePointData(q_spts,qtmp.data(),nvar,dataFlag);
  }
  else
  {
    for (int k = 0; k < nrecv; k++)
    {
      for (int i = 0; i < rcvPack[k].nints; i++)
      {
        mb->updateSolnData(rcvPack[k].intData[i],&rcvPack[k].realData[nvar*i],q_spts,nvar,dataFlag);
      }
    }
  }

  // release all memory
  pc->clearPackets2(sndPack,rcvPack);
  delete[] sndPack;
  delete[] rcvPack;
  free(integerRecords);
  free(realRecords);

  if (iartbnd && gpu)
    mb->sendFringeDataGPU(dataFlag);
}
#endif

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
