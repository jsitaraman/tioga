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
  nproc=nprocs;
  sendCount=(int *) malloc(sizeof(int)*nproc); /// Never used?
  recvCount=(int *) malloc(sizeof(int)*nproc);
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
  pc->numprocs=nproc;
 
  // instantiate the parallel communication class
  //   
  pc_cart=new parallelComm[1];
  pc_cart->myid=myid;
  pc_cart->scomm=scomm;
  pc_cart->numprocs=nproc;

#ifdef _GPU
//  mb_d = new dMeshBlock(); /// TODO
#endif
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

void tioga::registerFaceConnectivity(int gtype, int nftype, int *nf, int *nfv,
    int **fconn, int *f2c, int *c2f, int *iblank_face, int nOverFaces,
    int nWallFaces, int nMpiFaces, int *overFaces, int *wallFaces, int *mpiFaces, 
    int *mpiProcR, int *mpiFidR)
{
  mb->setFaceData(gtype, nftype, nf, nfv, fconn, f2c, c2f, iblank_face, nOverFaces,
                  nWallFaces, nMpiFaces, overFaces, wallFaces, mpiFaces, mpiProcR, mpiFidR);
}

#ifdef _GPU
void tioga::registerDeviceGridData(double *xyz, double *coords, int *ibc, int *ibf)
{
  mb->setDeviceData(xyz,coords,ibc,ibf);
}
#endif

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

    gridType = mb->gridType;

    // For the direct-cut method
    // gridType == 0: background;  gridType == 1: body grid
    gridIDs.resize(nproc);
    gridTypes.resize(nproc);
    MPI_Allgather(&mytag, 1, MPI_INT, gridIDs.data(), 1, MPI_INT, scomm);
    MPI_Allgather(&gridType, 1, MPI_INT, gridTypes.data(), 1, MPI_INT, scomm);

    nGrids = 0;
    for (auto g : gridIDs)
      nGrids = std::max(nGrids, g);
    nGrids++; // 0-index to size

    // Setup send/recv map based on grid ID
    std::vector<int> ptags(nproc);
    MPI_Allgather(&mytag, 1, MPI_INT, ptags.data(), 1, MPI_INT, scomm);

    // Enforcing symmetry in sends/recvs
    int nsend = 0;
    for (int i = 0; i < nproc; i++)
      if (ptags[i] != mytag)
        nsend++;

    std::vector<int> sndMap(nsend);
    for (int p = 0, m = 0; p < nproc; p++)
    {
      if (ptags[p] != mytag)
      {
        sndMap[m] = p;
        m++;
      }
    }

    pc->setMap(nsend, nsend, sndMap.data(), sndMap.data());

    mb->setupADT();

    mb->extraConn();
  }
}

void tioga::performConnectivity(void)
{
  doHoleCutting();

  doPointConnectivity();
}

void tioga::setIterIblanks(double dt, int nvar)
{
  // Determine blanking status for approximate grid at end of time step
  mb->calcNextGrid(dt);

  mb->preprocess();

  doHoleCutting();

  // Determine blanking status for current grid
  mb->resetCurrentGrid();

  mb->preprocess();

  doHoleCutting();

  // Determine final blanking status to use over time step
  int nunblank = mb->getIterIblanks();
if (nunblank > 0) printf("Rank %d unblanking %d elements\n",myid,nunblank);
  mb->calcFaceIblanks(meshcomm);

  MPI_Allreduce(MPI_IN_PLACE, &nunblank, 1, MPI_INT, MPI_SUM, scomm);

  if (nunblank > 0)
  {
    doPointConnectivity(); /// TODO: just do unblank cells only, no faces

    dataUpdate_artBnd(nvar, NULL, 0);

    mb->clearUnblanks();
  }
}

void tioga::unblankPart1(void)
{
  // Switch iblank_cell to separate memory range
  mb->swapPointers();

  mb->updateOBB();

  doHoleCutting();
}

void tioga::unblankPart2(int nvar)
{
  // Swap iblank_cell pointer back
  mb->resetCurrentGrid();

  mb->updateOBB();

  doHoleCutting();

  // Determine final blanking status to use over time step
  int nunblank = mb->getIterIblanks();

  if (nunblank > 0) printf("%d: Unblanking %d elements\n",myid,nunblank); /// DEBUGGING
  mb->calcFaceIblanks(meshcomm);

  MPI_Allreduce(MPI_IN_PLACE, &nunblank, 1, MPI_INT, MPI_SUM, scomm);

  if (nunblank > 0)
  {
    doPointConnectivity(true); /// TODO: just do unblank cells only, no faces

    dataUpdate_artBnd(nvar, NULL, 0);

    mb->clearUnblanks();
  }
}

#ifdef _GPU
#define TG_DIRECTCUT
#else
#define TG_NORMAL
#endif
void tioga::doHoleCutting(void)
{
#ifdef TG_NORMAL
  Timer tgTime("Normal Version: ");
  tgTime.startTimer();
  PUSH_NVTX_RANGE("TIOGA",2);
  // Generate structured map of solid boundary (hole) locations
  getHoleMap();

  // Send/Recv oriented bounding boxes to/from all ranks and setup sndMap / rcvMap
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
    mb->getCellIblanks(meshcomm);

    if (iartbnd)  // Only done by ranks with high-order Artificial Boundary grids
    {
      // Find artificial boundary faces
      mb->calcFaceIblanks(meshcomm);
    }
  }
  tgTime.stopTimer();
  POP_NVTX_RANGE;
#endif

#ifdef TG_DIRECTCUT
  PUSH_NVTX_RANGE("DirectCut", 2);
  Timer dcTime("Direct Cut: ");
  dcTime.startTimer();
  // Generate structured map of solid boundary (hole) locations
  if (holeMap == NULL || overMap.size() == 0)
  {
    getHoleMap();
    getOversetMap();
  }
  exchangeBoxes();
  directCut();
  mb->calcFaceIblanks(meshcomm);
  dcTime.stopTimer();
  POP_NVTX_RANGE;
#endif
}

void tioga::doPointConnectivity(bool unblanking)
{
  if (!ihighGlobal) return;

  // Get all fringe point locations (high-order face/cell points + low-order vertices)
  mb->getFringeNodes(unblanking);

  // Exchange new list of points, including high-order Artificial Boundary
  // face points or internal points (or fringe nodes for non-high order)
  exchangePointSearchData();

  // Search for donor cells for all given points
  mb->search();

  // Setup interpolation weights and such for final interp-point list
#ifdef _GPU
  mb->processPointDonorsGPU();
#else
  mb->processPointDonors();
#endif

#ifdef _GPU
  setupCommBuffersGPU();
#endif
}

#ifdef _GPU
void tioga::setupCommBuffersGPU(void)
{
  // Setup additional index arrars for GPU to automatically pack MPI buffers
  int nsend, nrecv, *sndMap, *rcvMap;
  pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);

  sndVPack.resize(nsend);
//  sndPack2.resize(nsend);
  rcvVPack.resize(nrecv);

//  pc->initPacketsV2(sndPack2,rcvVPack);
  pc->initPacketsV(sndVPack,rcvVPack);

//  mb->setupBuffersGPU(nsend, intData, sndPack2);
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

void tioga::directCut(void)
{
  int nDims = mb->nDims;

  /// TODO: Callbacks
  nCutFringe = mb->nCutFringe;
  nCutHole = mb->nCutHole;

  // Get cutting-group bounding-box data for this rank
  std::vector<double> faceNodesW; // wall face nodes
  std::vector<double> faceNodesO; // overset face nodes

  std::vector<double> bboxW, bboxO;
  int nvertf = mb->getCuttingFaces(faceNodesW, faceNodesO, bboxW, bboxO);

  // # nodes-per-face, and # of cutting faces, for each rank
  std::vector<int> nvertf_p(nproc);
  std::vector<int> nHoleFace_p(nproc), nOverFace_p(nproc);
  MPI_Allgather(&nvertf, 1, MPI_INT, nvertf_p.data(), 1, MPI_INT, scomm);
  MPI_Allgather(&nCutHole, 1, MPI_INT, nHoleFace_p.data(), 1, MPI_INT, scomm);
  MPI_Allgather(&nCutFringe, 1, MPI_INT, nOverFace_p.data(), 1, MPI_INT, scomm);

  // Setup buffers for receiving cutting surfaces from each grid

  std::vector<MPI_Request> sreqs, rreqs;
  sreqs.reserve(2*nproc); rreqs.reserve(2*nproc);
  std::vector<std::vector<double>> faceNodes_g(nGrids);
  std::vector<std::vector<double>> bbox_g(nGrids);
  std::vector<int> nFace_p(nproc), nFaceTot_g(nGrids);
  std::vector<double> bbox_tmp(2*nDims*nproc);
  std::vector<double> aabb_tmp(2*nDims*nproc);

  for (int g = 0; g < nGrids; g++)
  {
    if (g == mytag) continue;

    bbox_g[g].resize(2 * nDims);

    for (int d = 0; d < nDims; d++)
    {
      bbox_g[g][d]       =  BIG_DOUBLE;
      bbox_g[g][d+nDims] = -BIG_DOUBLE;
    }
  }

  // Send / recv cutting group bounding-box data

  for (int p = 0; p < nproc; p++)
  {
    int g = gridIDs[p];
    if (g == mytag) continue;

    nFace_p[p] = (gridType == 0) ? nOverFace_p[p] : nHoleFace_p[p];
    nFaceTot_g[g] += nFace_p[p];

    if (nFace_p[p] > 0)
    {
      // Receive the rank's cutting-face and overall bounding box
      rreqs.emplace_back();
      MPI_Irecv(&bbox_tmp[2*nDims*p], 2*nDims, MPI_DOUBLE, p, 1, scomm, &rreqs.back());
      sreqs.emplace_back();
      MPI_Isend(mb->aabb, 2*nDims, MPI_DOUBLE, p, 2, scomm, &sreqs.back());
    }
    
    // Get pointers to face node & bounding-box data based on grid type
    double *Bptr;
    int size = 0;
    if (nCutHole > 0 && gridTypes[p] > 0)
    {
      Bptr = bboxW.data();
      size = nCutHole;
    }
    else if (nCutFringe > 0 && gridTypes[p] == 0)
    {
      Bptr = bboxO.data();
      size = nCutFringe;
    }

    if (size > 0)
    {
      // Send bounding box
      sreqs.emplace_back();
      MPI_Isend(Bptr, 2*nDims, MPI_DOUBLE, p, 1, scomm, &sreqs.back());
      rreqs.emplace_back();
      MPI_Irecv(&aabb_tmp[2*nDims*p], 2*nDims, MPI_DOUBLE, p, 2, scomm, &rreqs.back());
    }
  }

  // Complete the communication
  MPI_Waitall((int)sreqs.size(), sreqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall((int)rreqs.size(), rreqs.data(), MPI_STATUSES_IGNORE);

  // Figure out which ranks cut faces must be sent/received from

  std::vector<int> dcSndMap;  // List of ranks to send cut faces to
  std::vector<int> dcRcvMap;  // List of ranks to recv cut faces from

  dcSndMap.resize(0);
  dcRcvMap.resize(0);

  for (int p = 0; p < nproc; p++)
  {
    int g = gridIDs[p];
    if (g == mytag) continue;

    // Check if our bbox overlaps with other rank's cut-group bbox
    if (nFace_p[p] > 0 && tg_funcs::boundingBoxCheck(mb->aabb, &bbox_tmp[2*nDims*p], 3))
      dcRcvMap.push_back(p);

    double *Bptr;
    int size = 0;
    if (nCutHole > 0 && gridTypes[p] > 0)
    {
      Bptr = bboxW.data();
      size = nCutHole;
    }
    else if (nCutFringe > 0 && gridTypes[p] == 0)
    {
      Bptr = bboxO.data();
      size = nCutFringe;
    }

    // Check if other rank's bbox overlaps with our cut-group bbox
    if (size > 0 && tg_funcs::boundingBoxCheck(Bptr, &aabb_tmp[2*nDims*p], 3))
      dcSndMap.push_back(p);
  }

  int nrecv = dcRcvMap.size();
  int nsend = dcSndMap.size();

  // Figure out displacements into per-grid storage (rather than per-rank)

  std::vector<int> nFace_g(nGrids), nFace_r(nrecv);
  std::vector<int> faceDisp_r(nrecv);
  std::vector<int> nProc_g(nGrids), nVertf_g(nGrids);

  for (int i = 0; i < nrecv; i++)
  {
    int p = dcRcvMap[i];
    int g = gridIDs[p];

    faceDisp_r[i] = nFace_g[g];
    if (gridType > 0)
      nFace_r[i] += nHoleFace_p[p];
    else
      nFace_r[i] += nOverFace_p[p];

    nFace_g[g] += nFace_r[i];
    nVertf_g[g] = nvertf_p[p];
    nProc_g[g]++;
  }

  for (int g = 0; g < nGrids; g++)
  {
    faceNodes_g[g].resize(nFace_g[g] * nVertf_g[g] * nDims);
  }

  // Reset the request buffers
  sreqs.resize(0); sreqs.reserve(nsend);
  rreqs.resize(0); rreqs.reserve(nrecv);

  // Send / recv face nodes data - only to required ranks

  for (int i = 0; i < nrecv; i++)
  {
    int p = dcRcvMap[i];
    int g = gridIDs[p];

    if (nFace_r[i] > 0)
    {
      // Receive the cutting faces
      rreqs.emplace_back();
      MPI_Irecv(faceNodes_g[g].data() + faceDisp_r[i]*nVertf_g[g]*nDims,
          nFace_r[i]*nvertf_p[p]*nDims, MPI_DOUBLE, p, 0, scomm, &rreqs.back());
    }
  }

  for (int i = 0; i < nsend; i++)
  {
    int p = dcSndMap[i];

    // Get pointers to face node & bounding-box data based on grid type
    double *Fptr;
    int size = 0;
    if (nCutHole > 0 && gridTypes[p] > 0)
    {
      Fptr = faceNodesW.data();
      size = nCutHole*nvertf*nDims;
    }
    else if (nCutFringe > 0 && gridTypes[p] == 0)
    {
      Fptr = faceNodesO.data();
      size = nCutFringe*nvertf*nDims;
    }

    if (size > 0)
    {
      // Send face nodes
      sreqs.emplace_back();
      MPI_Isend(Fptr, size, MPI_DOUBLE, p, 0, scomm, &sreqs.back());
    }
  }

  // Complete the communication
  MPI_Waitall((int)sreqs.size(), sreqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall((int)rreqs.size(), rreqs.data(), MPI_STATUSES_IGNORE);

  // Reduce the bounding box data to one box per grid
  //for (int p = 0; p < nproc; p++)
  for (int i = 0; i < nrecv; i++)
  {
    int p = dcRcvMap[i];
    int g = gridIDs[p];

    for (int d = 0; d < nDims; d++)
    {
      bbox_g[g][d] = std::min(bbox_g[g][d], bbox_tmp[2*nDims*p+d]);
      bbox_g[g][d+nDims] = std::max(bbox_g[g][d+nDims], bbox_tmp[2*nDims*p+d+nDims]);
    }
  }

  // Do the cutting
  std::vector<CutMap> cutMap(nGrids);

  Timer cutTime("Cutting Time: ");

  int ncut = 0;
  for (int g = 0; g < nGrids; g++)
  {
    if (g == mytag) continue;

    cutTime.startTimer();
    if (nFaceTot_g[g] > 0)
    {
      HOLEMAP hm = (gridType == 0) ? overMap[g] : holeMap[g];
#ifdef _GPU
      mb->directCut_gpu(faceNodes_g[g], nFace_g[g], nVertf_g[g], bbox_g[g], hm, cutMap[ncut], gridType);
#else
      mb->directCut(faceNodes_g[g], nFace_p[g], nVertf_g[g], bbox_g[g], cutMap[ncut], gridType);
#endif
      ncut++;
    }
    cutTime.stopTimer();
  }

//  cutTime.showTime(3);

  cutMap.resize(ncut);
  mb->unifyCutFlags(cutMap);

//  if (cutMap.size() > 0)
//    mb->writeCellFile(myid, cutMap.back().flag.data()); /// DEBUGGING
//  else
//    mb->writeCellFile(myid, NULL); /// DEBUGGING
//  MPI_Barrier(scomm);
//  exit(0); /// DEBUGGING
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
//  mb->writeCellFile(myid);
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
  dataUpdate_artBnd_send(nvar,dataFlag);

  dataUpdate_artBnd_recv(nvar,dataFlag);
}

void tioga::dataUpdate_artBnd_send(int nvar, int dataFlag)
{
  // initialize send and recv packets
  int nsend,nrecv;
  int *sndMap,*rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  if (nsend+nrecv == 0) return;

  // get the interpolated solution now
  int stride = nvar;
  if (iartbnd && dataFlag == 1) stride = 3*nvar;

  ubuf_d.resize(mb->ninterp2*stride);
  ubuf_h.resize(mb->ninterp2*stride);

  for (int k = 0; k < nsend; k++)
  {
    sndVPack[k].nreals = sndVPack[k].nints * stride;
    sndVPack[k].realData.resize(sndVPack[k].nreals);
//    sndPack2[k].nreals = sndPack2[k].nints * stride; /// newer method
  }

  interpTime.startTimer();

  if (dataFlag == 0)
    mb->interpSolution_gpu(ubuf_d.data(), nvar);
  else
    mb->interpGradient_gpu(ubuf_d.data(), nvar);

  ubuf_h.assign(ubuf_d.data(), ubuf_d.size(), &mb->stream_handle);

  interpTime.stopTimer();

  /* ------------------------------------------------------------------------ */
  /* Version 1: Wait for D2H transfer to complete and pack separate buffer */
  // Wait for device-to-host transfer to complete
  cudaStreamSynchronize(mb->stream_handle);
  PUSH_NVTX_RANGE("tg_packBuffers", 3);
  // Populate the packets [organize interp data by rank to send to]
  for (int p = 0; p < nsend; p++)
  {
    for (int i = 0; i < sndVPack[p].nints; i++)
      for (int j = 0; j < stride; j++)
        sndVPack[p].realData[i*stride+j] = ubuf_h[(mb->buf_disp[p]+i)*stride+j];
  }
  POP_NVTX_RANGE;
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Version 2: Let D2H be async and get pointers into buffer for later send
  PUSH_NVTX_RANGE("tg_packBuffers", 3);
  // Populate the packets [organize interp data by rank to send to]
  for (int p = 0; p < nsend; p++)
    sndPack2[p].realData = ubuf_h.data() + mb->buf_disp[p]*stride;
  POP_NVTX_RANGE;
  /* ------------------------------------------------------------------------ */
}

void tioga::dataUpdate_artBnd_recv(int nvar, int dataFlag)
{
  // initialize send and recv packets
  int nsend,nrecv;
  int *sndMap,*rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  if (nsend+nrecv == 0) return;

  // get the interpolated solution now
  int stride = nvar;
  if (iartbnd && dataFlag == 1) stride = 3*nvar;

  // communicate the data across all partitions
  PUSH_NVTX_RANGE("tg_pc_send", 0);

  /* ------------------------------------------------------------------------ */
  /* Version 1.2: Do the waiting and buffer packing here instead
  // Wait for device-to-host transfer to complete
  cudaStreamSynchronize(mb->stream_handle);
  PUSH_NVTX_RANGE("tg_packBuffers", 3);
  // Populate the packets [organize interp data by rank to send to]
  for (int p = 0; p < nsend; p++)
  {
    for (int i = 0; i < sndVPack[p].nints; i++)
      for (int j = 0; j < stride; j++)
        sndVPack[p].realData[i*stride+j] = ubuf_h[(mb->buf_disp[p]+i)*stride+j];
  }
  POP_NVTX_RANGE;
  /* ------------------------------------------------------------------------ */

  pc->sendPacketsV(sndVPack,rcvVPack);
//  pc->sendPacketsV2(sndPack2,rcvVPack); /// newer method
  POP_NVTX_RANGE;

  // Wait on all of the sends/recvs for the interpolated data
  PUSH_NVTX_RANGE("tg_pc_recv", 0);
  pc->recvPacketsV();
  POP_NVTX_RANGE;

  // Decode the packets and update the values in the solver's data array
  PUSH_NVTX_RANGE("tg_unpack_data", 1);
  if (ihigh)
  {
    fringebuf_h.resize(mb->ntotalPoints*stride);
    recv_itmp.assign(mb->ntotalPoints,0);

    for (int k = 0; k < nrecv; k++)
    {
      for (int i = 0; i < rcvVPack[k].nints; i++)
      {
        int ind = rcvVPack[k].intData[i];
        recv_itmp[ind] = 1;
        for (int j = 0; j < stride; j++)
        {
          fringebuf_h[ind*stride+j] = rcvVPack[k].realData[i*stride+j];
        }
      }
    }

    // Print any 'orphan' points which may exist
    std::ofstream fp;
    int norphanPoint = 0;
    for (int i = 0; i < mb->ntotalPoints; i++)
    {
      if (recv_itmp[i] == 0) {
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
    if (norphanPoint > 0) ThrowException("Orphan points found!");
    if (norphanPoint > 0 && iorphanPrint) {
      printf("Warning::number of orphans in rank %d = %d of %d\n",myid,norphanPoint,mb->ntotalPoints);
      iorphanPrint = 0;
    }

    // change the state of cells/nodes who are orphans
    if (!iartbnd)
      mb->clearOrphans(recv_itmp.data());

    if (iartbnd)
    {
      interpTime.startTimer();
      if (dataFlag == 0)
        mb->updateFringePointData(fringebuf_h.data(),nvar);
      else
        mb->updateFringePointGradient(fringebuf_h.data(),nvar);
      interpTime.stopTimer();
    }
    else
      ThrowException("Not written for non-artificial boundary codes right now");
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
  dataUpdate_artBnd_send(nvar,q_spts,dataFlag);

  dataUpdate_artBnd_recv(nvar,q_spts,dataFlag);
}

void tioga::dataUpdate_artBnd_send(int nvar, double *q_spts, int dataFlag)
{
  // initialize send and recv packets
  int nsend,nrecv;
  int *sndMap,*rcvMap;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  sndVPack.resize(nsend);
  rcvVPack.resize(nrecv);

  if (nsend+nrecv == 0) return;

  std::vector<int> icount(nsend), dcount(nsend);

  pc->initPacketsV(sndVPack,rcvVPack);

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
    sndVPack[k].nints++;
    sndVPack[k].nreals += stride;
  }

  for (int k = 0; k < nsend; k++)
  {
    sndVPack[k].intData.resize(sndVPack[k].nints);
    sndVPack[k].realData.resize(sndVPack[k].nreals);
    icount[k] = dcount[k] = 0;
  }

  for (int i = 0; i < nints; i++)
  {
    int k = integerRecords[2*i];
    sndVPack[k].intData[icount[k]++] = integerRecords[2*i+1];
    for (int j = 0; j < stride; j++)
      sndVPack[k].realData[dcount[k]++] = realRecords[stride*i+j];
  }

  // communicate the data across all partitions
  pc->sendPacketsV(sndVPack,rcvVPack);

  free(integerRecords);
  free(realRecords);
}

void tioga::dataUpdate_artBnd_recv(int nvar, double* q_spts, int dataFlag)
{
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  pc->recvPacketsV();

  int stride = nvar;
  if (iartbnd && dataFlag == 1) stride = 3*nvar;

  // Decode the packets and update the values in the solver's data array
  if (ihigh)
  {
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

    if (iartbnd)
    {
      interpTime.startTimer();
      if (dataFlag == 0)
        mb->updateFringePointData(qtmp.data(),nvar);
      else
        mb->updateFringePointGradient(qtmp.data(),nvar);
      interpTime.startTimer();
    }
    else
      mb->updatePointData(q_spts,qtmp.data(),nvar,dataFlag);
  }
  else
  {
    for (int k = 0; k < nrecv; k++)
    {
      for (int i = 0; i < rcvVPack[k].nints; i++)
      {
        mb->updateSolnData(rcvVPack[k].intData[i],&rcvVPack[k].realData[nvar*i],q_spts,nvar,dataFlag);
      }
    }
  }
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

#ifdef _GPU
void tioga::set_stream_handle(cudaStream_t handle, cudaEvent_t event)
{
  mb->set_stream_handle(handle, event);
}
#endif
