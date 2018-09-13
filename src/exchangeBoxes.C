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

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cassert>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
extern "C"{
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
			double vB[3][3],double xB[3],double dxB[3]);			   
}
void tioga::exchangeBoxes(void)
{
  int *sndMap;
  int *rcvMap;
  int nsend;
  int nrecv;
  PACKET *sndPack,*rcvPack;

  std::vector<int> nbPerProc(numprocs); // Number of chunks per processor
  MPI_Allgather(&nblocks, 1, MPI_INT, nbPerProc.data(), 1, MPI_INT, scomm);

  // Total number mesh chunks across all procs
  int ntotalblks = std::accumulate(nbPerProc.begin(),nbPerProc.end(),0);

  std::vector<int> alltags(ntotalblks); // Mesh tags for all blocks across all procs
  std::vector<int> displs(numprocs+1);  // Offsets for tags per proc
  std::vector<bool> sendFlag(numprocs,false); // Flag indicating send/recv from this proc

  displs[0] = 0;
  for (int i=1; i <= numprocs; i++) {
    displs[i] = displs[i-1] + nbPerProc[i-1];
  }

  MPI_Allgatherv(mtags.data(), nblocks, MPI_INT, alltags.data(),
                 nbPerProc.data(), displs.data(), MPI_INT, scomm);

  int maxtag = -1;
  //for (auto itag: alltags)
  for(int i=0;i<ntotalblks;i++) {
    int itag=alltags[i];
    if (maxtag < itag) maxtag = itag;
  }
  int mxtgsqr = maxtag * maxtag;

  //
  // count number of other processors to communicate to
  // in overset grid scenario, usually you do not communicate
  // to mesh blocks that carry the same mesh tag (i.e. you don't
  // talk to your sister partitions)
  //
  nsend = 0;
  for (int p=0; p < numprocs; p++) {

    for(int d=displs[p]; d < displs[p+1]; d++) {
      for (int ib=0; ib < nblocks; ib++) {
        if (alltags[d] != mtags[ib]) {
          nsend++;
          sendFlag[p] = true;
          break;
        }
      }
      if (sendFlag[p]) break;
    }
  }

  // In general we communicate forward
  // and backward, separate lists are maintained for
  // flexibility
  //
  nrecv = nsend;

  sndMap = (int*)malloc(sizeof(int) * nsend);
  rcvMap = (int*)malloc(sizeof(int) * nrecv);

  for (int i=0, ig=0; i < numprocs; i++) {
    if (sendFlag[i] == true) {
      sndMap[ig] = i;
      rcvMap[ig] = i;
      ig++;
    }
  }
  pc->setMap(nsend,nrecv,sndMap,rcvMap);

  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);

  // Populate sndPack data
  for (int k=0; k < nsend; k++) {
    std::vector<bool> bFlag(nblocks, false);
    int ip = sndMap[k];
    sndPack[k].nints = 0;

    // Determine buffer sizes
    for (int ib=0; ib < nblocks; ib++) {
      for (int d=displs[ip]; d < displs[ip+1]; d++) {
        if (( alltags[d] != mtags[ib] ) &&
            ( bFlag[ib] == false)) {
          sndPack[k].nints++;
          bFlag[ib] = true;
        }
      }
    }
    sndPack[k].nreals = 15 * sndPack[k].nints;
    sndPack[k].intData = (int*) malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double*) malloc(sizeof(double)*sndPack[k].nreals);

    // Populate int and real data arrays
    int im = 0;
    int m = 0;
    for (int ib=0; ib < nblocks; ib++) {
      if (bFlag[ib] == false) continue;

      sndPack[k].intData[im++] = mtags[ib];

      for (int i=0; i < 3; i++)
        for (int j=0; j< 3; j++)
          sndPack[k].realData[m++] = mblocks[ib]->obb->vec[i][j];

      for (int i=0; i< 3; i++)
        sndPack[k].realData[m++] = mblocks[ib]->obb->xc[i];

      for (int i=0; i< 3; i++)
        sndPack[k].realData[m++] = mblocks[ib]->obb->dxc[i];
    }
  }
  pc->sendRecvPackets(sndPack, rcvPack);

  // Determine total number of OBBs received 
  int nobb = 0;
  for (int k=0; k < nrecv; k++)
    nobb += rcvPack[k].nints;

  // Store all received OBBs in a temporary list
  std::vector<OBB> obbRecv(nobb);
  std::vector<int> obbID(nobb); // Mesh tag corresponding to the OBB
  std::vector<int> obbProc(nobb); // Proc ID corresponding to OBB

  for(int k=0, ix=0; k < nrecv; k++) {
    int m = 0;
    int im = 0;
    for (int n=0; n < rcvPack[k].nints; n++) {
      obbProc[ix] = rcvMap[k];
      obbID[ix] = rcvPack[k].intData[im++];
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          obbRecv[ix].vec[i][j] = rcvPack[k].realData[m++];

      for (int i=0; i<3; i++)
        obbRecv[ix].xc[i] = rcvPack[k].realData[m++];

      for (int i=0; i<3; i++)
        obbRecv[ix].dxc[i] = rcvPack[k].realData[m++];

      ix++;
    }
  }

  // Mapping of (local block_id, remote OBB block_id) for every intersected pair
  std::vector<std::pair<int, int>> intersectIDs;
  // Counter tracking number of intersected blocks per partition
  std::vector<int> obPerProc(numprocs,0);

  // Reset sendflags 
  std::fill(sendFlag.begin(), sendFlag.end(), false);
  //
  // Check for intersection of OBBs
  //
  for (int ob=0; ob < nobb; ob++) {
    for (int ib=0; ib < nblocks; ib++) {
      auto& mb = mblocks[ib];
      int meshtag = mb->getMeshTag();

        if (obbID[ob] == meshtag) continue;

        if ( obbIntersectCheck(
               mb->obb->vec, mb->obb->xc, mb->obb->dxc,
               obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc) ||
             obbIntersectCheck(
               obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc,
               mb->obb->vec, mb->obb->xc, mb->obb->dxc)) {
          // If there is an intersection, store the index pair, increment number
          // of intersections for this processor, and activate send flag
          intersectIDs.push_back(std::make_pair(ib, ob));
          obPerProc[obbProc[ob]]++;
          sendFlag[obbProc[ob]] = true;
      }
    }
  }
  int new_send = std::count(sendFlag.begin(), sendFlag.end(), true);
  assert(new_send <= nsend);

  // Populate send and recv maps
  std::map<int,int> invMap;
  nsend = nrecv = new_send;
  for (int p=0, ip=0; p < numprocs; p++)
    if (sendFlag[p]) {
      sndMap[ip] = p;
      rcvMap[ip] = p;
      invMap[p] = ip;
      ip++;
    }

  // clear packets before nsend and nrecv are modified in pc->setMap
  pc->clearPackets(sndPack,rcvPack);
  pc->setMap(nsend,nrecv,sndMap,rcvMap);


  // if (obblist) TIOGA_FREE(obblist);
  // obblist = (OBB*) malloc(sizeof(OBB) * intersectIDs.size());
  intBoxMap.clear();
  ibsPerProc.clear();
  ibProcMap.clear();
  obblist.clear();
  obblist.resize(intersectIDs.size());
  ibsPerProc.resize(nsend);
  ibProcMap.resize(nsend);

  // Determine packet sizes and reallocate arrays
  for(int k=0; k < nsend; k++) {
    sndPack[k].nints = 3 * obPerProc[sndMap[k]];
    sndPack[k].nreals = sndPack[k].nints * 6;
    sndPack[k].intData = (int*) malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double*) malloc(sizeof(double)*sndPack[k].nreals);
    ibsPerProc[k] = obPerProc[sndMap[k]];
    ibProcMap[k].resize(ibsPerProc[k]);
  }

  // Array tracking indices for populating reduced OBBs
  std::vector<int> idxOffset(nsend,0);
  for (size_t i=0; i<intersectIDs.size(); i++){
    auto ids = intersectIDs[i];
    int ib = ids.first;           // Block ID of the local mesh block
    int ob = ids.second;          // Index of the intersected block in OBB list 
    int k = invMap[obbProc[ob]];  // Index in sndMap for this proc ID
    auto& mb = mblocks[ib];       // Mesh block data object
    int ip = obbProc[ob];

    int ioff = idxOffset[k];      // Index to fill in sndPack
    int roff = idxOffset[k] * 6;

    int key_recv = mxtgsqr * ip + maxtag * (mtags[ib] - 1) + obbID[ob] - 1;
    int key_send = mxtgsqr * myid + maxtag * (obbID[ob]-1) + (mtags[ib]-1);
    intBoxMap[key_recv] = i;
    ibProcMap[k][ioff] = i;
    obblist[i].comm_idx = k;
    obblist[i].iblk_local = ib;
    // obblist[i].iblk_remote = obbID[ob];
    obblist[i].send_tag = key_send;
    obblist[i].recv_tag = key_recv;

    sndPack[k].intData[3*ioff] = key_send; // mb->getMeshTag();
    sndPack[k].intData[3*ioff+1] = ib;
    sndPack[k].intData[3*ioff+2] = mb->getMeshTag();
    mb->getReducedOBB(&obbRecv[ob], &(sndPack[k].realData[roff]));

    for(int ii=0; ii<3; ii++)
      for(int j=0; j<3; j++)
        obblist[i].vec[ii][j] = obbRecv[ob].vec[ii][j];

    // Increment index offset for next fill
    idxOffset[k]++;
    // std::cout << "# " << myid << " " << obbProc[ob] << " "
    //           << mtags[ib] << " " << obbID[ob] << " "
    //           << std::setw(3) << key_recv << " "
    //           << std::setw(3) << key_send << std::endl;
  }
  pc->sendRecvPackets(sndPack,rcvPack);

  for (int k=0; k<nrecv; k++) {
    int m=0;

    for (int n=0; n< rcvPack[k].nints; n+=3) {
      int key = rcvPack[k].intData[n];
      int ii = intBoxMap[key];
      obblist[ii].iblk_remote = rcvPack[k].intData[n+1];
      obblist[ii].tag_remote = rcvPack[k].intData[n+2];

      for (int i=0; i<3; i++)
        obblist[ii].xc[i] = rcvPack[k].realData[m++];

      for (int i=0; i<3; i++)
        obblist[ii].dxc[i] = rcvPack[k].realData[m++];

    }
  }

  pc->clearPackets(sndPack,rcvPack);
  //
  // Free local memory
  //
  TIOGA_FREE(sndMap);
  TIOGA_FREE(rcvMap);
  TIOGA_FREE(sndPack);
  TIOGA_FREE(rcvPack);
}
