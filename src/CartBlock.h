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

#ifndef CARTBLOCK_H
#define CARTBLOCK_H

#include <cassert>
#include <cstdlib>
#include <vector>
#include "codetypes.h"

struct INTERPLIST2;
struct DONORLIST;
struct HOLEMAP;

namespace TIOGA {
struct AMRMeshInfo;
}

class CartGrid;
class CartBlock
{
 private:
  int local_id;
  int global_id;
  int dims[3],nf,ncell,ncell_nf,nnode,nnode_nf;
  int nvar_cell,nvar_node;
  int myid;
  int *ibl_cell, *ibl_node;
  double *qcell, *qnode;
  double xlo[3]; 
  double dx[3];
  int ndonors;
  int interpListSize;
  int ninterp,nweights;
  INTERPLIST2 *interpList,*listptr;
  DONORLIST **donorList;
  std::vector<int> q_fringe_ind_full, q_fringe_ind_cell_nd;
  std::vector<double> q_fringe;

  // these are device pointers
  double *qcell_d,*qnode_d;
  int block_size{128};
  int *interpList_wcft{nullptr};
  double *interpList_weights{nullptr};
  int *interpList_inode{nullptr};

  void (*donor_frac) (int *,double *,int *,double *);
 public:
  CartBlock() { global_id=0;dims[0]=dims[1]=dims[2]=0;ibl_cell=NULL;ibl_node=NULL;qcell=NULL;qnode=NULL;interpListSize=0;donorList=NULL;interpList=NULL;
    donor_frac=nullptr;nvar_cell=0;nvar_node=0;};
  ~CartBlock() { clearLists();};

  void registerData(int lid, TIOGA::AMRMeshInfo* minfo);
  void registerSolution(int lid, TIOGA::AMRMeshInfo* minfo);

  void registerData(int local_id_in,int global_id_in,int *iblankin,int *iblanknin)
  {
    local_id=local_id_in;
    global_id=global_id_in;
    ibl_cell=iblankin;
    ibl_node=iblanknin;
  };
  void registerSolution(double *qin,int nq_cell,int nq_node) {
    assert ((nq_cell+nq_node > 0)
      && !((nq_cell > 0) && (nq_node > 0)));
    if (nq_node > 0) {
      nvar_node = nq_node;
      qnode = qin;
    }
    else if (nq_cell > 0){
      nvar_cell = nq_cell;
      qcell = qin;
    }
  };
  int num_cell_var() const { return nvar_cell; }
  int num_node_var() const { return nvar_node; }
  void preprocess(CartGrid *cg);
  void getInterpolatedData(int *nints,int *nreals,int **intData,
			   double **realData);
  void  pushInterpListsToDevice(void);
  void getInterpolatedDataDevice(double *realData,int nvar_cell,int nvar_node);
  void assembleFringeSolution(double *qval,int index);
  void update(double *qval,int index);
  void updateDevice();
  void getCancellationData(int *cancelledData, int *ncancel);
  void processDonors(HOLEMAP *holemap, int nmesh);
  void processIblank(HOLEMAP *holemap, int nmesh, bool isNodal);
  void insertInDonorList(int senderid,int index,int meshtagdonor,int remoteid,int remoteblockid,double cellRes);
  void insertInInterpList(int procid,int remoteid,int remoteblockid,double *xtmp);
  void writeCellFile(int bid);
  void clearLists(void);
  void initializeLists(void);
};

#endif /* CARTBLOCK_H */
