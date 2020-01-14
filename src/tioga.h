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

#ifndef TIOGA_H
#define TIOGA_H

#include <vector>
#include <map>
#include <memory>
#include <stdint.h>
#include "MeshBlock.h"
#include "CartGrid.h"
#include "CartBlock.h"
#include "parallelComm.h"

/** Define a macro entry flagging the versions that are safe to use with large
 *  meshes containing element and node IDs greater than what a 4-byte signed int
 *  can support
 */
#define TIOGA_HAS_UINT64T 1

/**
 * Topology Indpendent Overset Grid Assembler (TIOGA)
 * Base class and dependencies
 * The methods of this class are invoked from tiogaInterface.C
 *
 *  Jay Sitaraman 02/24/2014
 */

namespace TIOGA {

class tioga
{
 private :
  int nblocks;
  int nblocksq;
  int ncart;
  MeshBlock *mb;
  CartGrid *cg;
  CartBlock *cb;
  int nmesh;
  HOLEMAP *holeMap;
  MPI_Comm scomm;
  parallelComm *pc;
  parallelComm *pc_cart;
  int isym;
  int ierr;
  int myid,numprocs;
  int *sendCount;
  int *recvCount;
  //OBB *obblist;
  std::vector<OBB> obblist;
  int iorphanPrint;

  //! Mesh blocks in this processor 
  std::vector<std::unique_ptr<MeshBlock> > mblocks;
  //! Solver assigned mesh tags for the mesh blocks
  std::vector<int> mtags;
  std::vector<int> mytag;
  //! Mesh tag to local block index lookup mapping
  std::map<int, int> tag_iblk_map;

  //! Intersect block unique ID to index lookup mapping
  std::map<int, int> intBoxMap;

  //! Parallel comms to obblist indicies
  std::vector<int> ibsPerProc;
  std::vector<std::vector<int>> ibProcMap;
  //! q-variables registered
  double **qblock;


 public:
  int ihigh;
  int ihighGlobal;
  int iamrGlobal;
  int mexclude,nfringe;
  /** basic constuctor */
  tioga()
    /*
    : mblocks(0),
      mtags(0)
    */
    {
        mb = NULL; cg=NULL; cb=NULL;
        holeMap=NULL; pc=NULL; sendCount=NULL; recvCount=NULL;
        pc_cart = NULL;
        // obblist=NULL; isym=2;ihigh=0;nblocks=0;ncart=0;ihighGlobal=0;iamrGlobal=0;
        isym=3;ihigh=0;nblocks=0;ncart=0;ihighGlobal=0;iamrGlobal=0;
        mexclude=3,nfringe=1;
        qblock=NULL;
        mblocks.clear();
        mtags.clear();
    }
 
  /** basic destructor */
  ~tioga(); 
  
  /** set communicator */
  void setCommunicator(MPI_Comm communicator,int id_proc,int nprocs);

  /** registerGrid data */

  void registerGridData(int btag,int nnodes,double *xyz,int *ibl, int nwbc,int nobc,
                        int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn,
                        uint64_t* cell_gid=NULL, uint64_t* node_gid=NULL);

  void registerSolution(int btag,double *q);

  void profile(void);

  void exchangeBoxes(void);

  void exchangeSearchData(int at_points=0);

  void exchangeDonors(void);
    
  /** perform overset grid connectivity */

  void performConnectivity(void);
  void performConnectivityHighOrder(void);
  void performConnectivityAMR(void);

  /** update data */

  void dataUpdate(int nvar,int interptype,int at_points=0) ;

  void dataUpdate_AMR(int nvar,int interptype) ;
  
  void dataUpdate_highorder(int nvar,double *q,int interptype) ;

  /** get hole map for each mesh */
 
  void getHoleMap(void);

  /** output HoleMaps */
  
  void outputHoleMap(void);

  void writeData(int nvar,int interptype);

  void getDonorCount(int btag, int *dcount, int *fcount);
  
  void getDonorInfo(int btag, int *receptors,int *indices,double *frac,int *dcount);

  void getReceptorInfo(std::vector<int>&);

  /** set symmetry bc */
  void setSymmetry(int syminput) { isym=syminput;};
  /** set resolutions for nodes and cells */
  void setResolutions(double *nres,double *cres)
  { auto & mb = mblocks[0]; mb->setResolutions(nres,cres);}

  void setResolutions(int btag, double *nres,double *cres)
  {
    auto idxit = tag_iblk_map.find(btag);
    int iblk = idxit->second;
    auto& mb = mblocks[iblk];
    mb->setResolutions(nres, cres);
  }
  
  void setMexclude(int *mexclude_input)
  {
    mexclude=*mexclude_input;
  }

  void setNfringe(int *nfringe_input)
  {
    nfringe=*nfringe_input;
  }

  void set_cell_iblank(int *iblank_cell)
  {
   auto& mb = mblocks[0];
   mb->set_cell_iblank(iblank_cell);
  }

  void set_uniform_hex_flag(int btag, int flag)
  {
      auto idxit = tag_iblk_map.find(btag);
      int iblk = idxit->second;
      auto& mb = mblocks[iblk];
      mb->check_uniform_hex_flag = flag;
  }

  void set_cell_iblank(int btag, int* ib_cell)
  {
    auto idxit = tag_iblk_map.find(btag);
    int iblk = idxit->second;
    auto& mb = mblocks[iblk];
    mb->set_cell_iblank(ib_cell);
  }

  void setcallback(void (*f1)(int*, int*),
		    void (*f2)(int *,int *,double *),
		    void (*f3)(int *,double *,int *,double *),
		    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
		   void (*f5)(int *,int *,double *,int *,int*,double *))
  {
   for(int ib=0;ib<nblocks;ib++)
   {
    auto& mb = mblocks[ib];
    mb->setcallback(f1,f2,f3,f4,f5);
   }   
   ihigh=1;
  }

  void setp4estcallback(void (*f1)(double *,int *,int *,int *),
			void (*f2) (int *,int *))
  {
   for(int ib=0;ib<nblocks;ib++)
    { 
     auto& mb = mblocks[ib];   // TODO:this may have to based on unique tag of p4est blocks
     mb->setp4estcallback(f1,f2);
    } 
  }

  void set_p4est(void)
  {
    for(int ib=0;ib < nblocks;ib++)
    {
      mytag[ib]=-mytag[ib];
      auto& mb = mblocks[ib]; // TODO
      mb->resolutionScale=1000.0;
    }
  }
  
  void set_amr_callback(void (*f1)(int *,double *,int *,double *))
  {
    cg->setcallback(f1);
  }

  void register_amr_global_data(int, int, double *, int *,double *, int, int);
  void set_amr_patch_count(int);
  void register_amr_local_data(int, int ,int *, double *);  
  void exchangeAMRDonors(void);
  void checkComm(void);
  void outputStatistics(void);
  void myTimer(char const *, int);
  void reduce_fringes(void);
};
      
  
}

#endif /* TIOGA_H */
