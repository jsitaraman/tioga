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
/**
 * Topology Indpendent Overset Grid Assembler (TIOGA)
 * Base class and dependencies
 * The methods of this class are invoked from tiogaInterface.C
 *
 *  Jay Sitaraman 02/24/2014 
 */
#include "MeshBlock.h"
#include "CartGrid.h"
#include "CartBlock.h"
#include "parallelComm.h"

class tioga
{
 private :
  int nblocks;
  int ncart;
  MeshBlock *mb;
  CartGrid *cg;
  CartBlock *cb;
  int nmesh;
  HOLEMAP *holeMap;
  MPI_Comm scomm;
  MPI_Comm meshcomm;
  parallelComm *pc;
  parallelComm *pc_cart;
  int isym;
  int ierr;
  int mytag;
  int myid,numprocs;
  int *sendCount;
  int *recvCount;
  OBB *obblist;
  int iorphanPrint;

  int nprocMesh, meshRank;

 public:
  int ihigh;        /// High-Order flag for current rank
  int iartbnd;      /// Artificial-boundary flag for current rank
  int ihighGlobal;  /// Flag for whether high-order grids exist on any rank
  int iamrGlobal;   /// Flag for whether AMR cartesian grids exist on any rank
  int iabGlobal;    /// Flaag for whether high-order A.B.'s being used on any rank

  /** basic constuctor */
  tioga()
  {
    mb = NULL; cg = NULL; cb = NULL;
    pc = NULL; sendCount = NULL; recvCount = NULL;
    holeMap = NULL; obblist = NULL;
    nblocks = 0; ncart = 0;
    isym = 2; ihigh = 0; iartbnd = 0; ihighGlobal = 0; iamrGlobal = 0;
  }
 
  /** basic destructor */
  ~tioga(); 
  
  /** set communicator */
  void setCommunicator(MPI_Comm communicator,int id_proc,int nprocs);

  /** \brief Assign solver / grid / mesh data to TIOGA object
   *
   * @param[in] btag : body tag for rank
   * @param[in] nnodes : number of mesh nodes on rank
   * @param[in] xyz : x,y,z coordinates of each node
   * @param[inout] ibl : pointer to store nodal iblank data in (size: nnodes)
   * @param[in] nwbc : number of solid wall boundary nodes on rank
   * @param[in] nobc : number of overset boundary nodes on rank
   * @param[in] wbcnode : list of node IDs for wall boundary nodes
   * @param[in] obcnode : list of node IDs for overset boundary nodes
   * @param[in] ntypes : number of (volume) element types present on rank
   * @param[in] nv : number of vertices for each (volume) element type (size ntypes)
   * @param[in] nc : number of elements for each (volume) element type (size ntypes)
   * @param[in] vconn : element-to-vertex connectivity array */
  void registerGridData(int btag,int nnodes,double *xyz,int *ibl, int nwbc,int nobc,
			       int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn);

  /** Register additional face-connectivity arrays for Artificial Boundadry method */
  void registerFaceConnectivity(int nftype, int* nf, int *nfv, int** fconn, int *f2c, int *c2f);

  void profile(void);

  void exchangeBoxes(void);


  void exchangeSearchData(void);

  /** Send/Recv receptor points for each rank */
  void exchangePointSearchData(void);

  /** Perform overset connectivity: Determine receptor points & donor cells,
   *  and calculate interpolation weights for each receptor point */
  void exchangeDonors(void);

  /** Find donor element ID for given point */
  int findPointDonor(double* x_pt);

  /** Find all elements overlapping given bounding box */
  std::unordered_set<int> findCellDonors(double* bbox);

  /** perform overset grid connectivity (match receptor points to donor cells) */

  void performConnectivity(void);
  void performConnectivityHighOrder(void);
  void performConnectivityAMR(void);
  void performConnectivityArtificialBoundary(void);

  /** Perform overset interpolation and data communication */

  void dataUpdate(int nvar,double *q,int interptype) ;

  void dataUpdate_AMR(int nvar,double *q,int interptype) ;
  
  void dataUpdate_highorder(int nvar,double *q, int interptype) ;

  /** Perform data interpolation for artificial boundary method */
  void dataUpdate_artBnd(int nvar, double *q_spts, double* q_fpts, int interpType);

  /** get hole map for each mesh */
 
  void getHoleMap(void);

  /** output HoleMaps */
  
  void outputHoleMap(void);

  void writeData(int nvar,double *q,int interptype);

  void getDonorCount(int *dcount, int *fcount);
  
  void getDonorInfo(int *receptors,int *indices,double *frac,int *dcount);

  /** set symmetry bc */
  void setSymmetry(int syminput) { isym=syminput;}

  /** set resolutions for nodes and cells */
  void setResolutions(double *nres,double *cres)
  { mb->setResolutions(nres,cres);}

  void set_cell_iblank(int *iblank_cell)
  {
   mb->set_cell_iblank(iblank_cell);
  }

  //! Set callback functions for general high-order methods
  void setcallback(void (*f1)(int*, int*),
		    void (*f2)(int *,int *,double *),
		    void (*f3)(int *,double *,int *,double *),
		    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
		   void (*f5)(int *,int *,double *,int *,int*,double *))
  {
    mb->setcallback(f1,f2,f3,f4,f5);
    ihigh=1;
  }

  //! Set callback functions specific to Artificial Boundary method
  void set_ab_callback(void (*gnf)(int* id, int* npf),
                       void (*gfn)(int* id, int* npf, double* xyz),
                       void (*gqi)(int* id, int* fpt, int* ind, int* stride))
  {
    mb->setCallbackArtBnd(gnf, gfn, gqi);
    iartbnd = 1;
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
};
      
  


