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

#ifdef _GPU
#include <cuda_runtime.h>
#endif

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
class Timer
{
private:
  std::chrono::high_resolution_clock::time_point tStart;
  std::chrono::high_resolution_clock::time_point tStop;
  double duration = 0;
  std::string prefix = "Execution Time = ";

public:
  Timer(void) {}

  Timer(std::string prefix) { this->prefix = prefix; }

  void setPrefix(std::string prefix) { this->prefix = prefix; }

  void startTimer(void)
  {
    tStart = std::chrono::high_resolution_clock::now();
  }

  void stopTimer(void)
  {
    tStop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( tStop - tStart ).count();
    duration += (double)elapsed/1000.;
  }

  void resetTimer(void)
  {
    duration = 0;
    tStart = std::chrono::high_resolution_clock::now();
  }

  double getTime(void)
  {
    return duration;
  }

  void showTime(int precision = 2)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    if (duration > 60) {
      int minutes = floor(duration/60);
      double seconds = duration-(minutes*60);
      std::cout << "Rank " << rank << ": " << prefix << minutes << "min ";
      std::cout << std::setprecision(precision) << seconds << "s" << std::endl;
    }
    else
    {
      std::cout << "Rank " << rank << ": " << prefix;
      std::cout << std::setprecision(precision) << duration << "s" << std::endl;
    }
  }
};

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
  parallelComm *pc_ab;
  int isym;
  int ierr;
  int mytag;
  int myid,nproc;
  int *sendCount;
  int *recvCount;
  OBB *obblist;
  int iorphanPrint;

  int nprocMesh, meshRank;

  int nCutFringe, nCutHole;   //! # of fringe/hole-cutting faces on this rank
  int gridType;               //! Type of grid: background (0) or normal (1)
  std::vector<int> gridIDs;   //! Grid ID for each rank
  std::vector<int> gridTypes; //! Grid type for each rank

  //! NEW - attempting to speed up code...
  std::vector<VPACKET> sndVPack, rcvVPack;
  std::vector<PACKET> sndPack2;
  std::vector<int> intData;
  std::vector<double> dblData;

  /* ---- GPU-Related Variables ---- */
#ifdef _GPU
  int resizeFlag = 0;
  int resizeFringe = 0;
  int ninterp = 0;
  int ninterp_d = 0;
  int nfringe_h = 0;
  double *ubuf_d = NULL;
  double *gradbuf_d = NULL;
  double *ubuf_h = NULL;
  double *gradbuf_h = NULL;
  double *fringebuf_h = NULL;
//  std::vector<double> dbuffer;
  std::vector<int> buf_disp, buf_inds;
  std::vector<int> recv_itmp;
  int *buf_inds_d = NULL;
#endif

 public:
  int ihigh;        /// High-Order flag for current rank
  int iartbnd;      /// Artificial-boundary flag for current rank
  int ihighGlobal;  /// Flag for whether high-order grids exist on any rank
  int iamrGlobal;   /// Flag for whether AMR cartesian grids exist on any rank
  int iabGlobal;    /// Flaag for whether high-order A.B.'s being used on any rank

  bool gpu = false; /// Flag for whether GPUs are being used on the high-order code

  Timer waitTime;
  Timer interpTime;

  /** basic constuctor */
  tioga()
  {
    mb = NULL; cg = NULL; cb = NULL;
    pc = NULL; sendCount = NULL; recvCount = NULL;
    holeMap = NULL; obblist = NULL;
    nblocks = 0; ncart = 0;
    isym = 2; ihigh = 0; iartbnd = 0; ihighGlobal = 0; iamrGlobal = 0;

    waitTime.setPrefix("TIOGA: MPI Time: ");
    interpTime.setPrefix("TIOGA: Total Interp Time: ");
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
  void registerFaceConnectivity(int nftype, int* nf, int *nfv, int** fconn,
                                int *f2c, int *c2f, int* iblank_face,
                                int nOverFaces, int nMpi, int* overFaces,
                                int* mpiFaces, int *procR, int *idR);

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
  void directCut(void);

  void setTransform(double *mat, double *off, int nDims) { mb->setTransform(mat,off,nDims); }

  /** Perform overset interpolation and data communication */

  void dataUpdate(int nvar,double *q,int interptype) ;

  void dataUpdate_AMR(int nvar,double *q,int interptype) ;
  
  void dataUpdate_highorder(int nvar,double *q, int interptype) ;

  /** Perform data interpolation for artificial boundary method */
  void dataUpdate_artBnd(int nvar, double *q_spts, int dataFlag);

  void dataUpdate_artBnd_send(int nvar, int dataFlag);
  void dataUpdate_artBnd_recv(int nvar, int dataFlag);

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
                       double (*gqs)(int ic, int spt, int var),
                       double& (*gqf)(int ff, int fpt, int var),
                       double (*ggs)(int ic, int spt, int dim, int var),
                       double& (*ggf)(int ff, int fpt, int dim, int var),
                       double* (*gqss)(int& es, int& ss, int& vs),
                       double* (*gdqs)(int& es, int& ss, int& vs, int& ds))
  {
    mb->setCallbackArtBnd(gnf, gfn, gqs, gqf, ggs, ggf, gqss, gdqs);
    iartbnd = 1;
  }

  void set_ab_callback_gpu(void (*d2h)(int* ids, int nd, int grad),
                           void (*h2d)(int* ids, int nf, int grad, double* data),
                           double* (*gqd)(int&, int&, int&),
                           double* (*gdqd)(int&, int&, int&, int&))
  {
    mb->setCallbackArtBndGpu(d2h,h2d,gqd,gdqd);
    gpu = true;
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

#ifdef _GPU
  void setupCommBuffersGPU(void);

  void set_stream_handle(cudaStream_t handle, cudaEvent_t event);
#endif
};
      
  


