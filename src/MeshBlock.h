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
 * MeshBlock class - container and functions for generic unstructured grid partition in 3D
 *         
 * Jay Sitaraman
 * 02/20/2014
 */
#include <unordered_set>
#include <set>

#include "codetypes.h"
#include "funcs.hpp"
#include "points.hpp"
#include "ADT.h"

//! Helper struct for direct-cut method [Galbraith 2013]
typedef struct CutMap
{
  int type;                  //! Cut type: Solid wall (1) or overset bound (0)
  std::vector<int> flag;     //! Cut flag for all cells (essentially iblank)
  std::map<int,double> dist; //! Minimum distance to a cutting face
  std::map<int,int> nMin;    //! # of cut faces that are approx. 'dist' away
  std::map<int,Vec3> norm;   //! Normal vector of cutting face (or avg. of several)
  //std::map<int,Vec3> vec;    //! Vector from face to cell (between closest points)
} CutMap;

enum DIRECT_CUT_FLAG
{
  DC_HOLE, DC_UNASSIGNED, DC_CUT, DC_NORMAL
};

// forward declare to instantiate one of the methods
class parallelComm;
class CartGrid;

class MeshBlock
{
 private:
  int nnodes;  /** < number of grid nodes */
  int ncells;  /** < total number of cells */
  int nfaces;  /** < total number of faces (Art. Bnd.) */
  int ntypes;  /** < number of different types of cells */
  int nftype;  /** < number of different face types (triangle or quad) */
  int *nv;     /** < number of vertices for each type of cell */
  int *nc;     /** < number of each of different kinds of cells (tets, prism, pyramids, hex etc) */
  int *nf;     /** < number of faces for each cell type */
  int *nfv;    /** < number of vertices per face for each face type (3 or 4) */
  int nobc;    /** < number of overset boundary nodes */
  int nwbc;    /** < number of wall boundary nodes */
  //
  double *x;        /** < grid nodes x[3*nnodes] */
  int *iblank;      /** < iblank value for each grid node */
  int *iblank_cell; /** < iblank value at each grid cell */
  int *iblank_face; /** < iblank value at each grid face (Art. Bnd.) */
  //
  int **vconn;      /** < connectivity (cell to nodes) for each cell type */
  int **fconn;      /** < Connectivity (face to nodes) for each face type */
  int *f2c;         /** < Face to cell connectivity */
  int *c2f;         /** < Cell to face connectivity */
  int *wbcnode;     /** < wall boundary node indices */
  int *obcnode;     /** < overset boundary node indices */
  //
  double *nodeRes;  /** < node resolution  */
  double *userSpecifiedNodeRes;
  double *userSpecifiedCellRes;
  double *elementBbox; /** < bounding box of the elements */
  int *elementList;    /** < list of elements in */

  int nOverFaces;
  int nMpiFaces;
  int *overFaces; /** < List of explicitly-defined (by solver) overset faces */
  int *mpiFaces;  /** < List of MPI face IDs on rank */
  int *mpiFidR;   /** < Matching MPI face IDs on opposite rank */
  int *mpiProcR;  /** < Opposite rank for MPI face */

  std::vector<int> c2c;

  int nsend, nrecv;
  std::vector<int> sndMap, rcvMap;

  bool rrot = false;
  std::vector<double> Smat, offset;

  //
  // Alternating digital tree library
  //
  ADT *adt;   /** < Digital tree for searching this block */
  //
  DONORLIST **donorList;      /**< list of donors for the nodes of this mesh */
  //
  int ninterp;              /**< number of interpolations to be performed */
  int interpListSize;
  INTERPLIST *interpList;   /**< list of donor nodes in my grid, with fractions and information of
                                 who they donate to */ 
  int *interp2donor;

  INTEGERLIST *cancelList;  /** receptors that need to be cancelled because of */
  int ncancel;              /** conflicts with the state of their donors */

   /* ---- Callback functions for high-order overset connectivity ---- */

  /*!
   * \brief Get the number of solution points in given cell
   */
  void (*get_nodes_per_cell)(int* cellID, int* nNodes);

  /*!
   * \brief Get the number of flux points on given face
   *        For the artificial boundary method
   */
  void (*get_nodes_per_face)(int* faceID, int* nNodes);

  /*!
   * \brief Get the physical position of solution points in given cell
   *
   * input: cellID, nNodes
   * output: xyz [size: nNodes x 3, row-major]
   */
  void (*get_receptor_nodes)(int* cellID, int* nNodes, double* xyz);

  /*!
   * \brief Get the physical position of flux points on given face
   *        For the artificial boundary method
   *
   * @param[in]  faceID  Face ID within current rank
   * @param[in]  nNodes  Number of nodes expected on face
   * @param[out] xyz     Coordinates of each point on face [nNodes x 3, row-major]
   */
  void (*get_face_nodes)(int* faceID, int* nNodes, double* xyz);

  double (*get_q_spt)(int cellID, int spt, int var);

  double (*get_grad_spt)(int cellID, int spt, int dim, int var);
  double& (*get_grad_fpt)(int faceID, int fpt, int dim, int var);

  double& (*get_q_fpt)(int ff, int spt, int var);

  double* (*get_q_spts)(int& ele_stride, int& spt_stride, int& var_stride);
  double* (*get_dq_spts)(int& ele_stride, int& spt_stride, int& var_stride, int& dim_stride);

  // GPU-related functions
  double* (*get_q_spts_d)(int& ele_stride, int& spt_stride, int& var_stride);
  double* (*get_dq_spts_d)(int& ele_stride, int& spt_stride, int& var_stride, int& dim_stride);

  /*! Copy solution/gradient data for the donor elements from device to host */
  void (*data_from_device)(int* donorIDs, int nDonors, int gradFlag);

  /*! Copy updated solution/gradient for fringe faces from host to device */
  void (*data_to_device)(int* fringeIDs, int nFringe, int gradFlag, double *data);

  /*!
   * \brief Determine whether a point (x,y,z) lies within a cell
   *
   * Given a point's physical position, determine if it is contained in the
   * given cell; if so, return the reference coordinates of the point
   *
   * @param[in]  cellID    ID of cell within current mesh
   * @param[in]  xyz       Physical position of point to test
   * @param[out] passFlag  Is the point inside the cell? (no:0, yes:1)
   * @param[out] rst       Position of point within cell in reference coordinates
   */
  void (*donor_inclusion_test)(int* cellID, double* xyz, int* passFlag, double* rst);

  /*!
   * \brief Get interpolation weights for a point
   *
   * Get interpolation points & weights for current cell,
   * given a point in reference coordinates
   *
   * @param[in]  cellID    ID of cell within current mesh
   * @param[in]  xyz       Physical position of receptor point
   * @param[out] nweights  Number of interpolation points/weights to be used
   * @param[out] inode     Indices of donor points within global solution array
   * @param[out] weights   Interpolation weights for each donor point
   * @param[in]  rst       Reference coordinates of receptor point within cell
   * @param[in]  buffsize  Amount of memory allocated to 'weights' (# of doubles)
   */
  void (*donor_frac)(int* cellID, double* xyz, int* nweights, int* inode,
                     double* weights, double* rst, int* buffsize);

  void (*convert_to_modal)(int *,int *,double *,int *,int *,double *);

  int nreceptorCells;      /** number of receptor cells */
  int *ctag;               /** index of receptor cells */
  int *pointsPerCell;      /** number of receptor points per cell */
  int maxPointsPerCell;    /** max of pointsPerCell vector */

  /* ---- Artificial Boundary Variables ---- */
  int nreceptorFaces;      /** Number of artificial boundary faces */
  int *ftag;               /** Indices of artificial boundary faces */
  int *pointsPerFace;      /** number of receptor points per face */
  int maxPointsPerFace;    /** max of pointsPerFace vector */

  //std::unordered_set<int> overFaces; /** < List of Artificial Boundary face indices */

  double *rxyz;            /**  point coordinates */
  int ipoint; 
  int *picked;             /** < flag specifying if a node has been selected for high-order interpolation */

  int nreceptorCellsCart;
  int *ctag_cart;
  int *pickedCart;
 	
 public :
  int ntotalPointsCart;
  double *rxyzCart;
  int *donorIdCart;
  int donorListLength;

  int nfringe;
  int meshtag; /** < tag of the mesh that this block belongs to */
  //
  // oriented bounding box of this partition
  // 
  OBB *obb;
  //
  int nsearch;        /** < number of query points to search in this block */
  int *isearch;       /** < index of query points in the remote process */
  double *xsearch;    /** < coordinates of the query points */
  double *rst;            /**  natrural coordinates */
  int *donorId;       /** < donor indices for those found */
  int donorCount;
  int myid;
  double *cellRes;  /** < resolution for each cell */
  int ntotalPoints;        /**  total number of extra points to interpolate */
  int ihigh;               /** High-order flag for current rank */
  int iartbnd;             /** High-order artificial boundary flag for current rank */
  bool gpu = false;        /** Flag for GPUs being used on high-order solver */

  int ninterp2;            /** < number of interpolants for high-order points */
  int interp2ListSize;
  INTERPLIST *interpList2; /** < list for high-interpolation points */
  int ninterpCart;
  int interpListCartSize;
  INTERPLIST *interpListCart; 

  // Direct-Cut Method Variables - TO BE CLEANED UP
  int nDims = 3;
  int nGroups;
  int nCutFringe, nCutHole;   //! # of fringe/hole-cutting faces on this rank
  int gridType;               //! Type of grid: background (0) or normal (1)
  std::vector<int> cutFacesW, cutFacesO; //! Wall and overset cut face lists
  std::vector<std::vector<int>> cutFaces;  //! List of faces on each cut group
  std::vector<int> groupIDs;
  std::set<int> myGroups;
  std::vector<int> cutType_glob;
  std::vector<int> nGf; //! Number of faces on each cutting group

  /* ---- GPU-Related Variables ---- */
  int nSpts;  // Number of spts per ele on this rank
#ifdef _GPU
  double *weights_d = NULL;
  int *donors_d = NULL;
  int *buf_inds_d = NULL;
  std::vector<int> buf_inds, buf_disp;
#endif

  /** basic constructor */
  MeshBlock()
  {
    nv=NULL; nc=NULL; x=NULL;
    iblank=NULL; iblank_cell=NULL; iblank_face=NULL; vconn=NULL;
    wbcnode=NULL; obcnode=NULL;
    cellRes=NULL; nodeRes=NULL;
    elementBbox=NULL; elementList=NULL;
    adt=NULL; obb=NULL;
    donorList=NULL; interpList=NULL; interp2donor=NULL;
    nsearch=0; isearch=NULL; xsearch=NULL;
    donorId=NULL; cancelList=NULL;
    userSpecifiedNodeRes=NULL; userSpecifiedCellRes=NULL;
    nfringe=4;
    // new vars
    ninterp=ninterp2=interpListSize=interp2ListSize=0;
    ctag=NULL;
    ftag=NULL;
    pointsPerCell=NULL;
    pointsPerFace=NULL;
    maxPointsPerCell=0;
    rxyz=NULL;
    ntotalPoints=0;
    nreceptorFaces=0;
    rst=NULL;
    ihigh=0;
    ipoint=0;
    interpList2=NULL;
    picked=NULL;
    ctag_cart=NULL;
    rxyzCart=NULL;
    donorIdCart=NULL;
    pickedCart=NULL;
    ntotalPointsCart=0;
    nreceptorCellsCart=0;
    ninterpCart=0;
    interpListCartSize=0;
    interpListCart=NULL;
  }

  /** basic destructor */
  ~MeshBlock();
      
  void preprocess(void);

  void tagBoundary(void);
  
  void writeGridFile(int bid);

  void writeFlowFile(int bid,double *q,int nvar,int type);
  
  void setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
	       int *wbcnodei,int *obcnodei,
	       int ntypesi, int *nvi, int *nci, int **vconni);

  void setFaceData(int _nftype, int* _nf, int* _nfv, int** _f2v, int *_f2c,
                   int *_c2f, int* _ib_face, int nOver, int nMpi, int* oFaces,
                   int* mFaces, int* procR, int* idR);

  void setResolutions(double *nres,double *cres);    

  void setCommMap(int ns, int nr, int *sm, int *rm)
  {
    nsend = ns;  nrecv = nr;
    sndMap.assign(sm, sm+nsend);
    rcvMap.assign(rm, rm+nrecv);
  }

  void setTransform(double *mat, double* offset, int ndim);

  void search();

  /*! Given a 3D position, find the cell it lies within (-1 if not found) */
  int findPointDonor(double *x_pt);

  /*! Given a bounding box, find all elements which overlap with it */
  std::unordered_set<int> findCellDonors(double *bbox);

  void writeOBB(int bid);

  void writeOBB2(OBB *obc,int bid);

  void updateSolnData(int inode,double *qvar,double *q,int nvar,int interptype);

  int getNinterp(void) {return ninterp;}

  void getInterpolatedSolution(int *nints,int *nreals,int **intData,double **realData,double *q,
			       int nvar, int interptype);

  void getInterpolatedSolution2(int &nints,int &nreals,int *&intData, double *&realData,
                                double *q,int nvar, int interptype);

  void getInterpolatedSolutionAMR(int *nints,int *nreals,int **intData,double **realData,double *q,
				  int nvar, int interptype);
  
  void getInterpolatedSolutionArtBnd(int &nints, int &nreals,
           std::vector<int> &intData, std::vector<double> &realData, int nvar);

  void getInterpolatedGradientArtBnd(int& nints, int& nreals, std::vector<int>& intData, std::vector<double>& realData, int nvar);

  void checkContainment(int *cellIndex,int adtElement,double *xsearch);

  void getWallBounds(int *mtag,int *existWall, double wbox[6]);
  
  void markWallBoundary(int *sam,int nx[3],double extents[6]);

  void getQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  

  /** routines that do book keeping */

  void getDonorPacket(PACKET *sndPack, int nsend);

  void initializeDonorList();
  
  void insertAndSort(int pointid,int senderid,int meshtag, int remoteid, double donorRes);
  
  void processDonors(HOLEMAP *holemap, int nmesh,int **donorRecords,double **receptorResolution,
		     int *nrecords);

  void initializeInterpList(int ninterp_input);
  
  void findInterpData(int* recid,int irecord, double receptorRes);

  void findInterpListCart();

  void set_ninterp(int);

  void getCancellationData(int *nints, int **intData);

  void cancelDonor(int irecord);

  void getInterpData(int *nrecords, int **intData);

  void clearIblanks(void);

  void setIblanks(int inode);

  void getDonorCount(int *dcount,int *fcount);

  void getDonorInfo(int *receptors,int *indices, double *frac);
  //
  // routines for high order connectivity and interpolation
  //
  void getCellIblanks(const MPI_Comm meshComm);

  //! Find all artificial boundary faces using previously-set cell iblank values
  void calcFaceIblanks(const MPI_Comm &meshComm);

  void getCutGroupBoxes(std::vector<double> &cutBox, std::vector<std::vector<double>> &faceBox, int nGroups_glob);

  int getCuttingFaces(std::vector<double> &faceNodesW, std::vector<double> &faceNodesO);

  void getDirectCutCells(std::vector<std::unordered_set<int>> &cellList, std::vector<double> &cutBox_global, int nGroups_glob);

  //! Determine blanking status based upon given set of wall and overset faces
  void directCut(std::vector<double> &cutFaces, int nCut, int nvertf, CutMap& cutMap, int cutType = 1);

  //! Take the union of all cut flags
  void unifyCutFlags(std::vector<CutMap> &cutMap);

  void set_cell_iblank(int *iblank_cell_input)
  {
    iblank_cell=iblank_cell_input;
  }
  void setcallback(void (*f1)(int*, int*),
		    void (*f2)(int *,int *,double *),
		    void (*f3)(int *,double *,int *,double *),
		    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
		   void (*f5)(int *,int *,double *,int *,int*,double *))
  {
    get_nodes_per_cell=f1;
    get_receptor_nodes=f2;
    donor_inclusion_test=f3;
    donor_frac=f4;
    convert_to_modal=f5;

    ihigh = 1;
  }

  //! Set callback functions specific to Artificial Boundary method
  void setCallbackArtBnd(void (*gnf)(int* id, int* npf),
                         void (*gfn)(int* id, int* npf, double* xyz),
                         double (*gqs)(int ic, int spt, int var),
                         double& (*gqf)(int ff, int fpt, int var),
                         double (*ggs)(int ic, int spt, int dim, int var),
                         double& (*ggf)(int ff, int fpt, int dim, int var),
                         double* (*gqss)(int& es, int& ss, int& vs),
                         double* (*gdqs)(int& es, int& ss, int& vs, int& ds))
  {
    // See declaration of functions above for more details
    get_nodes_per_face = gnf;
    get_face_nodes = gfn;
    get_q_spt = gqs;
    get_q_fpt = gqf;
    get_grad_spt = ggs;
    get_grad_fpt = ggf;
    get_q_spts = gqss;
    get_dq_spts = gdqs;

    iartbnd = 1;
  }

  void setCallbackArtBndGpu(void (*d2h)(int* ids, int nd, int grad),
                            void (*h2d)(int* ids, int nf, int grad, double *data),
                            double* (*gqd)(int&, int&, int&),
                            double* (*gdqd)(int&, int&, int&, int&))
  {
    data_from_device = d2h;
    data_to_device = h2d;
    get_q_spts_d = gqd;
    get_dq_spts_d = gdqd;
    gpu = true;
  }


  void writeCellFile(int);

  /*! Gather a list of all receptor point locations (including for high-order) */
  void getInternalNodes(void);

  /*! Gather a list of all artificial boundary point locations (for high-order)
   * [Requires use of callback functions] */
  void getBoundaryNodes(void);

  void getExtraQueryPoints(OBB *obb, int &nints, int*& intData, int &nreals,
                           double*& realData);
  void processPointDonors(void);
  void getInterpolatedSolutionAtPoints(int *nints, int *nreals, int **intData,
      double **realData, double *q, int nvar, int interpType);

  void getInterpolatedGradientAtPoints(int &nints, int &nreals, int *&intData,
      double *&realData, double *q, int nvar);

  /*! Update high-order element data at internal degrees of freedom */
  void updatePointData(double *q,double *qtmp,int nvar,int interptype);

  /*! Update high-order element data at artificial boundary flux points */
  void updateFluxPointData(double *qtmp, int nvar);

  /*! Update solution gradient at artificial boundary flux points */
  void updateFluxPointGradient(double *dqtmp, int nvar);

  /*! Copy donor-cell data from device to host for use in interpolation */
  void getDonorDataGPU(int dataFlag = 0);

  /*! Copy fringe-face data back to device for computation in solver */
  void sendFringeDataGPU(int gradFlag);

  void outputOrphan(FILE *fp,int i) 
  {
    fprintf(fp,"%f %f %f\n",rxyz[3*i],rxyz[3*i+1],rxyz[3*i+2]);
  }

  void outputOrphan(std::ofstream &fp, int i)
  {
    fp << rxyz[3*i] << " " << rxyz[3*i+1] << " " << rxyz[3*i+2] << std::endl;
  }

  void clearOrphans(int *itmp);
  void getUnresolvedMandatoryReceptors();
  void getCartReceptors(CartGrid *cg, parallelComm *pc);
  void setCartIblanks(void);

  void setupADT(void);

  /*! Apply blanking algorithm (to nodal iblank?) to get cell & face iblanks */
  void setArtificialBoundaries(void);

  /*! Setup additional helpful connectivity structures */
  void extraConn(void);

  /* ---- GPU-Related Functions ---- */
#ifdef _GPU
  void setupBuffersGPU(int nsend, std::vector<int>& intData, std::vector<VPACKET>& sndPack);
  void interpSolution_gpu(double* q_out_d, int nvar);
  void interpGradient_gpu(double* dq_out_d, int nvar);
#endif
};
