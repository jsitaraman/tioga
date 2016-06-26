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
#include "ADT.h"

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
  int *nv;     /** < number of vertices for each type of cell */
  int *nc;     /** < number of each of different kinds of cells (tets, prism, pyramids, hex etc) */
  int *nf;     /** < number of faces for each cell type */
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

  std::vector<int> mpiFaces;  /** < List of MPI face IDs on rank */
  std::vector<int> mpiFidR;   /** < Matching MPI face IDs on opposite rank */
  std::vector<int> mpiProcR;  /** < Opposite rank for MPI face */

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

  /*!
   * \brief Get the location of 'q' (solution) data for the given face fpt
   *
   * 'ind' will be equivalent to calling std::distance between 'q_fpts(0,0,0)'
   * and 'q_fpts(faceID,fpt,0)'
   * 'stride' is equivalent to std::distance between 'q_fpts(faceID,fpt,var_0)'
   * and 'q_fpts(faceID,fpt,var_1)'
   *
   * @param[in] faceID   Rank-specific face ID
   * @param[in] fpt      Index of flux point on face [0 to nPointsFace-1]
   * @param[out] ind     Total starting index of solution data for point
   * @param[out] stride  Stride length to get from starting index to each
   *                     additional variable
   */
  void (*get_q_index_face)(int* faceID, int *fpt, int* ind, int* stride);

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
   * @param[in]  ndim      Amount of memory allocated to 'frac' (# of doubles)
   */
  void (*donor_frac)(int* cellID, double* xyz, int* nweights, int* inode, double* weights, double* rst, int* ndim);

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
  int ninterp2;            /** < number of interpolants for high-order points */
  int interp2ListSize;
  INTERPLIST *interpList2; /** < list for high-interpolation points */
  int ninterpCart;
  int interpListCartSize;
  INTERPLIST *interpListCart; 
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
    nfringe=2;
    // new vars
    ninterp=ninterp2=interpListSize=interp2ListSize=0;
    ctag=NULL;
    pointsPerCell=NULL;
    maxPointsPerCell=0;
    rxyz=NULL;
    ntotalPoints=0;
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

  void setFaceData(int* _nf, int *_f2v, int *_f2c, int *_c2f);

  void setResolutions(double *nres,double *cres);    
	       
  void search();

  /*! Given a 3D position, find the cell it lies within (-1 if not found) */
  int findPointDonor(double *x_pt);

  /*! Given a bounding box, find all elements which overlap with it */
  std::unordered_set<int> findCellDonors(double *bbox);

  void writeOBB(int bid);

  void writeOBB2(OBB *obc,int bid);

  void updateSolnData(int inode,double *qvar,double *q,int nvar,int interptype);

  int getNinterp(void) {return ninterp;};

  void getInterpolatedSolution(int *nints,int *nreals,int **intData,double **realData,double *q,
			       int nvar, int interptype);

  void getInterpolatedSolutionAMR(int *nints,int *nreals,int **intData,double **realData,double *q,
				  int nvar, int interptype);
  
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
  void getCellIblanks(void);

  //! Find all artificial boundary faces using previously-set cell iblank values
  void calcFaceIblanks(const MPI_Comm &meshComm);

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
  }

  //! Set callback functions specific to Artificial Boundary method
  void setCallbackArtBnd(void (*gnf)(int* id,int* npf),
                         void (*gfn)(int* id, int* npf, double* xyz),
                         void (*gqi)(int* id, int* fpt, int* ind, int* stride))
  {
    // See declaration of functions above for more details
    get_nodes_per_face = gnf;
    get_face_nodes = gfn;
    get_q_index_face = gqi;
  }

  void writeCellFile(int);

  /*! Gather a list of all receptor point locations (including for high-order) */
  void getInternalNodes(void);

  /*! Gather a list of all artificial boundary point locations (for high-order)
   * [Requires use of callback functions] */
  void getBoundaryNodes(void);

  void getExtraQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  void processPointDonors(void);
  void getInterpolatedSolutionAtPoints(int *nints,int *nreals,int **intData,
				       double **realData,
				       double *q,
				       int nvar, int interptype);

  /*! Update high-order element data at internal degrees of freedom */
  void updatePointData(double *q,double *qtmp,int nvar,int interptype);

  /*! Update high-order element data at artificial boundary flux points */
  void updateFluxPointData(double *q_fpts, double *qtmp, int nvar);

  void outputOrphan(FILE *fp,int i) 
  {
    fprintf(fp,"%f %f %f\n",rxyz[3*i],rxyz[3*i+1],rxyz[3*i+2]);
  }

  void outputOrphan(const std::fstream &fp, int i)
  {
    fp << rxyz[3*i] << " " << rxyz[3*i+1] << " " << rxyz[3*i+2] << std::endl;
  }

  void clearOrphans(int *itmp);
  void getUnresolvedMandatoryReceptors();
  void getCartReceptors(CartGrid *cg, parallelComm *pc);
  void setCartIblanks(void);

  /*! Apply blanking algorithm (to nodal iblank?) to get cell & face iblanks */
  void setArtificialBoundaries(void);
};
