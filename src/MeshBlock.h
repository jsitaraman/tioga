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

#ifndef MESHBLOCK_H
#define MESHBLOCK_H

#include <vector>
#include <stdint.h>
#include <assert.h>
#include "codetypes.h"
#include "ADT.h"
// forward declare to instantiate one of the methods
class parallelComm;
class CartGrid;

/**
 * MeshBlock class - container and functions for generic unstructured grid partition in 3D
 *
 * Jay Sitaraman
 * 02/20/2014
 */
class MeshBlock
{
 private:
  int nnodes;  /** < number of grid nodes */
  int ncells;  /** < total number of cells */
  int ntypes;  /** < number of different types of cells */
  int *nv;     /** < number of vertices for each types of cell */
  int *nc;     /**  < number of each of different kinds of cells (tets, prism, pyramids, hex etc) */
  int nobc;    /** < number of overset boundary nodes */
  int nwbc;    /** < number of wall boundary nodes */
  //
  double *x;        /** < grid nodes x[3*nnodes] */
  int *iblank_cell; /** < iblank value at each grid cell */
  //
  int **vconn;        /** < connectivity of each kind of cell */
  int *wbcnode;     /** < wall boundary node indices */
  int *obcnode;     /** < overset boundary node indices */
  uint64_t *cellGID;     /**< Global ID of the cell */
  uint64_t *nodeGID;     /**< Global ID for the nodes */
  //
  double *nodeRes;  /** < node resolution  */
  double *userSpecifiedNodeRes;
  double *userSpecifiedCellRes;
  double *elementBbox; /** < bounding box of the elements */
  int *elementList;    /** < list of elements in */
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
  void (*get_nodes_per_cell)(int*, int*);
  void (*get_receptor_nodes)(int *,int *,double *);
  void (*donor_inclusion_test)(int *,double *,int *,double *);
  void (*donor_frac)(int *,double *,int *,int *,double *,double *,int *);
  void (*convert_to_modal)(int *,int *,double *,int *,int *,double *);

  int nreceptorCells;      /** number of receptor cells */
  int *ctag;               /** index of receptor cells */
  int *pointsPerCell;      /** number of receptor points per cell */
  int maxPointsPerCell;     /** max of pointsPerCell vector */
  double *rxyz;            /**  point coordinates */
  int ipoint; 
  int *picked;             /** < flag specifying if a node has been selected for high-order interpolation */

  int nreceptorCellsCart;
  int *ctag_cart;
  int *pickedCart;
  int uniform_hex;
  double dx[3];
  double xlow[3];
  int idims[3];
  int *uindx;
  int *invmap;     // inverse map
  int *mapmask;    // mask
  int *icft; 	   // frequency table for nodal containment
  int mapdims[3];  // dimensions of the map
  double mapdx[3]; // sides of the map
 public :
  int *iblank;      /** < iblank value for each grid node */
  int *iblank_reduced;
  int ntotalPointsCart;
  double *rxyzCart;
  int *donorIdCart;
  int donorListLength;

  int nfringe;
  int mexclude;
  int meshtag; /** < tag of the mesh that this block belongs to */
  int check_uniform_hex_flag;
  double resolutionScale;
  //
  // oriented bounding box of this partition
  // 
  OBB *obb;
  OBB *obh;
  //
  int nsearch;        /** < number of query points to search in this block */
  int *isearch;       /** < index of query points in the remote process */
  int *tagsearch;       /** < index of query points in the remote process */
  double *res_search;   /** < resolution of search points */
  int *xtag;            /** < hash to determine if there are duplicates */
  double *xsearch;    /** < coordinates of the query points */
  double *rst;            /**  natrural coordinates */
  int *donorId;       /** < donor indices for those found */
  std::vector<uint64_t> gid_search; /**< Global node ID for the query points */
  int donorCount;
  int myid;
  double *cellRes;  /** < resolution for each cell */
  int ntotalPoints;        /**  total number of extra points to interpolate */
  int ihigh;
  int ninterp2;            /** < number of interpolants for high-order points */
  int interp2ListSize;
  INTERPLIST *interpList2; /** < list for high-interpolation points */
  int ninterpCart;
  int interpListCartSize;
  INTERPLIST *interpListCart; 
  int* receptorIdCart;

  //
  // call back functions to use p4est to search
  // its own internal data
  //
  void (*p4estsearchpt) (double *,int *,int *,int *);
  void (*check_intersect_p4est) (int *, int *);



  /** basic constructor */
  MeshBlock() { nv=NULL; nc=NULL; x=NULL;iblank=NULL;iblank_cell=NULL;vconn=NULL;wbcnode=NULL;
    obcnode=NULL; cellRes=NULL; nodeRes=NULL; elementBbox=NULL; elementList=NULL; adt=NULL; donorList=NULL;
    interpList=NULL; interp2donor=NULL; obb=NULL; nsearch=0; isearch=NULL; tagsearch=NULL;
    res_search=NULL;xsearch=NULL; donorId=NULL;xtag=NULL;
    adt=NULL; cancelList=NULL; userSpecifiedNodeRes=NULL; userSpecifiedCellRes=NULL; nfringe=1;
    mexclude=3;
    // new vars
    ninterp=ninterp2=interpListSize=interp2ListSize=0;
    ctag=NULL;pointsPerCell=NULL;maxPointsPerCell=0;rxyz=NULL;ntotalPoints=0;rst=NULL;ihigh=0;ipoint=0;
    interpList2=NULL;picked=NULL;ctag_cart=NULL;rxyzCart=NULL;donorIdCart=NULL;pickedCart=NULL;ntotalPointsCart=0;
    nreceptorCellsCart=0;ninterpCart=0;interpListCartSize=0;interpListCart=NULL;
    resolutionScale=1.0; receptorIdCart=NULL;

    cellGID = NULL;
    iblank_reduced=NULL;
    uniform_hex=0;
    check_uniform_hex_flag=0;
    uindx = NULL;
    obh   = NULL;
    invmap = NULL;
    icft   = NULL;
    mapmask= NULL;
  };

  /** basic destructor */
  ~MeshBlock();
      
  void preprocess(void);

  void tagBoundary(void);
  
  void writeGridFile(int bid);

  void writeFlowFile(int bid,double *q,int nvar,int type);
  
  void setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
	       int *wbcnodei,int *obcnodei,
               int ntypesi, int *nvi, int *nci, int **vconni,
               uint64_t* cell_gid=NULL, uint64_t* node_gid=NULL);


  void setResolutions(double *nres,double *cres);    
	       
  void search();
  void search_uniform_hex();
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

  void getQueryPoints2(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  void getQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  

  /** routines that do book keeping */

  void getDonorPacket(PACKET *sndPack, int nsend);

  void initializeDonorList();
  
  void insertAndSort(int pointid,int senderid,int meshtag, int remoteid, double donorRes,double receptorRes);
  
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
  
  void getStats(int mstat[2]);

  void setIblanks(int inode);

  void getDonorCount(int *dcount,int *fcount);

  void getDonorInfo(int *receptors,int *indices, double *frac);

  void getReceptorInfo(int *receptors);

  void getReducedOBB(OBB *,double *);
  void getReducedOBB2(OBB *,double *);

  void resetCoincident();
  //
  // routines for high order connectivity and interpolation
  //
  void getCellIblanks(void);
  void getCellIblanks2(void);
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

  void setp4estcallback(void (*f1)(double *,int *,int *,int *),
			void (*f2)(int *,int *))
  {
    p4estsearchpt=f1;
    check_intersect_p4est=f2;
  }

  void writeCellFile(int);
  void getInternalNodes(void);
  void getExtraQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  void processPointDonors(void);
  void getInterpolatedSolutionAtPoints(int *nints,int *nreals,int **intData,
				       double **realData,
				       double *q,
				       int nvar, int interptype);
  void updatePointData(double *q,double *qtmp,int nvar,int interptype);
  void outputOrphan(FILE *fp,int i) 
  {
    fprintf(fp,"%f %f %f\n",rxyz[3*i],rxyz[3*i+1],rxyz[3*i+2]);
  }
  void clearOrphans(HOLEMAP *holemap,int nmesh,int *itmp);
  void getUnresolvedMandatoryReceptors();
  void getCartReceptors(CartGrid *cg, parallelComm *pc);
  void setCartIblanks(void);
  
  // Getters
  inline int getMeshTag() const { return meshtag + (1 - BASE); }

  /**
   * Get donor packet for multi-block/partition setups
   *
   */
  void getMBDonorPktSizes(std::vector<int>&,
                          std::vector<int>&);

  void getMBDonorPackets(std::vector<int>&,
                         std::vector<int>&,
                         PACKET*);

  /** Reset interpolation list data structure
   *
   *  Reset the data structures in situations where the performConnectivity
   *  method is invoked at every timestep when meshes undergo relative motion.
   */
  void resetInterpData() {
    if (interpList) {
      for (int i = 0; i < interpListSize; i++) {
        if (interpList[i].inode)
          TIOGA_FREE(interpList[i].inode);
        if (interpList[i].weights)
          TIOGA_FREE(interpList[i].weights);
      }
      TIOGA_FREE(interpList);
    }
    ninterp = 0;
    interpListSize = 0;
  }
  void reduce_fringes() ;

  void check_for_uniform_hex();

  void create_hex_cell_map();
};

#endif /* MESHBLOCK_H */
