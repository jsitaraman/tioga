#ifndef DMESH_BLOCK_H
#define DMESH_BLOCK_H

#include "cuda_funcs.h"
#include <vector>

#define BIG_INT 2147483647

//! Helper struct for direct-cut method [Galbraith 2013]
typedef struct dCutMap
{
  int type;                  //! Cut type: Solid wall (1) or overset bound (0)
  dvec<int> flag;     //! Cut flag for all cells (essentially iblank)
  dvec<double> dist;  //! Minimum distance to a cutting face
  dvec<int> nMin;    //! # of cut faces that are approx. 'dist' away
  dvec<double> norm;   //! Normal vector of cutting face (or avg. of several)
  dvec<double> dot;   //! Dot prodcut of Normal vector with separation vector
} dCutMap;

class dMeshBlock
{
public:
  /* ------ Variables related to grid connectivity ------ */

  int nnodes;  //! Number of nodes in grid
  int ncells;  //! Number of cells/elements in grid
  int nc_adt;  //! Number of cells/elements in grid
  int nfaces;  //! Number of faces in grid (used for Art. Bnd.)
  int ntypes;  //! Number of different cell types present in grid
  int nftype;  //! Number of different face types present in grid
  int *nv;     //! Number of vertices for each type of cell
  int *nc;     //! Number of cells of each cell type
  int *nf;     //! Number of faces for each cell type
  int *nfv;    //! Number of vertices per face for each face type (3 or 4)
  int nobc;    //! Number of overset boundary nodes
  int nwbc;    //! Number of wall boundary nodes

  int rank;

  dvec<int> c2v;   //! Cell-to-vertex connectivity
  int** f2v;   //! Face-to-vertex connectivity
  int* f2c;    //! Face-to-cell connectivity
  int* c2f;    //! Cell-to-face connectivity
  int* wNodes; //! List of nodes on wall boundaries
  int* oNodes; //! List of nodes on pre-defined overset boundaries

  double *x;      //! Grid nodes coordinates [nnodes * ndim]
  double *coord;  //! Element node coordinates [ncells * nvert * ndim]

  int nDims = 3;

  int nvert;

  /* ------ Variables related to overset blanking ------ */

  dvec<int> iblank;
  int* iblank_cell;
  int* iblank_face;

  /* ------ Variables related to search operations ------ */

  dvec<int> eleList;     //! List of elements in d/ADT
  dvec<double> eleBBox;  //! Bounding box of elements in d/ADT

  int nsearch;
  int donorCount;
  dvec<int> isearch;
  dvec<double> xsearch;
  dvec<double> rst;
  dvec<int> donorId;

  /* ------ Auxiliary Geometry Representation Variables ------ */

  dvec<char> hm_sam;        //! Structured Auxiliary Map [hole map]
  dvec<int> hm_nx;          //! # of Cartesian cells in each direction of map
  dvec<double> hm_dx;       //! Size of Cartesian cells (extents / nx)
  dvec<double> hm_extents;  //! Hole map bounding extents [bounding box]

  /* ------ Misc. Variables ------ */

  cudaStream_t stream;

  dvec<int> ijk2gmsh;
  dvec<int> ijk2gmsh_quad; // For Direct Cut method
  dvec<double> xlist;
  dvec<float> xlistf;

  bool rrot = false;
  dvec<double> Rmat, offset;


  /* ------ Direct Cut Variables ------ */

  hvec<int> cutFlag_h;
  dvec<int> cutFlag_d;
  dvec<int> filt_eles;
  dvec<int> filt_faces;
  dvec<float> ele_bbox;
  dvec<float> face_bbox;

  /* ------ Member Functions ------ */

  dMeshBlock() { }

  void dataToDevice(int ndims, int nnodes, int ncells, int ncells_adt, int nsearch, int* nv,
      int* nc, int* eleList, double* eleBBox, int* isearch, double* xsearch, int rank);

  void extraDataToDevice(int* vconn);

  void assignHoleMap(bool hasWall, int* nx, int* sam, double* extents);

  void clearHoleMap(void);

  void setDeviceData(double* vx, double* ex, int* ibc, int* ibf);

  void setTransform(double *mat, double *off, int ndim);

  void updateADTData(int ncells_adt, int* eleList, double* eleBBox);

  void updateSearchPoints(int nsearch, int* isearch, double* xsearch);

  template<int ndim, int nside>
  __device__
  void checkContainment(int adtEle, int& cellID, const double* __restrict__ bbox,
      const double* __restrict__ xsearch, double* __restrict__ rst);

  template<int nSide>
  __device__ __forceinline__
  void calcDShape(double* __restrict__ shape,
      double* __restrict__ dshape, const double* loc);

  template<int nSide>
  __device__ __forceinline__
  bool getRefLoc(const double* __restrict__ coords, const double* __restrict__ bbox,
                 const double* __restrict__ xyz, double* __restrict__ rst);

  template<int ndim, int nside>
  __device__
  void checkContainment(int adtEle, int& cellID, const float* __restrict__ bbox,
      const float* __restrict__ xsearch, float* __restrict__ rst);

  template<int nSide>
  __device__ __forceinline__
  void calcDShape(float* __restrict__ shape,
      float* __restrict__ dshape, const float* loc);

  template<int nSide>
  __device__ __forceinline__
  int getRefLoc(const float* __restrict__ coords, const float* __restrict__ bbox,
                 const float* __restrict__ xyz, float* __restrict__ rst);

  __host__
  void directCut(double* cutFaces_h, int nCut, int nvertf, double* cutBbox, int* cutFlag, int cutType);
};

#ifdef __CUDACC__

template<int nSide>
__device__ __forceinline__
void dMeshBlock::calcDShape(double* __restrict__ shape, double* __restrict__ dshape,
                            const double* loc)
{
  double xi = loc[0];
  double eta = loc[1];
  double mu = loc[2];

  double lag_i[nSide];
  double lag_j[nSide];
  double lag_k[nSide];
  double dlag_i[nSide];
  double dlag_j[nSide];
  double dlag_k[nSide];

  for (int i = 0; i < nSide; i++)
  {
    lag_i[i] = cuda_funcs::Lagrange_gpu(xlist.data(), nSide,  xi, i);
    lag_j[i] = cuda_funcs::Lagrange_gpu(xlist.data(), nSide, eta, i);
    lag_k[i] = cuda_funcs::Lagrange_gpu(xlist.data(), nSide,  mu, i);
    dlag_i[i] = cuda_funcs::dLagrange_gpu(xlist.data(), nSide,  xi, i);
    dlag_j[i] = cuda_funcs::dLagrange_gpu(xlist.data(), nSide, eta, i);
    dlag_k[i] = cuda_funcs::dLagrange_gpu(xlist.data(), nSide,  mu, i);
  }

  for (int k = 0; k < nSide; k++)
    for (int j = 0; j < nSide; j++)
      for (int i = 0; i < nSide; i++)
      {
        int gnd = ijk2gmsh[i+nSide*(j+nSide*k)];
        shape[gnd] = lag_i[i] * lag_j[j] * lag_k[k];
        dshape[gnd*3+0] = dlag_i[i] *  lag_j[j] *  lag_k[k];
        dshape[gnd*3+1] =  lag_i[i] * dlag_j[j] *  lag_k[k];
        dshape[gnd*3+2] =  lag_i[i] *  lag_j[j] * dlag_k[k];
      }
}

template<int nSide>
__device__ __forceinline__
void dMeshBlock::calcDShape(float* __restrict__ shape, float* __restrict__ dshape,
                            const float* loc)
{
  float xi = loc[0];
  float eta = loc[1];
  float mu = loc[2];

  float lag_i[nSide];
  float lag_j[nSide];
  float lag_k[nSide];
  float dlag_i[nSide];
  float dlag_j[nSide];
  float dlag_k[nSide];

  for (int i = 0; i < nSide; i++)
  {
    lag_i[i] = cuda_funcs::Lagrange_gpu(xlistf.data(), nSide,  xi, i);
    lag_j[i] = cuda_funcs::Lagrange_gpu(xlistf.data(), nSide, eta, i);
    lag_k[i] = cuda_funcs::Lagrange_gpu(xlistf.data(), nSide,  mu, i);
    dlag_i[i] = cuda_funcs::dLagrange_gpu(xlistf.data(), nSide,  xi, i);
    dlag_j[i] = cuda_funcs::dLagrange_gpu(xlistf.data(), nSide, eta, i);
    dlag_k[i] = cuda_funcs::dLagrange_gpu(xlistf.data(), nSide,  mu, i);
  }

  for (int k = 0; k < nSide; k++)
    for (int j = 0; j < nSide; j++)
      for (int i = 0; i < nSide; i++)
      {
        int gnd = ijk2gmsh[i+nSide*(j+nSide*k)];
        shape[gnd] = lag_i[i] * lag_j[j] * lag_k[k];
        dshape[gnd*3+0] = dlag_i[i] *  lag_j[j] *  lag_k[k];
        dshape[gnd*3+1] = lag_i[i] * dlag_j[j] *  lag_k[k];
        dshape[gnd*3+2] = lag_i[i] *  lag_j[j] * dlag_k[k];
      }
}

template<int nSide>
__device__
bool dMeshBlock::getRefLoc(const double* __restrict__ coords,
    const double* __restrict__ bbox, const double* __restrict__ xyz,
    double* __restrict__ rst)
{
  const int nNodes = nSide*nSide*nSide;

  // Use a relative tolerance to handle extreme grids
  float h = ( (bbox[3]-bbox[0]) + (bbox[4]-bbox[1]) + (bbox[5]-bbox[2]) ) / 3.;

  float EPS = 1e-5f;
  float tol = EPS*h;

  int iter = 0;
  int iterMax = 10;
  float norm = 1;
  float norm_prev = 2;

  double shape[nNodes];
  double dshape[3*nNodes];

  rst[0] = 0.;
  rst[1] = 0.;
  rst[2] = 0.;

  while (norm > tol && iter < iterMax)
  {
    calcDShape<nSide>(shape, dshape, rst);

    float dx[3] = {(float)xyz[0], (float)xyz[1], (float)xyz[2]};
    float grad[3][3] = {{0.0}};

    for (int nd = 0; nd < nNodes; nd++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          grad[i][j] += coords[i+3*nd] * dshape[nd*3+j];

    for (int nd = 0; nd < nNodes; nd++)
      for (int i = 0; i < 3; i++)
        dx[i] -= shape[nd] * coords[i+3*nd];

    float idetJ = 1.f / cuda_funcs::det_3x3(&grad[0][0]);

    cuda_funcs::adjoint_3x3(&grad[0][0]);

    float delta[3] = {0.0};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        delta[i] += grad[i][j]*dx[j] * idetJ;

    norm = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int i = 0; i < 3; i++)
      rst[i] = max(min(rst[i]+delta[i],1.+1e-4),-1.-1e-4);

    if (iter > 1 && norm > .99*norm_prev) // If it's clear we're not converging
      break;

    norm_prev = norm;

    iter++;
  }

  if ( norm <= 10.*tol || (norm <= 20.*tol && (max(max(rst[0],rst[1]),rst[2])) <= 1.0) )
    return true;
  else
    return false;
}

template<int nSide>
__device__
int dMeshBlock::getRefLoc(const float* __restrict__ coords,
    const float* __restrict__ bbox, const float* __restrict__ xyz,
    float* __restrict__ rst)
{
  const int nNodes = nSide*nSide*nSide;

  // Use a relative tolerance to handle extreme grids
  const float h = ( (bbox[3]-bbox[0]) + (bbox[4]-bbox[1]) + (bbox[5]-bbox[2]) ) / 3.f;
  const double tol = 1e-4*h;

  int iter = 0;
  int iterMax = 10;
  double norm = 1;
  double norm_prev = 2;

  double shape[nNodes];
  double dshape[3*nNodes];

  rst[0] = 0.f;
  rst[1] = 0.f;
  rst[2] = 0.f;
  double rstd[3] = {0.0};

  while (norm > tol && iter < iterMax)
  {
    calcDShape<nSide>(shape, dshape, rstd);

    double dx[3] = {xyz[0], xyz[1], xyz[2]};
    double grad[3][3] = {{0.0}};

    for (int nd = 0; nd < nNodes; nd++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          grad[i][j] += coords[i+3*nd] * dshape[nd*3+j];

    for (int nd = 0; nd < nNodes; nd++)
      for (int i = 0; i < 3; i++)
        dx[i] -= shape[nd] * coords[i+3*nd];

    double idetJ = 1. / cuda_funcs::det_3x3(&grad[0][0]);

    cuda_funcs::adjoint_3x3(&grad[0][0]);

    double delta[3] = {0.0};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        delta[i] += grad[i][j]*dx[j] * idetJ;

    norm = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int i = 0; i < 3; i++)
      rstd[i] = max(min(rstd[i]+delta[i],1.0005),-1.0005);

    if (iter > 1 && norm > .99*norm_prev) // If it's clear we're not converging
      break;

    norm_prev = norm;

    iter++;
  }

  rst[0] = rstd[0];
  rst[1] = rstd[1];
  rst[2] = rstd[2];

  if (norm <= 10.*tol) // Return 'yes'
    return 1;
  else if (norm <= 100.*tol) // Return 'maybe'
    return -1;
  else
    return 0;
}

template<int ndim, int nside>
__device__
void dMeshBlock::checkContainment(int adtEle, int& cellID,
    const float* __restrict__ bbox, const float* __restrict__ xyz,
    float* __restrict__ rst)
{
  const int nNodes = nside*nside*nside;

  int ele = eleList[adtEle];
  cellID = -BIG_INT;

  float ecoord[nNodes*ndim];
  for (int i = 0; i < nNodes; i++)
    for (int d = 0; d < ndim; d++)
      ecoord[i*ndim+d] = coord[ele+ncells*(d+ndim*i)];

  int isInEle = 0;

  if (rrot) // Transform search point *back* to current *physical* location
  {
    float x2[ndim];
    for (int d1 = 0; d1 < ndim; d1++)
    {
      x2[d1] = offset[d1];
      for (int d2 = 0; d2 < ndim; d2++)
        x2[d1] += Rmat[d1*ndim+d2] * xyz[d2];
    }

    isInEle = getRefLoc<nside>(ecoord,bbox,x2,rst);
  }
  else
  {
    isInEle = getRefLoc<nside>(ecoord,bbox,xyz,rst);
  }

  if (isInEle > 0) cellID = ele;
  else if (isInEle < 0) cellID = -ele;
}

template<int ndim, int nside>
__device__
void dMeshBlock::checkContainment(int adtEle, int& cellID, const double* __restrict__ bbox,
    const double* __restrict__ xyz, double* __restrict__ rst)
{
  const int nNodes = nside*nside*nside;

  int ele = eleList[adtEle];
  cellID = -BIG_INT;

  double ecoord[nNodes*ndim];
  for (int i = 0; i < nNodes; i++)
    for (int d = 0; d < ndim; d++)
      ecoord[i*ndim+d] = coord[ele+ncells*(d+ndim*i)];

  int isInEle = 0;

  if (rrot) // Transform search point back to current physical location
  {
    double x2[ndim];
    for (int d1 = 0; d1 < ndim; d1++)
    {
      x2[d1] = offset[d1];
      for (int d2 = 0; d2 < ndim; d2++)
        x2[d1] += Rmat[d1*ndim+d2] * xyz[d2];
    }

    isInEle = getRefLoc<nside>(ecoord,bbox,x2,rst);
  }
  else
  {
    isInEle = getRefLoc<nside>(ecoord,bbox,xyz,rst);
  }

  if (isInEle) cellID = ele;
}

#endif

#endif
