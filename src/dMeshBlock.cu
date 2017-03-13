#include "dMeshBlock.h"
#include "funcs.hpp"

/* ------ dMeshBlock Member Functions ------ */

void dMeshBlock::dataToDevice(int ndims, int nnodes, int ncells, int ncells_adt,
    int nsearch, int* nv, int* nc, int* eleList, double* eleBBox, int* isearch,
    double* xsearch)
{
  this->nnodes = nnodes;
  this->ncells = ncells;
  this->nc_adt = ncells_adt;

  this->nv = nv;
  this->nc = nc;

  nvert = nv[0];

  this->eleBBox.assign(eleBBox, ncells_adt*ndims*2);
  this->eleList.assign(eleList, ncells_adt);

  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);

  auto ijk2gmsh_h = tg_funcs::structured_to_gmsh_hex(nvert);
  ijk2gmsh.assign(ijk2gmsh_h.data(), ijk2gmsh_h.size());

  int nSide = std::cbrt(nvert);
  std::vector<double> xlist_h(nSide);
  double dxi = 2./(nSide-1);

  for (int i = 0; i < nSide; i++)
    xlist_h[i] = -1. + i*dxi;

  xlist.assign(xlist_h.data(), xlist_h.size());
}

void dMeshBlock::updateSearchPoints(int nsearch, int *isearch, double *xsearch)
{
  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);
}

void dMeshBlock::setDeviceData(double* vx, double* ex, int* ibc, int* ibf)
{
  x = vx;
  iblank_cell = ibc;
  iblank_face = ibf;
  coord = ex;
}

void dMeshBlock::setTransform(double* mat, double* off, int ndim)
{
  if (ndim != nDims)
    ThrowException("dMeshBlock::set_transform: input ndim != nDims");

  rrot = true;
  Rmat.assign(mat, ndim*ndim);
  offset.assign(off, ndim);
}
