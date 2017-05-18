#include "codetypes.h"
#include "error.hpp"
#include "dADT.h"
#include "MeshBlock.h"
#include "funcs.hpp"

#define MAX_LEVEL 32

dADT::~dADT()
{
//  clearData();
}

void dADT::clearData(void)
{
  adtInts.free_data();
  adtReals.free_data();
  adtBBox.free_data();
  coord.free_data();
}

void dADT::copyADT(ADT *adt)
{
  ndim = adt->ndim;
  nelem = adt->nelem;
  adtInts.assign(adt->adtIntegers, nelem*4);
  adtReals.assign(adt->adtReals, nelem*ndim);
  adtBBox.assign(adt->adtExtents, ndim);
  coord.assign(adt->coord, ndim*nelem);
  rrot = adt->rrot;
  if (rrot)
  {
    offset.assign(adt->offset.data(), adt->offset.size());
    Rmat.assign(adt->Rmat.data(), adt->Rmat.size());
  }
}

void dADT::setTransform(double *mat, double *off, int nDims)
{
  if (nDims != ndim/2)
    FatalError("dADT:setTransform:nDims != dADT::ndim/2");

  rrot = true;
  Rmat.assign(mat, nDims*nDims);
  offset.assign(off, nDims);
}

template<int nDims, int nside>
__device__
void d_searchADTstack(dADT& adt, dMeshBlock& mb, int& cellIndex, double* xsearch, double* rst)
{
  const int ndim = 2*nDims;
  int stack[MAX_LEVEL] = {0};
  int size = 1;

  while (size > 0)
  {
    int node = stack[size-1];
    size--;

    int ele = adt.adtInts[4*node];
    double bbox[ndim];
    for (int i = 0; i < ndim; i++)
      bbox[i] = mb.eleBBox[ndim*ele+i];

    bool flag = true;
    for (int i = 0; i < nDims; i++)
    {
      flag = (flag && (xsearch[i] >= bbox[i]-TOL));
      flag = (flag && (xsearch[i] <= bbox[i+nDims]+TOL));
    }

    if (flag)
    {
      mb.checkContainment<nDims,nside>(ele,cellIndex,bbox,xsearch,rst);
      if (cellIndex > -1)
        return;
    }

    // check the left and right children now
    for (int d = 1; d < 3; d++)
    {
      int nodeChild = adt.adtInts[4*node+d];
      if (nodeChild > -1)
      {
        nodeChild = adt.adtInts[4*nodeChild+3];
        for (int i = 0; i < ndim; i++)
          bbox[i] = adt.adtReals[ndim*nodeChild+i];

        flag = true;
        for (int i = 0; i < nDims; i++)
        {
          flag = (flag && (xsearch[i] >= bbox[i]-TOL));
          flag = (flag && (xsearch[i] <= bbox[i+nDims]+TOL));
        }

        if (flag)
          stack[size++] = nodeChild;
      }
    }
  }
}

template<int nDims, int nside>
__global__
void searchADT_kernel(dADT adt, dMeshBlock mb)
{
  int pt = blockDim.x * blockIdx.x + threadIdx.x;

  if (pt >= mb.nsearch) return;

  //const int ndim_adt = 2*nDims;

  // check if the given point is in the bounds of the ADT
  //int rootNode = 0;
  int cellID = -1;

  double xsearch[nDims];
  for (int d = 0; d < nDims; d++)
    xsearch[d] = mb.xsearch[nDims*pt+d];

  if (adt.rrot) // Transform back to ADT's coordinate system
  {
    double x2[nDims] = {0.0};
    for (int d1 = 0; d1 < nDims; d1++)
      for (int d2 = 0; d2 < nDims; d2++)
        x2[d1] += adt.Rmat[d1+nDims*d2]*(xsearch[d2]-adt.offset[d2]);

    for (int d = 0; d < nDims; d++)
      xsearch[d] = x2[d];
  }

  bool flag = true;
  for (int d = 0; d < nDims; d++)
  {
    flag = (flag && (xsearch[d] >= adt.adtBBox[2*d]-TOL));
    flag = (flag && (xsearch[d] <= adt.adtBBox[2*d+1]+TOL));
  }

  // call recursive routine to check intersections with ADT nodes
  double rst[nDims] = {0.0};
  if (flag)
    d_searchADTstack<nDims,nside>(adt,mb,cellID,xsearch,rst);

  __syncthreads();

  mb.donorId[pt] = cellID;
  for (int d = 0; d < nDims; d++)
    mb.rst[nDims*pt+d] = rst[d];
}

void searchADT(dADT &adt, dMeshBlock &mb)
{
  int threads = 32;
  int blocks = (mb.nsearch + threads - 1) / threads;

  switch (mb.nvert)
  {
    case 8:
      searchADT_kernel<3,2><<<blocks, threads>>>(adt, mb); break;
    case 27:
      searchADT_kernel<3,3><<<blocks, threads>>>(adt, mb); break;
    case 64:
      searchADT_kernel<3,4><<<blocks, threads>>>(adt, mb); break;
    default:
      ThrowException("nvert case not implemented");
  }

  check_error();
}
