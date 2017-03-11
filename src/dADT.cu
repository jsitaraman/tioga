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
printf("Hey! You shouldn't be clearing my data!\n"); /// DEBUGGING
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
  //coord.assign(adt->coord, ndim*nelem);
}

template<int level, int nDims, int nside>
__device__
void d_searchADTrecursion(dMeshBlock mb, int& cellIndex, int* adtIntegers, double* adtReals,
    double* coord, int node, double* element, double* xsearch, double* rst, int nelem);

template<> __device__
void d_searchADTrecursion<MAX_LEVEL,3,2>(dMeshBlock, int&, int*, double*, double*, int, double*, double*, double*, int)
{ printf("ERROR: Max recursion depth for d_searchADTrecursion exceeded!\n"); }

template<> __device__
void d_searchADTrecursion<MAX_LEVEL,3,3>(dMeshBlock, int&, int*, double*, double*, int, double*, double*, double*, int)
{ printf("ERROR: Max recursion depth for d_searchADTrecursion exceeded!\n"); }

template<> __device__
void d_searchADTrecursion<MAX_LEVEL,3,4>(dMeshBlock, int&, int*, double*, double*, int, double*, double*, double*, int)
{ printf("ERROR: Max recursion depth for d_searchADTrecursion exceeded!\n"); }

template<int level, int nDims, int nside>
__device__
void d_searchADTrecursion(dMeshBlock mb, int& cellIndex, int* adtIntegers, double* adtReals,
    double* coord, int node, double* element, double* xsearch, double* rst, int nelem)
{
  const int ndim = 2*nDims;
  bool flag = true;

  int ele = adtIntegers[4*node];
  for (int i = 0; i < ndim; i++)
    element[i] = coord[ndim*ele+i];

  for (int i = 0; i < nDims; i++)
  {
    flag = (flag && (xsearch[i] >= element[i]-TOL));
    flag = (flag && (xsearch[i] <= element[i+nDims]+TOL));
  }

  if (flag)
  {
    mb.checkContainment<nDims,nside>(ele,cellIndex,element,xsearch,rst);
    if (cellIndex > -1)
      return;
  }

  // check the left and right children now
  for (int d = 1; d < 3; d++)
  {
    int nodeChild = adtIntegers[4*node+d];
    if (nodeChild > -1)
    {
      nodeChild = adtIntegers[4*nodeChild+3];
      if (nodeChild < 0) printf("ERROR!!! Invalid ADT node number!!\n"); /// DEBUGGING
      for (int i = 0; i < ndim; i++)
        element[i] = adtReals[ndim*nodeChild+i];

      flag = true;
      for (int i = 0; i < nDims; i++)
      {
        flag = (flag && (xsearch[i] >= element[i]-TOL));
        flag = (flag && (xsearch[i] <= element[i+nDims]+TOL));
      }

      if (flag)
      {
        d_searchADTrecursion<level+1,nDims,nside>(mb,cellIndex,adtIntegers,adtReals,coord,
                            nodeChild,element,xsearch,rst,nelem);
        if (cellIndex > -1)
          return;
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
  int rootNode = 0;
  int cellID = -1;

  double xsearch[3];
  for (int d = 0; d < nDims; d++) /// TODO: templatize ndim
    xsearch[d] = mb.xsearch[3*pt+d];

  bool flag = true;
  for (int d = 0; d < nDims; d++)
  {
    flag = (flag && (xsearch[d] >= adt.adtBBox[2*d]-TOL));
    flag = (flag && (xsearch[d] <= adt.adtBBox[2*d+1]+TOL));
  }

  // call recursive routine to check intersections with ADT nodes
  double rst[nDims] = {0.0};
  if (flag)
  {
    double element[2*nDims];
    d_searchADTrecursion<0,nDims,nside>(mb,cellID,adt.adtInts.data(),adt.adtReals.data(),
      mb.eleBBox.data(),rootNode,element,xsearch,rst,adt.nelem);
  }

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
//      searchADT_kernel<3,2><<<blocks, threads, 0, mb.stream>>>(*this, mb); break;
    case 27:
      searchADT_kernel<3,3><<<blocks, threads, 0, mb.stream>>>(adt, mb); break;
    case 64:
      searchADT_kernel<3,4><<<blocks, threads, 0, mb.stream>>>(adt, mb); break;
    default:
      ThrowException("nvert case not implemented");
  }

  check_error();
}
