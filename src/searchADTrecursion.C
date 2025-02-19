// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"
#include "MeshBlock.h"

void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,
			 double *coord,int level,int node,double *xsearch,int nelem,int ndim);

void ADT::searchADT(MeshBlock *mb, int *cellIndex,double *xsearch)
{
  int i;
  int flag;
  int rootNode;
  //
  // check if the given point is in the bounds of
  // the ADT
  //
  rootNode=0;
  cellIndex[0]=-1;
  cellIndex[1]=0;
  //
  flag=1;
  for(i=0;i<ndim/2;i++)
    flag = (flag && (xsearch[i] >= adtExtents[2*i]-TOL));
  for(i=0;i<ndim/2;i++)
    flag= (flag && (xsearch[i] <= adtExtents[2*i+1]+TOL));
  //
  // call recursive routine to check intersections with 
  // ADT nodes
  //
  if (flag) searchIntersections(mb,cellIndex,adtIntegers,adtReals,
				coord,0,rootNode,xsearch,nelem,ndim);
}

void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,
			 double *coord,int level,int node,double *xsearch,int nelem,int ndim)
{
  int i;
  int d,nodeChild,dimcut;
  double element[ndim];
  bool flag;
  //
  for(i=0;i<ndim;i++)
    element[i]=coord[ndim*(adtIntegers[4*node])+i];
  //
  flag=1;
  for(i=0;i<ndim/2;i++)
    flag = (flag && (xsearch[i] >=element[i]-TOL));
  for(i=ndim/2;i<ndim;i++)
    flag = (flag && (xsearch[i-ndim/2] <=element[i]+TOL));
  //
  if (flag)
    {
      mb->checkContainment(cellIndex,adtIntegers[4*node],xsearch);
      if (cellIndex[0] > -1 && cellIndex[1]==0) return;
    }
  //
  // check the left and right children
  // now
  //
  for(d=1;d<3;d++)
    {
      nodeChild=adtIntegers[4*node+d];
      if (nodeChild > -1) {
        nodeChild=adtIntegers[4*nodeChild+3];
	for(i=0;i<ndim;i++)
         {
	  element[i]=adtReals[ndim*nodeChild+i];
         }
	flag=1;
	for(i=0;i<ndim/2;i++)
	  flag = (flag && (xsearch[i] >=element[i]-TOL));
	for(i=ndim/2;i<ndim;i++)
	  flag = (flag && (xsearch[i-ndim/2] <=element[i]+TOL));	
	if (flag)
	  {
	    searchIntersections(mb,cellIndex,adtIntegers,adtReals,coord,level+1,
			       nodeChild,xsearch,nelem,ndim);
	    if (cellIndex[0] > -1 && cellIndex[1]==0) return; 
	  }
      }
    }
  return;
}
  
