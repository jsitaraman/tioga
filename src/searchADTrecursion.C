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
#include "MeshBlock.h"

void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,
			 double *coord,int level,int node,double *xsearch,int nelem,int ndim);

void searchBoxIntersections(int *elementList, std::unordered_set<int> &icells, int *adtIntegers,
  double *adtReals, double *coord, int level, int node, double *bbox, int nelem, int ndim);

void ADT::searchADT(MeshBlock *mb, int *cellIndex,double *xsearch)
{
  // check if the given point is in the bounds of the ADT
  int rootNode = 0;
  *cellIndex = -1;

  bool flag = true;
  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xsearch[i] >= adtExtents[2*i]-TOL));
  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xsearch[i] <= adtExtents[2*i+1]+TOL));

  // call recursive routine to check intersections with ADT nodes
  if (flag) searchIntersections(mb,cellIndex,adtIntegers,adtReals,
				coord,0,rootNode,xsearch,nelem,ndim);
}

void ADT::searchADT_box(int *elementList, std::unordered_set<int> &icells, double *bbox)
{
  int rootNode = 0;
  icells.clear();

  // Check if the given bounding box intersects with the the bounds of the ADT
  bool flag = true;
  for(int i = 0; i < ndim/2; i++)
  {
    flag = (flag && (bbox[i+ndim/2] >= adtExtents[2*i]  -TOL));
    flag = (flag && (bbox[i]   <= adtExtents[2*i+1]+TOL));
  }

  // Call recursive routine to check intersections with ADT nodes
  if (flag) searchBoxIntersections(elementList,icells,adtIntegers,adtReals,
                                   coord,0,rootNode,bbox,nelem,ndim);
}


void searchIntersections(MeshBlock *mb,int *cellIndex,int *adtIntegers,double *adtReals,
			 double *coord,int level,int node,double *xsearch,int nelem,int ndim)
{
  double element[ndim];
  bool flag = true;

  for (int i = 0; i < ndim; i++)
    element[i] = coord[ndim*(adtIntegers[4*node])+i];

  for (int i = 0; i < ndim/2; i++)
    flag = (flag && (xsearch[i] >= element[i]-TOL));

  for (int i = ndim/2; i < ndim; i++)
    flag = (flag && (xsearch[i-ndim/2] <= element[i]+TOL));

  if (flag)
  {
    mb->checkContainment(cellIndex,adtIntegers[4*node],xsearch);
    if (*cellIndex > -1) return;
  }

  // check the left and right children now
  for (int d = 1; d < 3; d++)
  {
    int nodeChild = adtIntegers[4*node+d];
    if (nodeChild > -1)
    {
      nodeChild = adtIntegers[4*nodeChild+3];
      for (int i = 0; i < ndim; i++)
        element[i] = adtReals[ndim*nodeChild+i];

      flag = true;
      for (int i = 0; i < ndim/2; i++)
        flag = (flag && (xsearch[i] >= element[i]-TOL));

      for (int i = ndim/2; i < ndim; i++)
        flag = (flag && (xsearch[i-ndim/2] <= element[i]+TOL));

      if (flag)
      {
        searchIntersections(mb,cellIndex,adtIntegers,adtReals,coord,level+1,
                            nodeChild,xsearch,nelem,ndim);
        if (*cellIndex > -1) return;
      }
    }
  }
  return;
}
  
void searchBoxIntersections(int *elementList, std::unordered_set<int> &icells, int *adtIntegers,
  double *adtReals, double *coord, int level, int node, double *bbox, int nelem, int ndim)
{
  double eleBox[ndim];
  for(int i = 0; i < ndim; i++)
    eleBox[i] = coord[ndim*(adtIntegers[4*node])+i];

  // Check if bbox intersects with bounding box of current mesh element in ADT
  bool flag = true;
  for (int i = 0; i < ndim/2; i++) {
    flag = (flag && (bbox[i+ndim/2] >= eleBox[i]  -TOL));
    flag = (flag && (bbox[i]   <= eleBox[i+ndim/2]+TOL));
  }

  if (flag)
  {
    int ind = elementList[adtIntegers[4*node]];
    icells.insert(ind);
  }

  // check the left and right children now
  for (int d = 1; d < 3; d++)
  {
    int nodeChild=adtIntegers[4*node+d];
    if (nodeChild > -1)
    {
      nodeChild = adtIntegers[4*nodeChild+3];

      double adtBox[ndim];
      for (int i=0; i<ndim; i++)
        adtBox[i] = adtReals[ndim*nodeChild+i];

      flag = true;
      for (int i = 0; i < ndim/2; i++)
      {
        flag = (flag && (bbox[i+ndim/2] >= adtBox[i] - TOL));
        flag = (flag && (bbox[i] <= adtBox[i+ndim/2] + TOL));
      }

      if (flag)
      {
        searchBoxIntersections(elementList,icells,adtIntegers,adtReals,coord,level+1,
                               nodeChild,bbox,nelem,ndim);
      }
    }
  }
}
