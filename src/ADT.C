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
 * Build an alternating digital tree 
 */
#include <stdio.h>
#include <stdlib.h>
#include "codetypes.h"
#include "ADT.h"
extern "C"{
void buildADTrecursion(double *coord,double *adtReals,double *adtWork,int *adtIntegers,
		       int *elementsAvailable,int *adtCount,int side,int parent,
		       int level,int ndim,int nelem, int nav);}

extern void median_(int *,double *,int *,double *);

void ADT::buildADT(int d, int nelements, double *elementBbox)
{
  /* set dimensions and number of elements */

  ndim = d;
  nelem = nelements;

  /* set element bbox pointer */

  coord = elementBbox;

  /* Allocate work arrays */

  int *elementsAvailable = (int *) malloc(sizeof(int)*nelem);
  double *adtWork = (double *) malloc(sizeof(double)*nelem);

  /* Allocate arrays in the class */

  if (adtExtents) free(adtExtents);
  adtExtents = (double *) malloc(sizeof(double)*ndim);
  if (adtIntegers) free(adtIntegers);
  adtIntegers = (int *) malloc(sizeof(int)*4*nelem);
  if (adtReals) free(adtReals);
  adtReals = (double *) malloc(sizeof(double)*nelem*ndim);

  /* Determine extent of elements */

  for (int i = 0; i < ndim/2; i++)
  {
    int i2=2*i;
    adtExtents[i2] = BIGVALUE;
    adtExtents[i2+1] = -BIGVALUE;
  }

  for (int j = 0; j < nelem; j++)
  {
    int jd = ndim*j;
    for (int i = 0; i < ndim/2; i++)
    {
      int i2 = 2*i;
      adtExtents[i2] = min(adtExtents[i2],coord[jd+i]);
    }
    for (int i = 0; i < ndim/2; i++)
    {
      int i2 = 2*i+1;
      adtExtents[i2] = max(adtExtents[i2],coord[jd+i+ndim/2]);
    }
  }

  // make the extents 1% larger

  double tolerance = 0.01;
  for (int i = 0; i < ndim/2; i++)
  {
    int i2 = 2*i;
    double delta=tolerance*(adtExtents[i2+1]-adtExtents[i2]);
    adtExtents[i2] -= delta;
    adtExtents[i2+1] += delta;
  }

  // Build ADT using a recursive process now

  for (int i = 0; i < nelem; i++)
    elementsAvailable[i] = i;

  // set initialvalues

  int adtCount = -1;
  int side = 0;
  int parent = 0;
  int level = 0;
  int nav = nelem;

  buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
                    &adtCount,side,parent,level,ndim,nelem,nav);

  // create Inverse map [ADT index <== original ele ID]
  // adtInt[eleID+3] = adtInd
  for (int i = 0; i < nelem; i++)
  {
    int eleID = 4*adtIntegers[4*i];
    adtIntegers[eleID+3] = i;
  }

  free(elementsAvailable);
  free(adtWork);
}
