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

#ifndef ADT_H
#define ADT_H

#include <cstdlib>
#include <memory>

// forward declaration for instantiation
class MeshBlock; 

/**
 * Generic Alternating Digital Tree For Search Operations
 */
class ADT
{
  private :
  
  int ndim;          /** < number of dimensions (usually 3 but can be more) */
  int nelem;         /** < number of elements */
  int *adtIntegers;  /** < integers that define the architecture of the tree */
  double *adtReals;  /** < real numbers that provide the extents of each box */
  double *adtExtents; /** < global extents */
  double *coord;          /** < bounding box of each element */

 public :
  ADT() {ndim=6;nelem=0;adtIntegers=NULL;adtReals=NULL;adtExtents=NULL;coord=NULL;};
  ~ADT() 
    {
      if (adtIntegers) free(adtIntegers);
      if (adtReals) free(adtReals);
      if (adtExtents) free(adtExtents);
      adtIntegers=NULL;
      adtReals=NULL;
      adtExtents=NULL;
    };
  void clearData(void)
    {
      if (adtIntegers) free(adtIntegers);
      if (adtReals) free(adtReals);
      if (adtExtents) free(adtExtents);
      adtIntegers=NULL;
      adtReals=NULL;
      adtExtents=NULL;
    };      
  void buildADT(int d,int nelements,double *elementBbox);  
  void searchADT(MeshBlock *mb,int *cellindx,double *xsearch);
};


#endif /* ADT_H */
