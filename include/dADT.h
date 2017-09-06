#ifndef DADT_H
#define DADT_H
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
#include "cuda_funcs.h"
#include "dMeshBlock.h"
#include "ADT.h"

/*! Alternating Digital Tree For Search Operations - GPU version */
class ADT;
class dADT
{
public:
  int ndim = 6;          /** < number of dimensions (usually 3 but can be more) */
  int nelem = 0;         /** < number of elements */

  dvec<int> adtInts;     //! Graph connectivity of ADT
  dvec<double> adtReals; //! Extents of each ADT node's region (layout: xmin,ymin,xmax,ymax)
  dvec<double> adtBBox;  //! Global bounding box of ADT (layout: xmin,xmax,ymin,ymax)
  dvec<double> coord;    //! Actual bounding box of each element in ADT

  dvec<double> offset;
  dvec<double> Rmat;

  bool rrot = false;         /** Flag for rigid-body rotation (apply transform to all search points) */

  dADT() { }

  ~dADT();

  void copyADT(ADT *adt);

  void clearData(void);

  void setTransform(double* mat, double* off, int nDims);
};

//! Search the ADT for the element containint the point xsearch
void searchADT(dADT &adt, dMeshBlock &mb);
#endif
