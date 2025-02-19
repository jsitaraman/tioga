// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

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
