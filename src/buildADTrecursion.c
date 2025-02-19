// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"

extern void median_(int *,double *,int *,double *);

void buildADTrecursion(double *coord,double *adtReals,double *adtWork,int *adtIntegers,
		       int *elementsAvailable,int *adtCount,int side,int parent,
		       int level,int ndim,int nelem, int nav)
{
  
  int nd=ndim/2;  
  double coordmid;
  int i,j;
  int dimcut;
  int nleft;
  int ii,iip,jj,jjp;
  int parentToChild;

  if (nav > 1) {
    //
    // find the dimension to create the cut
    //
    dimcut=(level%ndim);
    //
    // collect coordinates along the dimension dimcut
    //
    for(i=0;i<nav;i++)
      adtWork[i]=coord[ndim*elementsAvailable[i]+dimcut];
    //
    // reorder elements with nleft elements to
    // the left of median of adtWork
    //
    median_(elementsAvailable,adtWork,&nav,&coordmid);
    nleft=(nav+1)/2;
    (*adtCount)++;
    ii=(*adtCount)*4;
    adtIntegers[ii]=elementsAvailable[nleft-1];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;
    //
    // find minimum and maximum bounds of the elements
    // contained in this leaf
    //
    for(i=0;i<nd;i++)
      {
	adtReals[ndim*(*adtCount)+i]=BIGVALUE;
	adtReals[ndim*(*adtCount)+i+nd]=-BIGVALUE;
      }
    //
    for(i=0;i<nav;i++)
      for(j=0;j<nd;j++)
	{
	  ii=ndim*(*adtCount)+j;
	  iip=ii+nd;
	  jj=ndim*elementsAvailable[i]+j;
	  jjp=jj+nd;
	  //
	  adtReals[ii]=TIOGA_MIN(adtReals[ii],coord[jj]);
	  adtReals[iip]=TIOGA_MAX(adtReals[iip],coord[jjp]);
	}
    //
    // specify that the new element is the child of parent
    // unless root
    //
    if (side > 0) 
      {
	adtIntegers[4*parent+side]=elementsAvailable[nleft-1];
      }
    parentToChild=*adtCount;
    //
    // build the left side of the tree
    //
    if (nleft > 1) {
      buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
			adtCount,1,parentToChild,level+1,ndim,nelem,nleft-1);
    }
    //
    // build the right side of the tree
    //
    buildADTrecursion(coord,adtReals,adtWork,adtIntegers,&(elementsAvailable[nleft]),
		      adtCount,2,parentToChild,level+1,ndim,nelem,nav-nleft);
  }
  else if (nav==1) {
    (*adtCount)++;
    ii=4*(*adtCount);
    jj=ndim*(*adtCount);
    adtIntegers[ii]=elementsAvailable[0];
    adtIntegers[ii+1]=-1;
    adtIntegers[ii+2]=-1;
    adtIntegers[ii+3]=-1;
    for(j=0;j<ndim;j++)
      adtReals[jj+j]=coord[ndim*elementsAvailable[0]+j];
    if (side > 0) {
      adtIntegers[4*parent+side]=elementsAvailable[0];
    }
  }
}

