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

#ifndef CARTBLOCK_H
#define CARTBLOCK_H

#include <cstdlib>
#include "codetypes.h"

struct INTERPLIST2;
struct DONORLIST;
struct HOLEMAP;

class CartGrid;
class CartBlock
{
 private:
  int local_id;
  int global_id;
  int dims[3],nf,qstride,ndof,pdegree,p3;
  int d1,d2,d3,d3nf;
  int myid;
  int *ibl;
  double *q;
  double *qnode;
  double xlo[3]; 
  double dx[3];
  int ndonors;
  int interpListSize;
  INTERPLIST2 *interpList,*listptr;
  DONORLIST **donorList;
  void (*donor_frac) (int *,double *,int *,double *);
 public:
  CartBlock() { global_id=0;dims[0]=dims[1]=dims[2]=0;ibl=NULL;q=NULL;interpListSize=0;donorList=NULL;interpList=NULL;
    donor_frac=nullptr;};
  ~CartBlock() { clearLists();};
  void registerData(int local_id_in,int global_id_in,int *iblankin,double *qin)
  {
    local_id=local_id_in;
    global_id=global_id_in;
    ibl=iblankin;
    q=qin;
  };
  void preprocess(CartGrid *cg);
  void getInterpolatedData(int *nints,int *nreals,int **intData,
			   double **realData,
			   int nvar);
  void update(double *qval,int index,int nq);
  void getCancellationData(int *cancelledData, int *ncancel);
  void processDonors(HOLEMAP *holemap, int nmesh);
  void insertInDonorList(int senderid,int index,int meshtagdonor,int remoteid,int remoteblockid,double cellRes);
  void insertInInterpList(int procid,int remoteid,int remoteblockid,double *xtmp);
  void writeCellFile(int bid);
  void clearLists(void);
  void initializeLists(void);
};

#endif /* CARTBLOCK_H */
