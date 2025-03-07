// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

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
