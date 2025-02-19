// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#ifndef CARTGRID_H
#define CARTGRID_H

#include <cstdlib>

class CartGrid
{
 private:
  double xlosup[3];
  double *dxlvl;
  int *lcount;
  int maxlevel;

 public :
  int *global_id;
  int *level_num;
  int *proc_id;
  int *local_id;
  int *porder;
  int *ilo;
  int *ihi;
  int *dims;
  int myid;
  int nf,qstride;
  double *xlo;
  double *dx;
  double *qnode;
  int ngrids;
  void (*donor_frac) (int *,double *,int *,double *);
   
  CartGrid() { ngrids=0;global_id=NULL;level_num=NULL;local_id=NULL;porder=NULL;
    proc_id=NULL;local_id=NULL;ilo=NULL;ihi=NULL;dims=NULL;
               xlo=NULL;dx=NULL;dxlvl=NULL;lcount=NULL;qnode=NULL;donor_frac=nullptr;};
  ~CartGrid() { 
    if (global_id) free(global_id);
    if (level_num) free(level_num);
    if (proc_id) free(proc_id);
    if (local_id) free(local_id);
    if (porder) free(porder);
    if (ilo) free(ilo);
    if (ihi) free(ihi);
    if (dims) free(dims);
    if (xlo) free(xlo);
    if (dx) free(dx);
    if (lcount) free(lcount);
    if (dxlvl) free(dxlvl);
    if (qnode) free(qnode);
  };
  void registerData(int nf,int qstride,double *qnodein,
		    int *idata,double *rdata,
		    int ngridsin,int qnodesize);
  void preprocess(void);     
  void search(double *x,int *donorid,int nsearch);
  void setcallback(void (*f1)(int *,double *,int *,double *)) 
  {   donor_frac=f1;  }
};

#endif /* CARTGRID_H */
