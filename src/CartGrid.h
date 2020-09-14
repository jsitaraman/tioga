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

#ifndef CARTGRID_H
#define CARTGRID_H

#include <cstdlib>

class CartGrid
{
 private:
  double xlosup[3];
  double *dxlvl{nullptr};
  int *lcount{nullptr};
  int maxlevel;

 public :
  int *global_id;
  int *level_num;
  int *proc_id;
  int *local_id;
  int *ilo;
  int *ihi;
  int *dims;
  int myid;
  int nf;
  double *xlo;
  double *dx;
  int ngrids;
  void (*donor_frac) (int *,double *,int *,double *) = nullptr;
   
  CartGrid() { ngrids=0;global_id=NULL;level_num=NULL;local_id=NULL;
    proc_id=NULL;local_id=NULL;ilo=NULL;ihi=NULL;dims=NULL;
               xlo=NULL;dx=NULL;dxlvl=NULL;lcount=NULL;donor_frac=nullptr;};
  ~CartGrid() { 
    if (global_id) free(global_id);
    if (level_num) free(level_num);
    if (proc_id) free(proc_id);
    if (local_id) free(local_id);
    if (ilo) free(ilo);
    if (ihi) free(ihi);
    if (dims) free(dims);
    if (xlo) free(xlo);
    if (dx) free(dx);
    if (lcount) free(lcount);
    if (dxlvl) free(dxlvl);
  };
  void registerData(int nf,int *idata,double *rdata,int ngridsin);
  void preprocess(void);     
  void search(double *x,int *donorid,int nsearch);
  void setcallback(void (*f1)(int *,double *,int *,double *)) 
  {   donor_frac=f1;  }
};

#endif /* CARTGRID_H */
