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
//#include "codetypes.h"
class CartBlock
{
 private:
  int global_id;
  int jd,kd,ld,nf;
  int gsize1,gsize2,gsize3;
  int *ibl;
  double *q;
  double xlo[3]; 
  int ndonors;
  int interpListSize;
  INTERPLIST *interpList;

 public:
  CartBlock() { global_id=0;jd=kd=ld=0;ibl=NULL;q=NULL;interpListSize=0;interpList=NULL;};
  void registerData(int *idata,int *iblankin,double *qin,double *xloin,int nf)
  {
    global_id=idata[0];
    jd=idata[1];
    kd=idata[2];
    ld=idata[3];
    ibl=iblankin;
    q=qin;
    for(int n=0;n<3;n++) xlo[n]=xloin[n];
    gsize1=(jd+2*nf);
    gsize2=gsize1*(kd+2*nf);
    gsize3=gsize2*(ld+2*nf);
  };

  void getInterpolatedData(int *nints,int *nreals,int **intData,
			   double **realData,
			   double *q,
			   int nvar, int interptype);
  void update(double *qval,int index,int nq);
};
