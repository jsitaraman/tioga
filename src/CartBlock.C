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
#include "codetypes.h"
#include "CartBlock.h"
#include "CartGrid.h"
void CartBlock::getInterpolatedData(int *nints,int *nreals,int **intData,
				    double **realData,
				    double *q,
				    int nvar, int interptype)
{
}


void CartBlock::update(double *qval, int index,int nq)
{
  int i;
  for(i=0;i<nq;i++)
    q[index+gsize3*i]=qval[i];
}

  
void CartBlock::preprocess(CartGrid *cg)
  {
    int nf;    
    for(int n=0;n<3;n++) xlo[n]=cg->xlo[3*global_id+n];
    for(int n=0;n<3;n++) dx[n]=cg->dx[3*global_id+n];
    jd=cg->ihi[3*global_id]  -cg->ilo[3*global_id  ]+1;
    kd=cg->ihi[3*global_id+1]-cg->ilo[3*global_id+1]+1;
    ld=cg->ihi[3*global_id+2]-cg->ilo[3*global_id+2]+1;
    nf=cg->nf;
    gsize1=(jd+2*nf);
    gsize2=gsize1*(kd+2*nf);
    gsize3=gsize2*(ld+2*nf);    
  };
