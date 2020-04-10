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
#include "linCartInterp.h"

namespace cart_interp
{
void check_out_of_domain(int* dims, int nw, int** ijk_stencil, std::vector<double>& ref_coord)
{
  // adjust reference cooridnates to collapse dimensions if out of domain
  for(int i=0;i<nw;i++)
  {
    if(ijk_stencil[0][3*i]>=dims[0])
    {
      ijk_stencil[0][3*i]=dims[0]-1;
      ref_coord[0]=-1;
    }
    else if(ijk_stencil[0][3*i]<0)
    {
      ijk_stencil[0][3*i]=0;
      ref_coord[0]=1;
    }
    if(ijk_stencil[0][3*i+1]>=dims[1])
    {
      ijk_stencil[0][3*i+1]=dims[1]-1;
      ref_coord[1]=-1;
    }
    else if(ijk_stencil[0][3*i+1]<0)
    {
      ijk_stencil[0][3*i+1]=0;
      ref_coord[1]=1;
    }
    if(ijk_stencil[0][3*i+2]>=dims[2])
    {
      ijk_stencil[0][3*i+2]=dims[2]-1;
      ref_coord[2]=-1;
    }
    else if(ijk_stencil[0][3*i+2]<0)
    {
      ijk_stencil[0][3*i+2]=0;
      ref_coord[2]=1;
    }
  }
}

void compute_linear_weights(const std::vector<double>& ref_coord, double** weights)
{
  double eta = ref_coord[0];
  double mu  = ref_coord[1];
  double xi  = ref_coord[2];

  // assign weights based on shape functions for a FE brick element
  weights[0][0] = (1-eta)*(1+mu)*(1-xi)/8;
  weights[0][1] = (1-eta)*(1-mu)*(1-xi)/8;
  weights[0][2] = (1-eta)*(1-mu)*(1+xi)/8;
  weights[0][3] = (1-eta)*(1+mu)*(1+xi)/8;
  weights[0][4] = (1+eta)*(1+mu)*(1-xi)/8;
  weights[0][5] = (1+eta)*(1-mu)*(1-xi)/8;
  weights[0][6] = (1+eta)*(1-mu)*(1+xi)/8;
  weights[0][7] = (1+eta)*(1+mu)*(1+xi)/8;
}

void compute_ref_coords(double* ref_ratio, std::vector<double>& ref_coord)
{
  ref_coord[0] = (ref_ratio[0]-0.5>=0) ? ((ref_ratio[0]-0.5)*2-1) : ((ref_ratio[0]+0.5)*2-1);
  ref_coord[1] = (ref_ratio[1]-0.5>=0) ? ((ref_ratio[1]-0.5)*2-1) : ((ref_ratio[1]+0.5)*2-1);
  ref_coord[2] = (ref_ratio[2]-0.5>=0) ? ((ref_ratio[2]-0.5)*2-1) : ((ref_ratio[2]+0.5)*2-1);
}

void create_linear_donor_stencil(int* ijk_cell, double* ref_ratio, int** ijk_stencil)
{
  // determine if donor stencil is to right or left of the cell-center
   ijk_stencil[0][0] = (ref_ratio[0]-0.5>=0) ? (ijk_cell[0]) : (ijk_cell[0]-1);

   // determine if donor stencil is front of or behind the cell-center
   ijk_stencil[0][1] = (ref_ratio[1]-0.5>=0) ? (ijk_cell[1]+1) : (ijk_cell[1]);

   // determine if donor stencil is above or below the cell-center
   ijk_stencil[0][2] = (ref_ratio[2]-0.5>=0) ? (ijk_cell[2]) : (ijk_cell[2]-1);

   // node 2 of reference brick element
   ijk_stencil[0][3]=ijk_stencil[0][0];
   ijk_stencil[0][4]=ijk_stencil[0][1]-1;
   ijk_stencil[0][5]=ijk_stencil[0][2];

   // node 3 of reference brick element
   ijk_stencil[0][6]=ijk_stencil[0][0];
   ijk_stencil[0][7]=ijk_stencil[0][1]-1;
   ijk_stencil[0][8]=ijk_stencil[0][2]+1;

   // node 4 of reference brick element
   ijk_stencil[0][9] =ijk_stencil[0][0];
   ijk_stencil[0][10]=ijk_stencil[0][1];
   ijk_stencil[0][11]=ijk_stencil[0][2]+1;

   // node 5 of reference brick element
   ijk_stencil[0][12]=ijk_stencil[0][0]+1;
   ijk_stencil[0][13]=ijk_stencil[0][1];
   ijk_stencil[0][14]=ijk_stencil[0][2];

   // node 6 of reference brick element
   ijk_stencil[0][15]=ijk_stencil[0][0]+1;
   ijk_stencil[0][16]=ijk_stencil[0][1]-1;
   ijk_stencil[0][17]=ijk_stencil[0][2];

   // node 7 of reference brick element
   ijk_stencil[0][18]=ijk_stencil[0][0]+1;
   ijk_stencil[0][19]=ijk_stencil[0][1]-1;
   ijk_stencil[0][20]=ijk_stencil[0][2]+1;

   // node 8 of reference brick element
   ijk_stencil[0][21]=ijk_stencil[0][0]+1;
   ijk_stencil[0][22]=ijk_stencil[0][1];
   ijk_stencil[0][23]=ijk_stencil[0][2]+1;
}

void linear_interpolation(int* pdegree, int* ijk_cell, int* dims, double* ref_ratio,
  int* nw, int** ijk_stencil, double** weights)
{
  if(*pdegree != 0)
    throw std::runtime_error("#tioga: linear Cartesian interpolation not setup for pdegree > 0");

  std::vector<double> ref_coord(3,0); // array of reference coordinates

  // 8-node donor stencil where the nodes are neighboring cell-centers
  // ordering of cell centers in donor stencil is such that
  // left->right is along x-axis. back->front is along y-axis. bottom->top is along z-axis
  create_linear_donor_stencil(ijk_cell,ref_ratio,ijk_stencil);
  compute_ref_coords(ref_ratio,ref_coord);
  check_out_of_domain(dims,*nw,ijk_stencil,ref_coord);
  compute_linear_weights(ref_coord,weights);
}
}
