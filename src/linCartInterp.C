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
#include <stdexcept>

namespace cart_interp
{
void compute_1d_bases(const std::vector<double>& ref_coord,
  std::vector<double>& phi_x, std::vector<double>& phi_y, std::vector<double>& phi_z)
{
  double zeta = ref_coord[0];
  double mu   = ref_coord[1];
  double eta  = ref_coord[2];

  phi_x[0]=(1-zeta)/2;
  phi_x[1]=(1+zeta)/2;

  phi_y[0]=(1-mu)/2;
  phi_y[1]=(1+mu)/2;

  phi_z[0]=(1-eta)/2;
  phi_z[1]=(1+eta)/2;
}

void compute_weights(const std::vector<double>& ref_coord, double* weights)
{
  static constexpr int p = 1;
  std::vector<double> phi_x(p+1,0);
  std::vector<double> phi_y(p+1,0);
  std::vector<double> phi_z(p+1,0);

  compute_1d_bases(ref_coord,phi_x,phi_y,phi_z);

  int ind = 0;
  for(int k=0;k<p+1;k++) {
    for(int j=0;j<p+1;j++) {
      for(int i=0;i<p+1;i++) {
        weights[ind++] = phi_z[k]*phi_y[j]*phi_x[i];
      }
    }
  }
}

void compute_ref_coords_cell(double* ref_ratio, std::vector<double>& ref_coord)
{
  static constexpr int p = 1;
  // reference coordinates in -1 to 1. Assumes cells of uniform sizes
  ref_coord[0] = (ref_ratio[0]-0.5>=0) ?
      ((ref_ratio[0]-0.5)/p*2-1) : ((ref_ratio[0]+0.5*(2*p-1))/p*2-1);
  ref_coord[1] = (ref_ratio[1]-0.5>=0) ?
      ((ref_ratio[1]-0.5)/p*2-1) : ((ref_ratio[1]+0.5*(2*p-1))/p*2-1);
  ref_coord[2] = (ref_ratio[2]-0.5>=0) ?
      ((ref_ratio[2]-0.5)/p*2-1) : ((ref_ratio[2]+0.5*(2*p-1))/p*2-1);
}

void compute_ref_coords_node(double* ref_ratio, std::vector<double>& ref_coord)
{
  static constexpr int p = 1;
  // reference coordinates in -1 to 1. Assumes cells of uniform sizes
  ref_coord[0] = ref_ratio[0]/p*2-1;
  ref_coord[1] = ref_ratio[1]/p*2-1;
  ref_coord[2] = ref_ratio[2]/p*2-1;
}

void create_donor_stencil(const int nf, int* ijk_cell, int* dims, double* ref_ratio, int* ijk_stencil, bool isNodal)
{
  static constexpr int basis = 2;

  int start_x, start_y, start_z = 0;
  int nX, nY, nZ = 0;
  if(isNodal) {
    start_x = ijk_cell[0];
    start_y = ijk_cell[1];
    start_z = ijk_cell[2];
    nX = dims[0]+1;
    nY = dims[1]+1;
    nZ = dims[2]+1;
  }
  else {
    // determine start node if donor stencil is
    // right/left or front/behind or above/below of ijk_cell cell-center
     start_x = (ref_ratio[0]-0.5>=0) ? (ijk_cell[0]) : (ijk_cell[0]-1);
     start_y = (ref_ratio[1]-0.5>=0) ? (ijk_cell[1]) : (ijk_cell[1]-1);
     start_z = (ref_ratio[2]-0.5>=0) ? (ijk_cell[2]) : (ijk_cell[2]-1);
     nX = dims[0];
     nY = dims[1];
     nZ = dims[2];
  }

   int ind = 0;
   for(int k=0;k<basis;k++) {
     for(int j=0;j<basis;j++) {
       for(int i=0;i<basis;i++) {
         ijk_stencil[ind++] = std::max(-nf, std::min(start_x+i*1, nX+nf-1));
         ijk_stencil[ind++] = std::max(-nf, std::min(start_y+j*1, nY+nf-1));
         ijk_stencil[ind++] = std::max(-nf, std::min(start_z+k*1, nZ+nf-1));
       }
     }
   }
}

void linear_interpolation(const int nf, int* ijk_cell, int* dims, double* ref_ratio,
  int* nw, int* ijk_stencil, double* weights, bool isNodal)
{
  std::vector<double> ref_coord(3,0); // array of reference coordinates

  // 8-node donor stencil where the nodes are neighboring cell-centers
  // ordering of cell centers in donor stencil is such that
  // left->right is along x-axis. back->front is along y-axis. bottom->top is along z-axis
  create_donor_stencil(nf,ijk_cell,dims,ref_ratio,ijk_stencil,isNodal);

  if(isNodal)
    compute_ref_coords_node(ref_ratio,ref_coord);
  else
    compute_ref_coords_cell(ref_ratio,ref_coord);

  compute_weights(ref_coord,weights);
}
}
