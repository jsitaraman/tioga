// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "linCartInterp.h"
#include <stdexcept>

namespace cart_interp
{
void compute_1d_bases(const int& p, const std::vector<double>& ref_coord,
  std::vector<double>& phi_x, std::vector<double>& phi_y, std::vector<double>& phi_z)
{
  double zeta = ref_coord[0];
  double mu   = ref_coord[1];
  double eta  = ref_coord[2];

  switch (p) {
      case 1: // FE shape functions in a (-1,1) reference element
        phi_x[0]=(1-zeta)/2;
        phi_x[1]=(1+zeta)/2;

        phi_y[0]=(1-mu)/2;
        phi_y[1]=(1+mu)/2;

        phi_z[0]=(1-eta)/2;
        phi_z[1]=(1+eta)/2;
      break;

      default:
        throw std::runtime_error("#tioga: Cartesian 1d bases not setup for chosen interpolation order");
  }
}

void compute_weights(const int& p, const std::vector<double>& ref_coord, double* weights)
{
  std::vector<double> phi_x(p+1,0);
  std::vector<double> phi_y(p+1,0);
  std::vector<double> phi_z(p+1,0);

  compute_1d_bases(p,ref_coord,phi_x,phi_y,phi_z);

  int ind = 0;
  for(int k=0;k<p+1;k++) {
    for(int j=0;j<p+1;j++) {
      for(int i=0;i<p+1;i++) {
        weights[ind++] = phi_z[k]*phi_y[j]*phi_x[i];
      }
    }
  }
}

void compute_ref_coords(const int& p, double* ref_ratio, std::vector<double>& ref_coord)
{
  // reference coordinates in -1 to 1. Assumes cells of uniform sizes
  ref_coord[0] = (ref_ratio[0]-0.5>=0) ?
      ((ref_ratio[0]-0.5)/p*2-1) : ((ref_ratio[0]+0.5*(2*p-1))/p*2-1);
  ref_coord[1] = (ref_ratio[1]-0.5>=0) ?
      ((ref_ratio[1]-0.5)/p*2-1) : ((ref_ratio[1]+0.5*(2*p-1))/p*2-1);
  ref_coord[2] = (ref_ratio[2]-0.5>=0) ?
      ((ref_ratio[2]-0.5)/p*2-1) : ((ref_ratio[2]+0.5*(2*p-1))/p*2-1);
}

void create_donor_stencil(const int& p, int* ijk_cell, int* dims, double* ref_ratio, int* ijk_stencil)
{
  // determine start node if donor stencil is
  // right/left or front/behind or above/below of ijk_cell cell-center
   int start_x = (ref_ratio[0]-0.5>=0) ? (ijk_cell[0]) : (ijk_cell[0]-1);
   int start_y = (ref_ratio[1]-0.5>=0) ? (ijk_cell[1]) : (ijk_cell[1]-1);
   int start_z = (ref_ratio[2]-0.5>=0) ? (ijk_cell[2]) : (ijk_cell[2]-1);

   int ind = 0;
   for(int k=0;k<p+1;k++) {
     for(int j=0;j<p+1;j++) {
       for(int i=0;i<p+1;i++) {
         ijk_stencil[ind++] = std::max(0, std::min(start_x+i*1, dims[0]-1));
         ijk_stencil[ind++] = std::max(0, std::min(start_y+j*1, dims[1]-1));
         ijk_stencil[ind++] = std::max(0, std::min(start_z+k*1, dims[2]-1));
       }
     }
   }
}

void linear_interpolation(const int& p, int* ijk_cell, int* dims, double* ref_ratio,
  int* nw, int* ijk_stencil, double* weights)
{
  std::vector<double> ref_coord(3,0); // array of reference coordinates

  // 8-node donor stencil where the nodes are neighboring cell-centers
  // ordering of cell centers in donor stencil is such that
  // left->right is along x-axis. back->front is along y-axis. bottom->top is along z-axis
  create_donor_stencil(p,ijk_cell,dims,ref_ratio,ijk_stencil);
  compute_ref_coords(p,ref_ratio,ref_coord);
  compute_weights(p,ref_coord,weights);
}
}
