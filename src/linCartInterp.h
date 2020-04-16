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
#include <vector>

#ifndef LINCARTINTERP_H
#define LINCARTINTERP_H

namespace cart_interp
{
void compute_1d_bases(const int& p, const std::vector<double>& ref_coord,
  std::vector<double>& phi_x, std::vector<double>& phi_y, std::vector<double>& phi_z);

void compute_linear_weights(const std::vector<double>& ref_coord, double* weights);

void compute_ref_coords(double* ref_ratio, std::vector<double>& ref_coord);

void create_donor_stencil(const int& p, int* ijk_cell, int* dims, double* ref_ratio, int* ijk_stencil);

void linear_interpolation(const int& p, int* ijk_cell, int* dims, double* ref_ratio,
  int* nw, int* ijk_stencil, double* weights);
} // namespacee cart_interp

#endif /* LINCARTINTERP_H */
