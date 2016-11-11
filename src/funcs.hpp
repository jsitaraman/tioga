/*!
 * \file funcs.hpp
 * \brief Miscellaneous helper functions (header)
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * This file is part of the Tioga software library
 *
 * Tioga  is a tool for overset grid assembly on parallel distributed systems
 * Copyright (C) 2015-2016 Jay Sitaraman
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "error.hpp"
#include "points.hpp"

#define ThrowException(msg) \
{ std::stringstream s; s << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": " << msg; \
  throw std::runtime_error(s.str());}\

namespace tg_funcs
{

/* ---------------------------- Helpful Objects ---------------------------- */

point operator/(point a, double b);
point operator*(point a, double b);

bool operator<(const point &a, const point &b); // Just a sort of 'dummy' function for sorting purposes

std::ostream& operator<<(std::ostream &os, const point &pt);

std::ostream& operator<<(std::ostream &os, const std::vector<int> &vec);
std::ostream& operator<<(std::ostream &os, const std::vector<double> &vec);

double getDist(point a, point b);

/*! Evaluate the 1D Lagrange polynomial mode based on points x_lag at point y */
double Lagrange(std::vector<double> &x_lag, double y, uint mode);

/*! Evaluate the first derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double dLagrange(std::vector<double> &x_lag, double y, uint mode);

/* ---------------------------- Misc. Functions ---------------------------- */

/*! Calculate the adjoint of a 'size x size' matrix stored row-major in 'mat' */
std::vector<double> adjoint(const std::vector<double> &mat, unsigned int size);

/*! In-place adjoint calculation */
void adjoint(const std::vector<double> &mat, std::vector<double> &adj, unsigned int size);

/*! Calculate the determinant of a 'size x size' matrix stored row-major in 'mat' */
double determinant(const std::vector<double> &mat, unsigned int size);

void getBoundingBox(double *pts, int nPts, int nDims, double *bbox);

/*! Get reference location out_rst of point in_xyz within an element defined by the points in xv */
bool getRefLocNewton(double *xv, double *in_xyz, double *out_rst, int nNodes, int nDims);

/*! Compute the volume of a high-order quad or hex */
double computeVolume(double *xv, int nNodes, int nDims);

std::vector<int> get_int_list(int N, int start = 0);
std::vector<uint> get_int_list(uint N, uint start = 0);

std::vector<int> reverse_map(const std::vector<int> &map1);

//! Map a structured ijk-type index to the equivalent Gmsh node index
std::vector<int> structured_to_gmsh_quad(unsigned int nNodes);
std::vector<int> structured_to_gmsh_hex(unsigned int nNodes);

//! Map a Gmsh node index to the equivalent structured ijk-type index
std::vector<int> gmsh_to_structured_quad(unsigned int nNodes);
std::vector<int> gmsh_to_structured_hex(unsigned int nNodes);

template<typename T>
int findFirst(const std::vector<T>& vec, T val)
{
  for (int i = 0; i < vec.size(); i++)
    if (vec[i] == val)
      return i;

  return -1;
}

/* ---------------------------- Shape Functions ---------------------------- */

//! Shape function for linear or quadratic quad (TODO: Generalize to N-noded quad)
void shape_quad(const point &in_rs, std::vector<double> &out_shape, int nNodes);
void shape_quad(const point &in_rs, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic quad
void dshape_quad(const std::vector<point> loc_pts, double* out_dshape, int nNodes);
void dshape_quad(const point &in_rs, double* out_dshape, int nNodes);

//! Shape function for linear or quadratic hexahedron
void shape_hex(const point &in_rst, std::vector<double> &out_shape, int nNodes);
void shape_hex(const point &in_rst, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic hexahedron
void dshape_hex(const std::vector<point>& loc_pts, double* out_dshape, int nNodes);
void dshape_hex(const point &in_rst, double* out_dshape, int nNodes);

} // namespace tg_funcs

#endif // FUNCS_HPP

