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

#include <algorithm>
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

/* ------------------------ Mathematical Functions ------------------------- */

/*! Evaluate the 1D Lagrange polynomial mode based on points x_lag at point y */
double Lagrange(std::vector<double> &x_lag, double y, uint mode);

/*! Evaluate the first derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double dLagrange(std::vector<double> &x_lag, double y, uint mode);

/*! Calculate the adjoint of a 'size x size' matrix stored row-major in 'mat' */
std::vector<double> adjoint(const std::vector<double> &mat, unsigned int size);

/*! In-place adjoint calculation */
void adjoint(const std::vector<double> &mat, std::vector<double> &adj, unsigned int size);

/*! Calculate the determinant of a 'size x size' matrix stored row-major in 'mat' */
double determinant(const std::vector<double> &mat, unsigned int size);

Vec3 faceNormal(double* xv, int nDims);

/* ---------------------------- Misc. Functions ---------------------------- */

/*! Return the bounding box of a collection of points [min x,y,z; max x,y,z] */
void getBoundingBox(double *pts, int nPts, int nDims, double *bbox);

/*! Create an N-dimensional simplex of size L centered at x0 */
void getSimplex(int nDims, const std::vector<double> &x0, double L, std::vector<double> &X);

/*! Get reference location out_rst of point in_xyz within an element defined by the points in xv */
bool getRefLocNewton(double *xv, double *in_xyz, double *out_rst, int nNodes, int nDims);

/*! Compute the volume of a high-order quad or hex */
double computeVolume(double *xv, int nNodes, int nDims);

/*! Determine whether a given face and cell intersect */
Vec3 intersectionCheck(double *fxv, int nvf, double *exv, int nve, int nDims);

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

void shape_line(double xi, std::vector<double> &out_shape, int nNodes);
void shape_line(double xi, double* out_shape, int nNodes);

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

/*!
 * Nelder-Mead Minimzation Routine
 *
 * Returns: the coordinates of the point found to be the minimum
 *
 * \param[in] U0: Starting coordinates for search
 * \param[in] minFunc: a normal or lambda function accepting a
 *            std::vector<double> and returning a double
 */
template<typename Func>
std::pair<double,std::vector<double>> NelderMead(const std::vector<double> &U0, Func minFunc, double L = 1.)
{
  /// TODO: Optimize the crap out of this

  int nVars = U0.size();
  int nPts = nVars+1;
  std::vector<std::pair<double,std::vector<double>>> FX(nPts);

  // Starting location for search
  std::vector<double> X;
  getSimplex(nVars, U0, L, X);

  for (int i = 0; i < nPts; i++)
    for (int j = 0; j < nVars; j++)
      FX[i].second[j] = X[i*nVars+j];

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++)
    FX[i].first = minFunc(FX[i].second);

  std::sort(FX.begin(),FX.end());

  std::vector<double> Xn(nVars);  // Point with the highest value of F
  std::vector<double> X0(nVars);  // Centroid of all other points
  std::vector<double> Xr(nVars);  // Reflected point
  std::vector<double> Xe(nVars);  // Expanded point
  std::vector<double> Xc(nVars);  // Contracted point

  // Use a relative tolerance...?
  double tol = 1e-8;
  int iter = 0;
  while (iter < 200 && FX[0].first > tol) {
    Xn = FX[nVars].second;

    // Take centroid of all points besides Xn
    for (int j=0; j<nPts-1; j++)
      for (int k=0; k<nVars; k++)
        X0[k] += FX[j].second[k]/(nPts-1);

    // Reflect Xn around X0
    for (int k=0; k<nVars; k++)
      Xr[k] = X0[k] + (X0[k]-Xn[k]);

    double Fr = minFunc(Xr);

    // Determine what to do with the new point
    if (Fr < FX[nPts-2].first) {
      // We will be keeping this point
      if (Fr < FX[0].first) {
        // This one's good; keep going! Expand from Xr
        for (int i=0; i<nVars; i++)
          Xe[i] = Xr[i] + (X0[i]-Xn[i]);
        double Fe = minFunc(Xe);

        if (Fe < Fr) {
          // This one's even better; use it instead
          FX[nPts-1].first = Fe;
          FX[nPts-1].second = Xe;
        }
        else {
          // Xe/Fe was no better; stick with Fr, Xr
          FX[nPts-1].first = Fr;
          FX[nPts-1].second = Xr;
        }
      }
      else {
        // This one's somewhere in the middle; replace Xn with Xr
        FX[nPts-1].first = Fr;
        FX[nPts-1].second = Xr;
      }
    }
    else {
      // Try reducing the size of the simplex
      for (int i=0; i<nVars; i++)
        Xc[i] = X0[i] - (X0[i]-Xn[i])*.5;
      double Fc = minFunc(Xc);
      if (Fc < FX[nPts-1].first) {
        // Bringing this point in is better; use it
        FX[nPts-1].first = Fc;
        FX[nPts-1].second = Xc;
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        Xc = FX[0].second;
        for (int i=1; i<nPts; i++) {
          for (int j=0; j<nVars; j++) {
            FX[i].second[j] = Xc[j] + 0.5*(FX[i].second[j]-Xc[j]);
          }
          FX[i].first = minFunc(FX[i].second);
        }
      }
    }

    std::sort(FX.begin(),FX.end());

    // Continue to iterate
    iter++;
  }

  return FX[0];
}

} // namespace tg_funcs
#endif // FUNCS_HPP

