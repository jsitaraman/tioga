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
#include <tuple>

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

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/*! Evaluate the 1D Lagrange polynomial mode based on points x_lag at point y */
double Lagrange(std::vector<double> &x_lag, double y, uint mode);

/*! Evaluate the first derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double dLagrange(std::vector<double> &x_lag, double y, uint mode);

void adjoint_3x3(double *mat, double *adj);
void adjoint_4x4(double *mat, double *adj);

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

/*! Return bounding box of a collection of points after applying linear transform */
void getBoundingBox(double *pts, int nPts, int nDims, double *bbox, double *Smat);

/*! Get the centroid of a collection of points */
void getCentroid(double *pts, int nPts, int nDims, double *xc);

/*! Create an N-dimensional simplex of size L centered at x0 */
void getSimplex(int nDims, const std::vector<double> &x0, double L, std::vector<double> &X);

/*! Get reference location out_rst of point in_xyz within an element defined by the points in xv */
bool getRefLocNewton(double *xv, double *in_xyz, double *out_rst, int nNodes, int nDims);

/*! Compute the volume of a high-order quad or hex */
double computeVolume(double *xv, int nNodes, int nDims);

/*! Determine whether a given face and cell intersect */
Vec3 intersectionCheck2(double *fxv, int nfv, double *exv, int nev, int nDims);
Vec3 intersectionCheck(double *fxv, int nfv, double *exv, int nev, int nDims);

std::vector<int> get_int_list(int N, int start = 0);
std::vector<uint> get_int_list(uint N, uint start = 0);

std::vector<int> reverse_map(const std::vector<int> &map1);

//! Map a structured ijk-type index to the equivalent Gmsh node index
std::vector<int> structured_to_gmsh_quad(unsigned int nNodes);
std::vector<int> structured_to_gmsh_hex(unsigned int nNodes);

//! Map a Gmsh node index to the equivalent structured ijk-type index
std::vector<int> gmsh_to_structured_quad(unsigned int nNodes);
std::vector<int> gmsh_to_structured_hex(unsigned int nNodes);

double quick_select(int* inds, double* arr, int n);

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

/* ------------------------- Optimization Functions ------------------------ */

typedef struct NM_FVAL
{
  double f;
  double g;
  std::vector<double> x;
} NM_FVAL;

inline bool operator<  (const NM_FVAL &a, const NM_FVAL &b) { return a.f <  b.f; }
inline bool operator<= (const NM_FVAL &a, const NM_FVAL &b) { return a.f <= b.f; }
inline bool operator>  (const NM_FVAL &a, const NM_FVAL &b) { return a.f >  b.f; }
inline bool operator>= (const NM_FVAL &a, const NM_FVAL &b) { return a.f >= b.f; }

inline std::vector<double> operator+ (const std::vector<double> &a,
                               const std::vector<double> &b)
{
  std::vector<double> c(a);
  for (uint i = 0; i < c.size(); i++)
    c[i] += b[i];

  return c;
}

inline std::vector<double> operator* (double a, const std::vector<double> &b)
{
  std::vector<double> c(b);
  for (uint i = 0; i < c.size(); i++)
    c[i] *= a;

  return c;
}

inline std::vector<double> operator/ (const std::vector<double> &a, double b)
{
  std::vector<double> c(a);
  for (uint i = 0; i < c.size(); i++)
    c[i] /= b;

  return c;
}

inline std::vector<double> operator- (const std::vector<double> &a,
                               const std::vector<double> &b)
{
  std::vector<double> c(a);
  for (uint i = 0; i < c.size(); i++)
    c[i] -= b[i];

  return c;
}

inline unsigned getIndAbsMax(const std::vector<double> &vec)
{
  unsigned ind = 0;
  double maxVal = -1e15;
  for (unsigned i = 0; i < vec.size(); i++)
  {
    if (std::abs(vec[i]) > maxVal)
    {
      maxVal = std::abs(vec[i]);
      ind = i;
    }
  }

  return ind;
}

/*! ----------------------------------- Optimization Routines ----------------------------------- */

static std::vector<double> Xn;  // Point with the highest value of F
static std::vector<double> X0;  // Centroid of all other points
static std::vector<double> Xr;  // Reflected point
static std::vector<double> Xe;  // Expanded point
static std::vector<double> Xc;  // Contracted point

static std::vector<double> Dx;

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
NM_FVAL NelderMead(const std::vector<double> &U0, Func minFunc, double L = 1.)
{
  /// TODO: Optimize the crap out of this

  int nVars = U0.size();
  int nPts = nVars+1;
  std::vector<NM_FVAL> FX(nPts);

  // Starting location for search
  if (nVars == 5) // 3D intersection
  {
    double x1 = L*.5;
    double x2 = L*std::sqrt(3)/2.;
    double x3 = L*cos(M_PI/3.);
    double x4 = L*sin(M_PI/3.);
    FX[0].x = {L,    0,  x1,  x1,   0};
    FX[1].x = {-x1, x2,  x1,  x3,  x4};
    FX[2].x = {-x1,-x2,  x1,  x3, -x4};
    FX[3].x = {-L,   0, -x1, -x1,   0};
    FX[4].x = {x1, -x2, -x1, -x3, -x4};
    FX[5].x = {x1,  x2, -x1, -x3,  x4};
  }
  else
  {
    std::vector<double> X;
    getSimplex(nVars, U0, L, X);

    for (int i = 0; i < nPts; i++)
    {
      FX[i].x.resize(nVars);
      for (int j = 0; j < nVars; j++)
        FX[i].x[j] = X[i*nVars+j];
    }
  }

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++)
    FX[i].f = minFunc(FX[i].x);

  std::sort(FX.begin(),FX.end());

  Xn.resize(nVars);  // Point with the highest value of F
  X0.resize(nVars);  // Centroid of all other points
  Xr.resize(nVars);  // Reflected point
  Xe.resize(nVars);  // Expanded point
  Xc.resize(nVars);  // Contracted point

  // Use a relative tolerance...?
  double tol = 1e-8;
  int iter = 0;
  while (iter < 200 && FX[0].f > tol) {
    Xn = FX[nVars].x;

    // Take centroid of all points besides Xn
    X0.assign(nVars,0);
    for (int j=0; j<nPts-1; j++)
      X0 = X0 + FX[j].x / (nPts-1.);

    // Reflect Xn around X0
    Xr = X0 + X0-Xn;

    double Fr = minFunc(Xr);

    // Determine what to do with the new point
    if (Fr < FX[nPts-2].f) {
      // We will be keeping this point
      if (Fr < FX[0].f) {
        // This one's good; keep going! Expand from Xr
        Xe = Xr + X0-Xn;

        double Fe = minFunc(Xe);

        if (Fe < Fr) {
          // This one's even better; use it instead
          FX[nPts-1].f = Fe;
          FX[nPts-1].x = Xe;
        }
        else {
          // Xe/Fe was no better; stick with Fr, Xr
          FX[nPts-1].f = Fr;
          FX[nPts-1].x = Xr;
        }
      }
      else {
        // This one's somewhere in the middle; replace Xn with Xr
        FX[nPts-1].f = Fr;
        FX[nPts-1].x = Xr;
      }
    }
    else {
      // Try reducing the size of the simplex
      Xc = X0 + X0-Xn;

      double Fc = minFunc(Xc);
      if (Fc < FX[nPts-1].f) {
        // Bringing this point in is better; use it
        FX[nPts-1].f = Fc;
        FX[nPts-1].x = Xc;
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        Xc = FX[0].x;
        for (int i=1; i<nPts; i++) {
          FX[i].x = Xc + .5*(FX[i].x-Xc);
          FX[i].f = minFunc(FX[i].x);
        }
      }
    }

    std::sort(FX.begin(),FX.end());

    // Continue to iterate
    iter++;
  }

  return FX[0];
}

/*!
 * Nelder-Mead Minimzation Routine With Constraints
 *
 * Returns: the coordinates of the point found to be the minimum
 *
 * \param[in] U0: Starting coordinates for search
 * \param[in] minFunc: a normal or lambda function accepting a
 *            std::vector<double> and returning a double
 * \param[in] G: Constraint function of form G(x) < 0
 */
template<int nVars, typename MinFunc, typename Constraint>
NM_FVAL NelderMead_constrained(const std::vector<double> &U0,
    MinFunc minFunc, Constraint G, double L = 1.)
{
  /// TODO: Optimize the crap out of this

//  int nVars = U0.size();
  int nPts = nVars+1;
  std::vector<NM_FVAL> FX(nPts);

//  std::ofstream fout("nm-output.csv");

  // Starting location for search
  if (nVars == 5)
  {
    double x1 = L*.6;
    double x2 = L*std::sqrt(3.)/2.;
    double x3 = x1*cos(M_PI/3.);
    double x4 = x1*sin(M_PI/3.);
    FX[0].x = {L,    0,  x1,  x1,   0};
    FX[1].x = {-x1, x2,  x1,  x3,  x4};
    FX[2].x = {-x1,-x2,  x1,  x3, -x4};
    FX[3].x = {-L,   0, -x1, -x1,   0};
    FX[4].x = {x1, -x2, -x1, -x3, -x4};
    FX[5].x = {x1,  x2, -x1, -x3,  x4};
  }
  else
  {
    std::vector<double> X;
    getSimplex(nVars, U0, L, X);

    for (int i = 0; i < nPts; i++)
    {
      FX[i].x.resize(nVars);
      for (int j = 0; j < nVars; j++)
        FX[i].x[j] = X[i*nVars+j];
    }
  }

  // Evaluate the 'function' at the initial 'points'
  double maxG = -1;
  for (int i=0; i<nPts; i++)
  {
    FX[i].f = minFunc(FX[i].x);
    FX[i].g = G(FX[i].x);
    maxG = std::max(FX[i].g, maxG);
  }

  if (nVars == 5)
  {
    // Sort the repsective face & cell pts so that the closest are paired, then
    // the next closest, etc.
    std::vector<std::vector<double>> XI(nPts), XJ(nPts);
    for (int i = 0; i < nPts; i++)
    {
      XI[i] = {FX[i].x[0], FX[i].x[1]};
      XJ[i] = {FX[i].x[2], FX[i].x[3], FX[i].x[4]};
    }

    for (int k = 0; k < nPts; k++)
    {
      double minVal = 1e15;
      int minIndI = -1;
      int minIndJ = -1;
      for (int i = 0; i < nPts-k; i++)
      {
        for (int j = 0; j < nPts-k; j++)
        {
          std::vector<double> xk = {XI[i][0], XI[i][1], XJ[j][0], XJ[j][1], XJ[j][2]};
          double fk = minFunc(xk);
          if (fk < minVal)
          {
            minVal = fk;
            minIndJ = j;
            minIndI = i;
          }
        }
      }

      for (int m = 0; m < 2; m++)
        FX[k].x[m] = XI[minIndI][m];
      for (int m = 0; m < 3; m++)
        FX[k].x[m+2] = XJ[minIndJ][m];
      FX[k].f = minVal;
      FX[k].g = G(FX[k].x);

      XI.erase(XI.begin()+minIndI);
      XJ.erase(XJ.begin()+minIndJ);
    }
  }

  std::sort(FX.begin(),FX.end());

  Xn.resize(nVars);  // Point with the highest value of F
  X0.resize(nVars);  // Centroid of all other points
  Xr.resize(nVars);  // Reflected point
  Xe.resize(nVars);  // Expanded point
  Xc.resize(nVars);  // Contracted point

  Dx.resize(nVars);

  // Use a relative tolerance...?
  double tol = 1e-8;
  int iter = 0;
  while (iter < 200 && FX[0].f > tol)
  {
    Xn = FX[nVars].x;

    // Take centroid of all points besides Xn
    X0.assign(nVars,0);
    for (int j = 0; j < nPts-1; j++)
      X0 = X0 + FX[j].x / (nPts-1.);

    // Reflect Xn around X0
    Dx = X0 - Xn;
    Xr = X0 + Dx;

    double Gr = G(Xr);
    if (Gr > 0)
    {
      int ind = getIndAbsMax(Xr);
      double fac = .8*(1.-(std::abs(Xr[ind])-1.)/(std::abs(Dx[ind])));
      Xr = X0 + fac*Dx;
    }

    double Fr = minFunc(Xr);
    Gr = G(Xr);

    // Determine what to do with the new point
    if (Fr < FX[nPts-2].f and Gr < 0)
    {
      // We will be keeping this point
      if (Fr < FX[0].f)
      {
        // This one's good; keep going! Expand from Xr
        Xe = Xr + X0-Xn;

        double Fe = minFunc(Xe);
        double Ge = G(Xe);

        if (Fe < Fr and Ge < 0)
        {
          // This one's even better; use it instead
          FX[nPts-1].f = Fe;
          FX[nPts-1].x = Xe;
        }
        else
        {
          // Xe/Fe was no better; stick with Fr, Xr
          FX[nPts-1].f = Fr;
          FX[nPts-1].x = Xr;
        }
      }
      else
      {
        // This one's somewhere in the middle; replace Xn with Xr
        FX[nPts-1].f = Fr;
        FX[nPts-1].x = Xr;
      }
    }
    else
    {
      // Try reducing the size of the simplex
      Xc = X0 + .5*(X0-Xn);

      double Fc = minFunc(Xc);
      double Gc = G(Xc);
      if (Fc < FX[nPts-1].f and Gc < 0)
      {
        // Bringing this point in is better; use it
        FX[nPts-1].f = Fc;
        FX[nPts-1].x = Xc;
      }
      else
      {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        Xc = FX[0].x;
        for (int i = 1; i < nPts; i++)
        {
          auto xn = FX[i].x;
          xn = Xc + .5*(FX[i].x - Xc);

          double Gn = G(xn);
          if (Gn < 0)
          {
            FX[i].x = xn;
            FX[i].f = minFunc(FX[i].x);
            FX[i].g = G(FX[i].x);
          }
        }
      }
    }

    std::sort(FX.begin(),FX.end());

//    for (auto &fx : FX)
//    {
//      fout << iter << ", " << fx.f;
//      for (auto &x : fx.x)
//        fout << ", " << x;
//      fout << std::endl;
//    }
    // Continue to iterate
    iter++;
  }

  return FX[0];
}


template<typename MinFunc, typename Constraint>
double NelderMeadStep_constrained(std::vector<NM_FVAL> &FX, MinFunc minFunc,
    Constraint G)
{
  /// TODO: Optimize the crap out of this

  int nPts = FX.size();
  int nVars = nPts-1;

  std::sort(FX.begin(),FX.end());

  auto Xn = FX[nVars].x;

  // Take centroid of all points besides Xn
  std::vector<double> X0(nVars);
  for (int j = 0; j < nPts-1; j++)
    X0 = X0 + FX[j].x / (nPts-1.);

  // Reflect Xn around X0
  auto Dx = X0 - Xn;
  auto Xr = X0 + Dx;

  double Gr = G(Xr);
  if (Gr > 0)
  {
    int ind = getIndAbsMax(Xr);
    double fac = .8*(1.-(std::abs(Xr[ind])-1.)/(std::abs(Dx[ind])));
    Xr = X0 + fac*Dx;
  }

  double Fr = minFunc(Xr);
  Gr = G(Xr);

  // Determine what to do with the new point
  if (Fr < FX[nPts-2].f and Gr < 0) {
    // We will be keeping this point
    if (Fr < FX[0].f) {
      // This one's good; keep going! Expand from Xr
      auto Xe = Xr + Dx;

      double Fe = minFunc(Xe);
      double Ge = G(Xe);

      if (Fe < Fr and Ge < 0) {
        // This one's even better; use it instead
        FX[nPts-1].f = Fe;
        FX[nPts-1].x = Xe;
      } else {
        // Xe/Fe was no better; stick with Fr, Xr
        FX[nPts-1].f = Fr;
        FX[nPts-1].x = Xr;
      }
    } else {
      // This one's somewhere in the middle; replace Xn with Xr
      FX[nPts-1].f = Fr;
      FX[nPts-1].x = Xr;
    }
  } else {
    // Try reducing the size of the simplex
    auto Xc = X0 + .5*Dx;

    double Fc = minFunc(Xc);
    double Gc = G(Xc);

    if (Fc < FX[nPts-1].f and Gc < 0) {
      // Bringing this point in is better; use it
      FX[nPts-1].f = Fc;
      FX[nPts-1].x = Xc;
    } else {
      // Bringing this point in didn't work; shrink the simplex onto
      // the smallest-valued vertex (while adhering to constraints)
      Xc = FX[0].x;
      for (int i = 1; i < nPts; i++) {
        auto xn = Xc + .5*(FX[i].x - Xc);

        double Gn = G(xn);
        if (Gn < 0) {
          FX[i].x = xn;
          FX[i].f = minFunc(FX[i].x);
        }
      }
    }
  }

  std::sort(FX.begin(),FX.end());

  return FX[0].f;
}

} // namespace tg_funcs
#endif // FUNCS_HPP

