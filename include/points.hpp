/*!
 * \file points.cpp
 * \brief Functions related to quadrature point locations and weights (header)
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
#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "error.hpp"

/*! Useful 3D point object with simple geometric functions */
struct point
{
  double x, y, z;

  point() {
    x = 0;
    y = 0;
    z = 0;
  }

  point (double _x, double _y, double _z) {
    x = _x;
    y = _y;
    z = _z;
  }

  point(const double* pt, int nDims=3) {
    x = pt[0];
    y = pt[1];
    if (nDims==3)
      z = pt[2];
    else
      z = 0;
  }

  inline
  void zero() {
    x = 0;
    y = 0;
    z = 0;
  }

  inline
  double& operator[](int ind) {
    switch(ind) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        std::cout << "ind = " << ind << ": " << std::flush;
        FatalError("Invalid index for point struct.");
    }
  }

  inline
  double operator[](int ind) const {
    switch(ind) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        std::cout << "ind = " << ind << ": " << std::flush;
        FatalError("Invalid index for point struct.");
    }
  }

  point operator=(double* a) {
    struct point pt;
    pt.x = a[0];
    pt.y = a[1];
    pt.z = a[2];
    return pt;
  }

  point operator-(point b) {
    struct point c;
    c.x = x - b.x;
    c.y = y - b.y;
    c.z = z - b.z;
    return c;
  }

  point operator+(point b) {
    struct point c;
    c.x = x + b.x;
    c.y = y + b.y;
    c.z = z + b.z;
    return c;
  }

  point operator/(point b) {
    struct point c;
    c.x = x / b.x;
    c.y = y / b.y;
    c.z = z / b.z;
    return c;
  }

  point& operator+=(point b) {
    x += b.x;
    y += b.y;
    z += b.z;
    return *this;
  }

  point& operator-=(point b) {
    x -= b.x;
    y -= b.y;
    z -= b.z;
    return *this;
  }

  point& operator+=(double* b) {
    x += b[0];
    y += b[1];
    z += b[2];
    return *this;
  }

  point& operator-=(double* b) {
    x -= b[0];
    y -= b[1];
    z -= b[2];
    return *this;
  }

  point& operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  point& operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  double operator*(point b) {
    return x*b.x + y*b.y + z*b.z;
  }

  void abs(void) {
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);
  }

  double norm(void) {
    return std::sqrt(x*x+y*y+z*z);
  }

  point cross(point b) {
    point v;
    v.x = y*b.z - z*b.y;
    v.y = z*b.x - x*b.z;
    v.z = x*b.y - y*b.x;
    return v;
  }

};

//! For clearer notation when a vector is implied, rather than a point
typedef struct point Vec3;

enum ETYPE
{
  QUAD, HEX
};

/*! Get the locations of the Consistent Grid Points for the CSC metrics
 *  (just an equi-spaced grid containing the corners of the element; see Abe et al 2016) */
std::vector<point> getLocCGPts(int order, int nDims);

//! Get the reference-domain location of the solution points for the given element & polynomial order
std::vector<point> getLocSpts(int eType, int order, std::string sptsType);

//! Get the reference-domain location of the flux points for the given element & polynomial order
std::vector<point> getLocFpts(int eType, int order, std::string sptsType);

//! Get the reference-domain location of the plot points for the given element & polynomial order
std::vector<point> getLocPpts(int eType, int order, std::string sptsType);

//! Get the point locations of the requested type (i.e. Gauss, Lobatto) for the given order
std::vector<double> getPts1D(std::string ptsType, int order);

//! Get the Gauss quadrature weights for the Gauss points of the given order [2D]
std::vector<double> getQptWeights(int order, int nDims);

//! Get the Gauss quadrature weights for the Gauss points of the given order [1D]
std::vector<double> getQptWeights1D(int order);

//! Get quadrature rule (points & weights) for a tetrahedron for a given order
void getQuadRuleTet(int order, std::vector<point> &locQpts, std::vector<double> &weights);

//! Get quadrature rule (points & weights) for a triangle for a given order
void getQuadRuleTri(int order, std::vector<point> &locQpts, std::vector<double> &weights);
