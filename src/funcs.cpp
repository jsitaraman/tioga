/*!
 * \file funcs.cpp
 * \brief Miscellaneous helper functions
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
#include <cmath>
#include <limits.h>
#include <omp.h>

#include "funcs.hpp"

namespace tg_funcs
{

point operator/(point a, double b) { return a/=b; }
point operator*(point a, double b) { return a*=b; }

bool operator<(const point &a, const point &b) { return a.x < b.x; }

std::ostream& operator<<(std::ostream &os, const point &pt)
{
  os << "(x,y,z) = " << pt.x << ", " << pt.y << ", " << pt.z;
  return os;
}

std::ostream& operator<<(std::ostream &os, const std::vector<int> &vec)
{
  for (auto &val:vec) std::cout << val << ", ";
  return os;
}

std::ostream& operator<<(std::ostream &os, const std::vector<double> &vec)
{
  for (auto &val:vec) std::cout << val << ", ";
  return os;
}

double Lagrange(std::vector<double> &x_lag, double y, uint mode)
{
  double lag = 1.0;

  for (uint i = 0; i < x_lag.size(); i++) {
    if (i != mode) {
      lag = lag*((y-x_lag[i])/(x_lag[mode]-x_lag[i]));
    }
  }

  return lag;
}

double dLagrange(std::vector<double> &x_lag, double y, uint mode)
{
  uint i, j;
  double dLag, dLag_num, dLag_den;

  dLag = 0.0;

  for (i=0; i<x_lag.size(); i++) {
    if (i!=mode) {
      dLag_num = 1.0;
      dLag_den = 1.0;

      for (j=0; j<x_lag.size(); j++) {
        if (j!=mode && j!=i) {
          dLag_num = dLag_num*(y-x_lag[j]);
        }

        if (j!=mode) {
          dLag_den = dLag_den*(x_lag[mode]-x_lag[j]);
        }
      }

      dLag = dLag+(dLag_num/dLag_den);
    }
  }

  return dLag;
}

std::vector<double> adjoint(const std::vector<double> &mat, unsigned int size)
{
  std::vector<double> adj(size*size);

  int signRow = -1;
  std::vector<double> Minor((size-1)*(size-1));
  for (int row=0; row<size; row++) {
    signRow *= -1;
    int sign = -1*signRow;
    for (int col=0; col<size; col++) {
      sign *= -1;
      // Setup the minor matrix (expanding along row, col)
      int i0 = 0;
      for (int i=0; i<size; i++) {
        if (i==row) continue;
        int j0 = 0;
        for (int j=0; j<size; j++) {
          if (j==col) continue;
          Minor[i0*(size-1)+j0] = mat[i*size+j];
          j0++;
        }
        i0++;
      }
      // Recall: adjoint is TRANSPOSE of cofactor matrix
      adj[col*size+row] = sign*determinant(Minor,size-1);
    }
  }

  return adj;
}

double determinant(const std::vector<double> &data, unsigned int size)
{
  if (size == 1) {
    return data[0];
  }
  else if (size == 2) {
    // Base case
    return data[0]*data[3] - data[1]*data[2];
  }
  else {
    // Use minor-matrix recursion
    double Det = 0;
    int sign = -1;
    std::vector<double> Minor((size-1)*(size-1));
    for (int row=0; row<size; row++) {
      sign *= -1;
      // Setup the minor matrix (expanding along first column)
      int i0 = 0;
      for (int i=0; i<size; i++) {
        if (i==row) continue;
        for (int j=1; j<size; j++) {
          Minor[i0*(size-1)+j-1] = data[i*size+j];
        }
        i0++;
      }
      // Add in the minor's determinant
      Det += sign*determinant(Minor,size-1)*data[row*size+0];
    }

    return Det;
  }
}


void getBoundingBox(double *pts, int nPts, int nDims, double *bbox)
{
  for (int i=0; i<nDims; i++) {
    bbox[i]       =  INFINITY;
    bbox[nDims+i] = -INFINITY;
  }

  for (int i=0; i<nPts; i++) {
    for (int dim=0; dim<nDims; dim++) {
      bbox[dim]       = std::min(bbox[dim],      pts[i*nDims+dim]);
      bbox[nDims+dim] = std::max(bbox[nDims+dim],pts[i*nDims+dim]);
    }
  }
}

std::vector<int> gmsh_to_structured_quad(unsigned int nNodes)
{
  std::vector<int> gmsh_to_ijk(nNodes,0);

  /* Lagrange Elements (or linear serendipity) */
  if (nNodes != 8)
  {
    int nNodes1D = sqrt(nNodes);

    if (nNodes1D * nNodes1D != nNodes)
      ThrowException("nNodes must be a square number.");

    int nLevels = nNodes1D / 2;

    /* Set shape values via recursive strategy (from Flurry) */
    int node = 0;
    for (int i = 0; i < nLevels; i++)
    {
      /* Corner Nodes */
      int i2 = (nNodes1D - 1) - i;
      gmsh_to_ijk[node]     = i  + nNodes1D * i;
      gmsh_to_ijk[node + 1] = i2 + nNodes1D * i;
      gmsh_to_ijk[node + 2] = i2 + nNodes1D * i2;
      gmsh_to_ijk[node + 3] = i  + nNodes1D * i2;

      node += 4;

      int nEdgeNodes = nNodes1D - 2 * (i + 1);
      for (int j = 0; j < nEdgeNodes; j++)
      {
        gmsh_to_ijk[node + j]                = i+1+j  + nNodes1D * i;
        gmsh_to_ijk[node + nEdgeNodes + j]   = i2     + nNodes1D * (i+1+j);
        gmsh_to_ijk[node + 2*nEdgeNodes + j] = i2-1-j + nNodes1D * i2;
        gmsh_to_ijk[node + 3*nEdgeNodes + j] = i      + nNodes1D * (i2-1-j);
      }

      node += 4 * nEdgeNodes;
    }

    /* Add center node in odd case */
    if (nNodes1D % 2 != 0)
    {
      gmsh_to_ijk[nNodes - 1] = nNodes1D/2 + nNodes1D * (nNodes1D/2);
    }
  }

  /* 8-node Serendipity Element */
  else
  {
    gmsh_to_ijk[0] = 0; gmsh_to_ijk[1] = 2;  gmsh_to_ijk[2] = 7;
    gmsh_to_ijk[3] = 5; gmsh_to_ijk[4] = 1;  gmsh_to_ijk[5] = 3;
    gmsh_to_ijk[6] = 4; gmsh_to_ijk[7] = 6;
  }

  return gmsh_to_ijk;
}

std::vector<int> structured_to_gmsh_quad(unsigned int nNodes)
{
  auto gmsh2ijk = gmsh_to_structured_quad(nNodes);

  return reverse_map(gmsh2ijk);
}

std::vector<int> structured_to_gmsh_hex(unsigned int nNodes)
{
  auto gmsh2ijk = gmsh_to_structured_hex(nNodes);

  return reverse_map(gmsh2ijk);
}

std::vector<int> gmsh_to_structured_hex(unsigned int nNodes)
{
  std::vector<int> gmsh_to_ijk(nNodes,0);

  int nSide = cbrt(nNodes);

  if (nSide*nSide*nSide != nNodes)
  {
    std::cout << "nNodes = " << nNodes << std::endl;
    ThrowException("For Lagrange hex of order N, must have (N+1)^3 shape points.");
  }

  std::vector<double> xlist(nSide);
  double dxi = 2./(nSide-1);

  for (int i=0; i<nSide; i++)
    xlist[i] = -1. + i*dxi;

  int nLevels = nSide / 2;
  int isOdd = nSide % 2;

  /* Recursion for all high-order Lagrange elements:
             * 8 corners, each edge's points, interior face points, volume points */
  int nPts = 0;
  for (int i = 0; i < nLevels; i++) {
    // Corners
    int i2 = (nSide-1) - i;
    gmsh_to_ijk[nPts+0] = i  + nSide * (i  + nSide * i);
    gmsh_to_ijk[nPts+1] = i2 + nSide * (i  + nSide * i);
    gmsh_to_ijk[nPts+2] = i2 + nSide * (i2 + nSide * i);
    gmsh_to_ijk[nPts+3] = i  + nSide * (i2 + nSide * i);
    gmsh_to_ijk[nPts+4] = i  + nSide * (i  + nSide * i2);
    gmsh_to_ijk[nPts+5] = i2 + nSide * (i  + nSide * i2);
    gmsh_to_ijk[nPts+6] = i2 + nSide * (i2 + nSide * i2);
    gmsh_to_ijk[nPts+7] = i  + nSide * (i2 + nSide * i2);
    nPts += 8;

    // Edges
    int nSide2 = nSide - 2 * (i+1);
    for (int j = 0; j < nSide2; j++) {
      // Edges around 'bottom'
      gmsh_to_ijk[nPts+0*nSide2+j] = i+1+j  + nSide * (i     + nSide * i);
      gmsh_to_ijk[nPts+3*nSide2+j] = i2     + nSide * (i+1+j + nSide * i);
      gmsh_to_ijk[nPts+5*nSide2+j] = i2-1-j + nSide * (i2    + nSide * i);
      gmsh_to_ijk[nPts+1*nSide2+j] = i      + nSide * (i+1+j + nSide * i);

      // 'Vertical' edges
      gmsh_to_ijk[nPts+2*nSide2+j] = i  + nSide * (i  + nSide * (i+1+j));
      gmsh_to_ijk[nPts+4*nSide2+j] = i2 + nSide * (i  + nSide * (i+1+j));
      gmsh_to_ijk[nPts+6*nSide2+j] = i2 + nSide * (i2 + nSide * (i+1+j));
      gmsh_to_ijk[nPts+7*nSide2+j] = i  + nSide * (i2 + nSide * (i+1+j));

      // Edges around 'top'
      gmsh_to_ijk[nPts+ 8*nSide2+j] = i+1+j  + nSide * (i     + nSide * i2);
      gmsh_to_ijk[nPts+10*nSide2+j] = i2     + nSide * (i+1+j + nSide * i2);
      gmsh_to_ijk[nPts+11*nSide2+j] = i2-1-j + nSide * (i2    + nSide * i2);
      gmsh_to_ijk[nPts+ 9*nSide2+j] = i      + nSide * (i+1+j + nSide * i2);
    }
    nPts += 12*nSide2;

    /* --- Faces [Use recursion from quadrilaterals] --- */

    int nLevels2 = nSide2 / 2;
    int isOdd2 = nSide2 % 2;

    // --- Bottom face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = j  + nSide * (j  + nSide * i);
      gmsh_to_ijk[nPts+1] = j  + nSide * (j2 + nSide * i);
      gmsh_to_ijk[nPts+2] = j2 + nSide * (j2 + nSide * i);
      gmsh_to_ijk[nPts+3] = j2 + nSide * (j  + nSide * i);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = j      + nSide * (j+1+k  + nSide * i);
        gmsh_to_ijk[nPts+1*nSide3+k] = j+1+k  + nSide * (j2     + nSide * i);
        gmsh_to_ijk[nPts+2*nSide3+k] = j2     + nSide * (j2-1-k + nSide * i);
        gmsh_to_ijk[nPts+3*nSide3+k] = j2-1-k + nSide * (j      + nSide * i);
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = nSide/2 +  nSide*(nSide/2) +  nSide*nSide*i;
      nPts += 1;
    }

    // --- Front face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = j  + nSide * (i + nSide * j);
      gmsh_to_ijk[nPts+1] = j2 + nSide * (i + nSide * j);
      gmsh_to_ijk[nPts+2] = j2 + nSide * (i + nSide * j2);
      gmsh_to_ijk[nPts+3] = j  + nSide * (i + nSide * j2);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = j+1+k  + nSide * (i + nSide * j);
        gmsh_to_ijk[nPts+1*nSide3+k] = j2     + nSide * (i + nSide * (j+1+k));
        gmsh_to_ijk[nPts+2*nSide3+k] = j2-1-k + nSide * (i + nSide * j2);
        gmsh_to_ijk[nPts+3*nSide3+k] = j      + nSide * (i + nSide * (j2-1-k));
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = nSide/2 + nSide*(i + nSide*(nSide/2));
      nPts += 1;
    }

    // --- Left face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = i + nSide * (j  + nSide * j);
      gmsh_to_ijk[nPts+1] = i + nSide * (j  + nSide * j2);
      gmsh_to_ijk[nPts+2] = i + nSide * (j2 + nSide * j2);
      gmsh_to_ijk[nPts+3] = i + nSide * (j2 + nSide * j);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = i + nSide * (j      + nSide * (j+1+k));
        gmsh_to_ijk[nPts+1*nSide3+k] = i + nSide * (j+1+k  + nSide * j2);
        gmsh_to_ijk[nPts+2*nSide3+k] = i + nSide * (j2     + nSide * (j2-1-k));
        gmsh_to_ijk[nPts+3*nSide3+k] = i + nSide * (j2-1-k + nSide * j);
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = i + nSide * (nSide/2 + nSide * (nSide/2));
      nPts += 1;
    }

    // --- Right face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = i2 + nSide * (j  + nSide * j);
      gmsh_to_ijk[nPts+1] = i2 + nSide * (j2 + nSide * j);
      gmsh_to_ijk[nPts+2] = i2 + nSide * (j2 + nSide * j2);
      gmsh_to_ijk[nPts+3] = i2 + nSide * (j  + nSide * j2);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = i2 + nSide * (j+1+k  + nSide * j);
        gmsh_to_ijk[nPts+1*nSide3+k] = i2 + nSide * (j2     + nSide * (j+1+k));
        gmsh_to_ijk[nPts+2*nSide3+k] = i2 + nSide * (j2-1-k + nSide * j2);
        gmsh_to_ijk[nPts+3*nSide3+k] = i2 + nSide * (j      + nSide * (j2-1-k));
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = i2 + nSide * (nSide/2 + nSide * (nSide/2));
      nPts += 1;
    }

    // --- Back face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = j2 + nSide * (i2 + nSide * j);
      gmsh_to_ijk[nPts+1] = j  + nSide * (i2 + nSide * j);
      gmsh_to_ijk[nPts+2] = j  + nSide * (i2 + nSide * j2);
      gmsh_to_ijk[nPts+3] = j2 + nSide * (i2 + nSide * j2);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = j2-1-k + nSide * (i2 + nSide*j);
        gmsh_to_ijk[nPts+1*nSide3+k] = j      + nSide * (i2 + nSide*(j+1+k));
        gmsh_to_ijk[nPts+2*nSide3+k] = j+1+k  + nSide * (i2 + nSide*j2);
        gmsh_to_ijk[nPts+3*nSide3+k] = j2     + nSide * (i2 + nSide*(j2-1-k));
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = nSide/2 + nSide * (i2 + nSide * (nSide/2));
      nPts += 1;
    }

    // --- Top face ---
    for (int j0 = 0; j0 < nLevels2; j0++) {
      // Corners
      int j = j0 + i + 1;
      int j2 = i + 1 + (nSide2-1) - j0;
      gmsh_to_ijk[nPts+0] = j  + nSide * (j  + nSide * i2);
      gmsh_to_ijk[nPts+1] = j2 + nSide * (j  + nSide * i2);
      gmsh_to_ijk[nPts+2] = j2 + nSide * (j2 + nSide * i2);
      gmsh_to_ijk[nPts+3] = j  + nSide * (j2 + nSide * i2);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide3 = nSide2 - 2 * (j0+1);
      for (int k = 0; k < nSide3; k++) {
        gmsh_to_ijk[nPts+0*nSide3+k] = j+1+k  + nSide * (j      + nSide * i2);
        gmsh_to_ijk[nPts+1*nSide3+k] = j2     + nSide * (j+1+k  + nSide * i2);
        gmsh_to_ijk[nPts+2*nSide3+k] = j2-1-k + nSide * (j2     + nSide * i2);
        gmsh_to_ijk[nPts+3*nSide3+k] = j      + nSide * (j2-1-k + nSide * i2);
      }
      nPts += 4*nSide3;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd2) {
      gmsh_to_ijk[nPts] = nSide/2 + nSide * (nSide/2 +  nSide * i2);
      nPts += 1;
    }
  }

  // Center node for even-ordered Lagrange quads (odd value of nSide)
  if (isOdd) {
    gmsh_to_ijk[nNodes-1] = nSide/2 + nSide * (nSide/2 + nSide * (nSide/2));
  }

  return gmsh_to_ijk;
}

std::vector<int> reverse_map(const std::vector<int> &map1)
{
  auto map2 = map1;
  for (int i = 0; i < map1.size(); i++)
    map2[i] = findFirst(map1, i);

  return map2;
}

std::vector<int> get_int_list(int N, int start)
{
  std::vector<int> list(N);
  for (int i = 0; i < N; i++)
    list[i] = start + i;

  return list;
}

std::vector<uint> get_int_list(uint N, uint start)
{
  std::vector<uint> list(N);
  for (uint i = 0; i < N; i++)
    list[i] = start + i;

  return list;
}

bool getRefLocNewton(double *xv, double *in_xyz, double *out_rst, int nNodes, int nDims)
{
  // First, do a quick check to see if the point is even close to being in the element
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  xmin = ymin = zmin =  1e15;
  xmax = ymax = zmax = -1e15;
  double eps = 1e-10;

  std::vector<double> box(6);
  getBoundingBox(xv, nNodes, nDims, box.data());
  xmin = box[0];  ymin = box[1];  zmin = box[2];
  xmax = box[3];  ymax = box[4];  zmax = box[5];

  point pos = point(in_xyz);
  if (pos.x < xmin-eps || pos.y < ymin-eps || pos.z < zmin-eps ||
      pos.x > xmax+eps || pos.y > ymax+eps || pos.z > zmax+eps) {
    // Point does not lie within cell - return an obviously bad ref position
    for (int i = 0; i < nDims; i++) out_rst[i] = 99.;
    return false;
  }

  // Use a relative tolerance to handle extreme grids
  double h = std::min(xmax-xmin,ymax-ymin);
  if (nDims==3) h = std::min(h,zmax-zmin);

  double tol = 1e-12*h;

  std::vector<double> shape(nNodes);
  std::vector<double> dshape(nNodes*nDims);
  std::vector<double> grad(nDims*nDims);

  int iter = 0;
  int iterMax = 20;
  double norm = 1;

  // Starting location: {0,0,0}
  for (int i = 0; i < nDims; i++) out_rst[i] = 0;
  point loc = point(out_rst,nDims);

  while (norm > tol && iter<iterMax) {
    if (nDims == 2) {
      shape_quad(loc,shape.data(),nNodes);
      dshape_quad(loc,dshape.data(),nNodes);
    } else {
      shape_hex(loc,shape.data(),nNodes);
      dshape_hex(loc,dshape.data(),nNodes);
    }

    point dx = pos;
    grad.assign(nDims*nDims,0.);

    for (int n=0; n<nNodes; n++) {
      for (int i=0; i<nDims; i++) {
        for (int j=0; j<nDims; j++) {
          grad[i*nDims+j] += xv[n*nDims+i]*dshape[n*nDims+j];
        }
        dx[i] -= shape[n]*xv[n*nDims+i];
      }
    }

    double detJ = determinant(grad,nDims);

    auto ginv = adjoint(grad,nDims);

    point delta = {0,0,0};
    for (int i=0; i<nDims; i++)
      for (int j=0; j<nDims; j++)
        delta[i] += ginv[i*nDims+j]*dx[j]/detJ;

    norm = 0;
    for (int i=0; i<nDims; i++) {
      norm += dx[i]*dx[i];
      loc[i] += delta[i];
      loc[i] = std::max(std::min(loc[i],1.01),-1.01);
    }

    iter++;
  }

  if (std::max( std::abs(loc[0]), std::max( std::abs(loc[1]), std::abs(loc[2]) ) ) <= 1. + 1e-10)
    return true;
  else
    return false;
}


double computeVolume(double *xv, int nNodes, int nDims)
{
  int order;
  std::vector<point> locSpts;

  if (nDims == 2)
  {
    order = std::max((int)std::sqrt(nNodes)-1, 0);
    locSpts = getLocSpts(QUAD,order,std::string("Legendre"));
  }
  else
  {
    order = std::max((int)cbrt(nNodes)-1, 0);
    locSpts = getLocSpts(HEX,order,std::string("Legendre"));
  }

  auto weights = getQptWeights(order, nDims);

  uint nSpts = locSpts.size();

  std::vector<double> shape(nSpts*nNodes);
  std::vector<double> dshape(nSpts*nNodes*nDims);

  if (nDims == 2)
  {
    for (uint spt = 0; spt < nSpts; spt++)
    {
      shape_quad(locSpts[spt], &shape[spt*nNodes], nNodes);
      dshape_quad(locSpts[spt], &dshape[spt*nNodes*nDims], nNodes);
    }
  }
  else
  {
    for (uint spt = 0; spt < nSpts; spt++)
    {
      shape_hex(locSpts[spt], &shape[spt*nNodes], nNodes);
      dshape_hex(locSpts[spt], &dshape[spt*nNodes*nDims], nNodes);
    }
  }

  std::vector<double> jaco(nDims*nDims);
  double vol = 0.;

  for (uint spt = 0; spt < nSpts; spt++)
  {
    jaco.assign(jaco.size(), 0);
    for (uint n = 0; n < nNodes; n++)
      for (uint d1 = 0; d1 < nDims; d1++)
        for (uint d2 = 0; d2 < nDims; d2++)
          jaco[d1*nDims+d2] += dshape[(spt*nNodes+n)*nDims+d2] * xv[n*nDims+d1];

    double detJac = 0;
    if (nDims == 2)
    {
      detJac = jaco[0*nDims+0]*jaco[1*nDims+1] - jaco[1*nDims+0]*jaco[0*nDims+1];
    }
    else
    {
      double xr = jaco[0*nDims+0];   double xs = jaco[0*nDims+1];   double xt = jaco[0*nDims+2];
      double yr = jaco[1*nDims+0];   double ys = jaco[1*nDims+1];   double yt = jaco[1*nDims+2];
      double zr = jaco[2*nDims+0];   double zs = jaco[2*nDims+1];   double zt = jaco[2*nDims+2];
      detJac = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);
    }
    if (detJac<0) FatalError("Negative Jacobian at quadrature point.");

    vol += detJac * weights[spt];
  }

  return vol;
}

void shape_quad(const point &in_rs, std::vector<double> &out_shape, int nNodes)
{
  out_shape.resize(nNodes);
  shape_quad(in_rs, out_shape.data(), nNodes);
}

void shape_quad(const point &in_rs, double* out_shape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;

  int nSide = sqrt(nNodes);

  if (nSide*nSide != nNodes)
    FatalError("For Lagrange quad of order N, must have (N+1)^2 shape points.");

  std::vector<double> xlist(nSide);
  double dxi = 2./(nSide-1);

  for (int i=0; i<nSide; i++)
    xlist[i] = -1. + i*dxi;

  int nLevels = nSide / 2;
  int isOdd = nSide % 2;

  // Pre-compute Lagrange values to avoid redundant calculations
  std::vector<double> lag_i(nSide), lag_j(nSide);
  for (int i = 0; i < nSide; i++)
  {
    lag_i[i] = Lagrange(xlist,  xi, i);
    lag_j[i] = Lagrange(xlist, eta, i);
  }

  /* Recursion for all high-order Lagrange elements:
   * 4 corners, each edge's points, interior points */
  int nPts = 0;
  for (int i = 0; i < nLevels; i++) {
    // Corners
    int i2 = (nSide-1) - i;
    out_shape[nPts+0] = lag_i[i]  * lag_j[i];
    out_shape[nPts+1] = lag_i[i2] * lag_j[i];
    out_shape[nPts+2] = lag_i[i2] * lag_j[i2];
    out_shape[nPts+3] = lag_i[i]  * lag_j[i2];
    nPts += 4;

    // Edges: Bottom, right, top, left
    int nSide2 = nSide - 2 * (i+1);
    for (int j = 0; j < nSide2; j++) {
      out_shape[nPts+0*nSide2+j] = lag_i[i+1+j]  * lag_j[i];
      out_shape[nPts+1*nSide2+j] = lag_i[i2]     * lag_j[i+1+j];
      out_shape[nPts+2*nSide2+j] = lag_i[i2-1-j] * lag_j[i2];
      out_shape[nPts+3*nSide2+j] = lag_i[i]      * lag_j[i2-1-j];
    }
    nPts += 4*nSide2;
  }

  // Center node for even-ordered Lagrange quads (odd value of nSide)
  if (isOdd) {
    out_shape[nNodes-1] = lag_i[nSide/2] * lag_j[nSide/2];
  }
}

void shape_hex(const point &in_rst, std::vector<double> &out_shape, int nNodes)
{
  out_shape.resize(nNodes);
  shape_hex(in_rst, out_shape.data(), nNodes);
}

void shape_hex(const point &in_rst, double* out_shape, int nNodes)
{
  double xi  = in_rst.x;
  double eta = in_rst.y;
  double mu = in_rst.z;

  if (nNodes == 20) {
    double XI[8]  = {-1,1,1,-1,-1,1,1,-1};
    double ETA[8] = {-1,-1,1,1,-1,-1,1,1};
    double MU[8]  = {-1,-1,-1,-1,1,1,1,1};
    // Corner nodes
    for (int i=0; i<8; i++) {
      out_shape[i] = .125*(1+xi*XI[i])*(1+eta*ETA[i])*(1+mu*MU[i])*(xi*XI[i]+eta*ETA[i]+mu*MU[i]-2);
    }
    // Edge nodes, xi = 0
    out_shape[8]  = .25*(1-xi*xi)*(1-eta)*(1-mu);
    out_shape[10] = .25*(1-xi*xi)*(1+eta)*(1-mu);
    out_shape[16] = .25*(1-xi*xi)*(1-eta)*(1+mu);
    out_shape[18] = .25*(1-xi*xi)*(1+eta)*(1+mu);
    // Edge nodes, eta = 0
    out_shape[9]  = .25*(1-eta*eta)*(1+xi)*(1-mu);
    out_shape[11] = .25*(1-eta*eta)*(1-xi)*(1-mu);
    out_shape[17] = .25*(1-eta*eta)*(1+xi)*(1+mu);
    out_shape[19] = .25*(1-eta*eta)*(1-xi)*(1+mu);
    // Edge Nodes, mu = 0
    out_shape[12] = .25*(1-mu*mu)*(1-xi)*(1-eta);
    out_shape[13] = .25*(1-mu*mu)*(1+xi)*(1-eta);
    out_shape[14] = .25*(1-mu*mu)*(1+xi)*(1+eta);
    out_shape[15] = .25*(1-mu*mu)*(1-xi)*(1+eta);
  }
  else {
    int nSide = cbrt(nNodes);

    if (nSide*nSide*nSide != nNodes)
    {
      std::cout << "nNodes = " << nNodes << std::endl;
      ThrowException("For Lagrange hex of order N, must have (N+1)^3 shape points.");
    }

    std::vector<double> xlist(nSide);
    double dxi = 2./(nSide-1);

    for (int i=0; i<nSide; i++)
      xlist[i] = -1. + i*dxi;

    auto ijk2gmsh = structured_to_gmsh_hex(nNodes);

    // Pre-compute Lagrange function values to avoid redundant calculations
    std::vector<double> lag_i(nSide), lag_j(nSide), lag_k(nSide);

#pragma omp parallel for
    for (int i = 0; i < nSide; i++)
    {
      lag_i[i] = Lagrange(xlist,  xi, i);
      lag_j[i] = Lagrange(xlist, eta, i);
      lag_k[i] = Lagrange(xlist,  mu, i);
    }

#pragma omp parallel for collapse(3)
    for (int k = 0; k < nSide; k++)
      for (int j = 0; j < nSide; j++)
        for (int i = 0; i < nSide; i++)
          out_shape[ijk2gmsh[i+nSide*(j+nSide*k)]] = lag_i[i] * lag_j[j] * lag_k[k];
  }
}

void dshape_quad(const std::vector<point> loc_pts, double* out_dshape, int nNodes)
{
  for (int i = 0; i < loc_pts.size(); i++)
    dshape_quad(loc_pts[i], &out_dshape[i*nNodes*2], nNodes);
}

void dshape_quad(const point &in_rs, double* out_dshape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;

  int nSide = sqrt(nNodes);

  if (nSide*nSide != nNodes)
    FatalError("For Lagrange quad of order N, must have (N+1)^2 shape points.");

  std::vector<double> xlist(nSide);
  double dxi = 2./(nSide-1);

  for (int i=0; i<nSide; i++)
    xlist[i] = -1. + i*dxi;

  int nLevels = nSide / 2;
  int isOdd = nSide % 2;

  /* Recursion for all high-order Lagrange elements:
     * 4 corners, each edge's points, interior points */
  int nPts = 0;
  for (int i = 0; i < nLevels; i++) {
    // Corners
    int i2 = (nSide-1) - i;
    out_dshape[2*(nPts+0)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i);
    out_dshape[2*(nPts+1)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i);
    out_dshape[2*(nPts+2)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2);
    out_dshape[2*(nPts+3)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2);

    out_dshape[2*(nPts+0)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i);
    out_dshape[2*(nPts+1)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i);
    out_dshape[2*(nPts+2)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i2);
    out_dshape[2*(nPts+3)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i2);
    nPts += 4;

    // Edges
    int nSide2 = nSide - 2 * (i+1);
    for (int j = 0; j < nSide2; j++) {
      out_dshape[2*(nPts+0*nSide2+j)] = dLagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i);
      out_dshape[2*(nPts+1*nSide2+j)] = dLagrange(xlist, xi, i2)   * Lagrange(xlist, eta, i+1+j);
      out_dshape[2*(nPts+2*nSide2+j)] = dLagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2);
      out_dshape[2*(nPts+3*nSide2+j)] = dLagrange(xlist, xi, i)    * Lagrange(xlist, eta, i2-1-j);

      out_dshape[2*(nPts+0*nSide2+j)+1] = Lagrange(xlist, xi, i+1+j)  * dLagrange(xlist, eta, i);
      out_dshape[2*(nPts+1*nSide2+j)+1] = Lagrange(xlist, xi, i2)   * dLagrange(xlist, eta, i+1+j);
      out_dshape[2*(nPts+2*nSide2+j)+1] = Lagrange(xlist, xi, i2-1-j) * dLagrange(xlist, eta, i2);
      out_dshape[2*(nPts+3*nSide2+j)+1] = Lagrange(xlist, xi, i)    * dLagrange(xlist, eta, i2-1-j);
    }
    nPts += 4*nSide2;
  }

  // Center node for even-ordered Lagrange quads (odd value of nSide)
  if (isOdd) {
    out_dshape[2*(nNodes-1)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2);
    out_dshape[2*(nNodes-1)+1] = Lagrange(xlist, xi, nSide/2) * dLagrange(xlist, eta, nSide/2);
  }
}

void dshape_hex(const std::vector<point> &loc_pts, double* out_dshape, int nNodes)
{
  for (int i = 0; i < loc_pts.size(); i++) {
    point pt = loc_pts[i];
    dshape_hex(pt, &out_dshape[i*nNodes*3], nNodes);
  }
}

void dshape_hex(const point &in_rst, double* out_dshape, int nNodes)
{
  double xi  = in_rst.x;
  double eta = in_rst.y;
  double mu = in_rst.z;

  if (nNodes == 20) {
    /* Quadratic Serendiptiy Hex */
    double XI[8]  = {-1,1,1,-1,-1,1,1,-1};
    double ETA[8] = {-1,-1,1,1,-1,-1,1,1};
    double MU[8]  = {-1,-1,-1,-1,1,1,1,1};
    // Corner Nodes
    for (int i=0; i<8; i++) {
      out_dshape[3*i+0] = .125*XI[i] *(1+eta*ETA[i])*(1 + mu*MU[i])*(2*xi*XI[i] +   eta*ETA[i] +   mu*MU[i]-1);
      out_dshape[3*i+1] = .125*ETA[i]*(1 + xi*XI[i])*(1 + mu*MU[i])*(  xi*XI[i] + 2*eta*ETA[i] +   mu*MU[i]-1);
      out_dshape[3*i+2] = .125*MU[i] *(1 + xi*XI[i])*(1+eta*ETA[i])*(  xi*XI[i] +   eta*ETA[i] + 2*mu*MU[i]-1);
    }
    // Edge Nodes, xi = 0
    out_dshape[ 3*8+0] = -.5*xi*(1-eta)*(1-mu);  out_dshape[ 3*8+1] = -.25*(1-xi*xi)*(1-mu);  out_dshape[ 3*8+2] = -.25*(1-xi*xi)*(1-eta);
    out_dshape[3*10+0] = -.5*xi*(1+eta)*(1-mu);  out_dshape[3*10+1] =  .25*(1-xi*xi)*(1-mu);  out_dshape[3*10+2] = -.25*(1-xi*xi)*(1+eta);
    out_dshape[3*16+0] = -.5*xi*(1-eta)*(1+mu);  out_dshape[3*16+1] = -.25*(1-xi*xi)*(1+mu);  out_dshape[3*16+2] =  .25*(1-xi*xi)*(1-eta);
    out_dshape[3*18+0] = -.5*xi*(1+eta)*(1+mu);  out_dshape[3*18+1] =  .25*(1-xi*xi)*(1+mu);  out_dshape[3*18+2] =  .25*(1-xi*xi)*(1+eta);
    // Edge Nodes, eta = 0
    out_dshape[ 3*9+1] = -.5*eta*(1+xi)*(1-mu);  out_dshape[ 3*9+0] =  .25*(1-eta*eta)*(1-mu);  out_dshape[ 3*9+2] = -.25*(1-eta*eta)*(1+xi);
    out_dshape[3*11+1] = -.5*eta*(1-xi)*(1-mu);  out_dshape[3*11+0] = -.25*(1-eta*eta)*(1-mu);  out_dshape[3*11+2] = -.25*(1-eta*eta)*(1-xi);
    out_dshape[3*17+1] = -.5*eta*(1+xi)*(1+mu);  out_dshape[3*17+0] =  .25*(1-eta*eta)*(1+mu);  out_dshape[3*17+2] =  .25*(1-eta*eta)*(1+xi);
    out_dshape[3*19+1] = -.5*eta*(1-xi)*(1+mu);  out_dshape[3*19+0] = -.25*(1-eta*eta)*(1+mu);  out_dshape[3*19+2] =  .25*(1-eta*eta)*(1-xi);
    // Edge Nodes, mu = 0;
    out_dshape[3*12+2] = -.5*mu*(1-xi)*(1-eta);  out_dshape[3*12+0] = -.25*(1-mu*mu)*(1-eta);  out_dshape[3*12+1] = -.25*(1-mu*mu)*(1-xi);
    out_dshape[3*13+2] = -.5*mu*(1+xi)*(1-eta);  out_dshape[3*13+0] =  .25*(1-mu*mu)*(1-eta);  out_dshape[3*13+1] = -.25*(1-mu*mu)*(1+xi);
    out_dshape[3*14+2] = -.5*mu*(1+xi)*(1+eta);  out_dshape[3*14+0] =  .25*(1-mu*mu)*(1+eta);  out_dshape[3*14+1] =  .25*(1-mu*mu)*(1+xi);
    out_dshape[3*15+2] = -.5*mu*(1-xi)*(1+eta);  out_dshape[3*15+0] = -.25*(1-mu*mu)*(1+eta);  out_dshape[3*15+1] =  .25*(1-mu*mu)*(1-xi);
  }
  else {
    int nSide = cbrt(nNodes);

    if (nSide*nSide*nSide != nNodes)
      ThrowException("For Lagrange hex of order N, must have (N+1)^3 shape points.");

    std::vector<double> xlist(nSide);
    double dxi = 2./(nSide-1);

    for (int i=0; i<nSide; i++)
      xlist[i] = -1. + i*dxi;

    auto ijk2gmsh = structured_to_gmsh_hex(nNodes);

    // Pre-compute Lagrange function values to save redundant calculations
    std::vector<double>  lag_i(nSide),  lag_j(nSide),  lag_k(nSide);
    std::vector<double> dlag_i(nSide), dlag_j(nSide), dlag_k(nSide);

#pragma omp parallel for
    for (int i = 0; i < nSide; i++)
    {
      lag_i[i] = Lagrange(xlist,  xi, i);  dlag_i[i] = dLagrange(xlist,  xi, i);
      lag_j[i] = Lagrange(xlist, eta, i);  dlag_j[i] = dLagrange(xlist, eta, i);
      lag_k[i] = Lagrange(xlist,  mu, i);  dlag_k[i] = dLagrange(xlist,  mu, i);
    }

#pragma omp parallel for collapse(3)
    for (int k = 0; k < nSide; k++)
    {
      for (int j = 0; j < nSide; j++)
      {
        for (int i = 0; i < nSide; i++)
        {
          int pt = i + nSide * (j + nSide * k);
          out_dshape[3*ijk2gmsh[pt]+0] = dlag_i[i] *  lag_j[j] *  lag_k[k];
          out_dshape[3*ijk2gmsh[pt]+1] =  lag_i[i] * dlag_j[j] *  lag_k[k];
          out_dshape[3*ijk2gmsh[pt]+2] =  lag_i[i] *  lag_j[j] * dlag_k[k];
        }
      }
    }
  }
}

} // namespace tg_funcs
