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

#include "funcs.hpp"

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
  for (auto &val:vec) cout << val << ", ";
  return os;
}

std::ostream& operator<<(std::ostream &os, const std::vector<double> &vec)
{
  for (auto &val:vec) cout << val << ", ";
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
          Minor(i0*(size-1)+j0) = mat(i*size+j);
          j0++;
        }
        i0++;
      }
      // Recall: adjoint is TRANSPOSE of cofactor matrix
      adj(col*size+row) = sign*determinant(Minor,size-1);
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
    for (int row=0; row<this->dims[0]; row++) {
      sign *= -1;
      // Setup the minor matrix (expanding along first column)
      int i0 = 0;
      for (int i=0; i<size; i++) {
        if (i==row) continue;
        for (int j=1; j<size; j++) {
          Minor(i0,j-1) = data[i*size+j];
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
      bbox[dim]       = min(bbox[dim],      pts[i*nDims+dim]);
      bbox[nDims+dim] = max(bbox[nDims+dim],pts[i*nDims+dim]);
    }
  }
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
  double h = min(xmax-xmin,ymax-ymin);
  if (nDims==3) h = min(h,zmax-zmin);

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

    double detJ = grad.det();

    auto ginv = grad.adjoint();

    point delta = {0,0,0};
    for (int i=0; i<nDims; i++)
      for (int j=0; j<nDims; j++)
        delta[i] += ginv[i*nDims+j]*dx[j]/detJ;

    norm = 0;
    for (int i=0; i<nDims; i++) {
      norm += dx[i]*dx[i];
      loc[i] += delta[i];
      loc[i] = max(min(loc[i],1.01),-1.01);
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
    jaco.initializeToZero();
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

  /* Recursion for all high-order Lagrange elements:
     * 4 corners, each edge's points, interior points */
  int nPts = 0;
  for (int i = 0; i < nLevels; i++) {
    // Corners
    int i2 = (nSide-1) - i;
    out_shape[nPts+0] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i);
    out_shape[nPts+1] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i);
    out_shape[nPts+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2);
    out_shape[nPts+3] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2);
    nPts += 4;

    // Edges: Bottom, right, top, left
    int nSide2 = nSide - 2 * (i+1);
    for (int j = 0; j < nSide2; j++) {
      out_shape[nPts+0*nSide2+j] = Lagrange(xlist, xi, i+1+j) * Lagrange(xlist, eta, i);
      out_shape[nPts+1*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i+1+j);
      out_shape[nPts+2*nSide2+j] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2);
      out_shape[nPts+3*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2-1-j);
    }
    nPts += 4*nSide2;
  }

  // Center node for even-ordered Lagrange quads (odd value of nSide)
  if (isOdd) {
    out_shape[nNodes-1] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2);
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
      FatalError("For Lagrange hex of order N, must have (N+1)^3 shape points.");

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
      out_shape[nPts+0] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i);
      out_shape[nPts+1] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i);
      out_shape[nPts+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_shape[nPts+3] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_shape[nPts+4] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i2);
      out_shape[nPts+5] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i2);
      out_shape[nPts+6] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);
      out_shape[nPts+7] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);
      nPts += 8;

      // Edges
      int nSide2 = nSide - 2 * (i+1);
      for (int j = 0; j < nSide2; j++) {
        // Edges around 'bottom'
        out_shape[nPts+0*nSide2+j] = Lagrange(xlist, xi, i+1+j) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i);
        out_shape[nPts+3*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i);
        out_shape[nPts+5*nSide2+j] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
        out_shape[nPts+1*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i);

        // 'Vertical' edges
        out_shape[nPts+2*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i+1+j);
        out_shape[nPts+4*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i+1+j);
        out_shape[nPts+6*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);
        out_shape[nPts+7*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);

        // Edges around 'top'
        out_shape[nPts+8*nSide2+j] = Lagrange(xlist, xi, i+1+j) * Lagrange(xlist, eta, i) * Lagrange(xlist, mu, i2);
        out_shape[nPts+10*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);
        out_shape[nPts+11*nSide2+j] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);
        out_shape[nPts+9*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);
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
        out_shape[nPts+0] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i);
        out_shape[nPts+1] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_shape[nPts+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_shape[nPts+3] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j+1+k) * Lagrange(xlist, mu, i);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, xi, j+1+k) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, i);
        nPts += 1;
      }

      // --- Front face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_shape[nPts+0] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i);
        out_shape[nPts+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i);
        out_shape[nPts+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
        out_shape[nPts+3] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, xi, j+1+k) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j+1+k) * Lagrange(xlist, eta, i);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, eta, i);
        nPts += 1;
      }

      // --- Left face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_shape[nPts+0] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i);
        out_shape[nPts+1] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_shape[nPts+2] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_shape[nPts+3] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j+1+k) * Lagrange(xlist, xi, i);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, eta, j+1+k) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, xi, i);
        nPts += 1;
      }

      // --- Right face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_shape[nPts+0] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i2);
        out_shape[nPts+1] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i2);
        out_shape[nPts+2] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
        out_shape[nPts+3] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, eta, j+1+k) * Lagrange(xlist, mu, j) * Lagrange(xlist, xi, i2);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j+1+k) * Lagrange(xlist, xi, i2);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, eta, j) * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, xi, i2);
        nPts += 1;
      }

      // --- Back face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_shape[nPts+0] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i2);
        out_shape[nPts+1] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i2);
        out_shape[nPts+2] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
        out_shape[nPts+3] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j) * Lagrange(xlist, eta, i2);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, xi, j) * Lagrange(xlist, mu, j+1+k) * Lagrange(xlist, eta, i2);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, xi, j+1+k) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, eta, i2);
        nPts += 1;
      }

      // --- Top face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_shape[nPts+0] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i2);
        out_shape[nPts+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i2);
        out_shape[nPts+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);
        out_shape[nPts+3] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_shape[nPts+0*nSide3+k] = Lagrange(xlist, xi, j+1+k) * Lagrange(xlist, eta, j) * Lagrange(xlist, mu, i2);
          out_shape[nPts+1*nSide3+k] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j+1+k) * Lagrange(xlist, mu, i2);
          out_shape[nPts+2*nSide3+k] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);
          out_shape[nPts+3*nSide3+k] = Lagrange(xlist, xi, j) * Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_shape[nPts] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, i2);
        nPts += 1;
      }
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd) {
      out_shape[nNodes-1] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2);
    }
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
    dshape_hex(loc_pts[i], dshape_tmp[i*nNodes*3], nNodes);
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
      FatalError("For Lagrange hex of order N, must have (N+1)^3 shape points.");

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
      out_dshape[3*(nPts+0)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+1)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+2)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+3)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+4)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+5)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+6)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+7)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);

      out_dshape[3*(nPts+0)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+1)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+2)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+3)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i);
      out_dshape[3*(nPts+4)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+5)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+6)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);
      out_dshape[3*(nPts+7)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i2);

      out_dshape[3*(nPts+0)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i);
      out_dshape[3*(nPts+1)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i);
      out_dshape[3*(nPts+2)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i);
      out_dshape[3*(nPts+3)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i);
      out_dshape[3*(nPts+4)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i2);
      out_dshape[3*(nPts+5)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i2);
      out_dshape[3*(nPts+6)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i2);
      out_dshape[3*(nPts+7)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i2);
      nPts += 8;

      // Edges
      int nSide2 = nSide - 2 * (i+1);
      for (int j = 0; j < nSide2; j++) {
        // Edges around 'bottom'
        out_dshape[3*(nPts+0*nSide2+j)+0] = dLagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i)      * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+3*nSide2+j)+0] = dLagrange(xlist, xi, i2)     * Lagrange(xlist, eta, i+1+j)  * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+5*nSide2+j)+0] = dLagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2)     * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+1*nSide2+j)+0] = dLagrange(xlist, xi, i)      * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i);

        out_dshape[3*(nPts+0*nSide2+j)+1] = Lagrange(xlist, xi, i+1+j)  * dLagrange(xlist, eta, i)      * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+3*nSide2+j)+1] = Lagrange(xlist, xi, i2)     * dLagrange(xlist, eta, i+1+j)  * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+5*nSide2+j)+1] = Lagrange(xlist, xi, i2-1-j) * dLagrange(xlist, eta, i2)     * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+1*nSide2+j)+1] = Lagrange(xlist, xi, i)      * dLagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i);

        out_dshape[3*(nPts+0*nSide2+j)+2] = Lagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i)      * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+3*nSide2+j)+2] = Lagrange(xlist, xi, i2)     * Lagrange(xlist, eta, i+1+j)  * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+5*nSide2+j)+2] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2)     * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+1*nSide2+j)+2] = Lagrange(xlist, xi, i)      * Lagrange(xlist, eta, i+1+j) * dLagrange(xlist, mu, i);

        // 'Vertical' edges
        out_dshape[3*(nPts+2*nSide2+j)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+4*nSide2+j)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+6*nSide2+j)+0] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+7*nSide2+j)+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);

        out_dshape[3*(nPts+2*nSide2+j)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+4*nSide2+j)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i)  * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+6*nSide2+j)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+7*nSide2+j)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i2) * Lagrange(xlist, mu, i+1+j);

        out_dshape[3*(nPts+2*nSide2+j)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+4*nSide2+j)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i)  * dLagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+6*nSide2+j)+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i+1+j);
        out_dshape[3*(nPts+7*nSide2+j)+2] = Lagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2) * dLagrange(xlist, mu, i+1+j);

        // Edges around 'top'
        out_dshape[3*(nPts+ 8*nSide2+j)+0] = dLagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i)     * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+10*nSide2+j)+0] = dLagrange(xlist, xi, i2)     * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+11*nSide2+j)+0] = dLagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2)    * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+ 9*nSide2+j)+0] = dLagrange(xlist, xi, i)      * Lagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);

        out_dshape[3*(nPts+ 8*nSide2+j)+1] = Lagrange(xlist, xi, i+1+j)  * dLagrange(xlist, eta, i)     * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+10*nSide2+j)+1] = Lagrange(xlist, xi, i2)     * dLagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+11*nSide2+j)+1] = Lagrange(xlist, xi, i2-1-j) * dLagrange(xlist, eta, i2)    * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+ 9*nSide2+j)+1] = Lagrange(xlist, xi, i)      * dLagrange(xlist, eta, i+1+j) * Lagrange(xlist, mu, i2);

        out_dshape[3*(nPts+ 8*nSide2+j)+2] = Lagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i)     * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+10*nSide2+j)+2] = Lagrange(xlist, xi, i2)     * Lagrange(xlist, eta, i+1+j) * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+11*nSide2+j)+2] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2)    * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+ 9*nSide2+j)+2] = Lagrange(xlist, xi, i)      * Lagrange(xlist, eta, i+1+j) * dLagrange(xlist, mu, i2);
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
        out_dshape[3*(nPts+0)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+1)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+2)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+3)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, i);

        out_dshape[3*(nPts+0)+1] = Lagrange(xlist, xi, j)  * dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+1)+1] = Lagrange(xlist, xi, j)  * dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+2)+1] = Lagrange(xlist, xi, j2) * dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts+3)+1] = Lagrange(xlist, xi, j2) * dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, i);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, xi, j)  * Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, xi, j)  * Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, i);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = dLagrange(xlist, xi, j)      * Lagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+1*nSide3+k)+0] = dLagrange(xlist, xi, j+1+k)  * Lagrange(xlist, eta, j2)     * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+2*nSide3+k)+0] = dLagrange(xlist, xi, j2)     * Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+3*nSide3+k)+0] = dLagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j)      * Lagrange(xlist, mu, i);

          out_dshape[3*(nPts+0*nSide3+k)+1] = Lagrange(xlist, xi, j)      * dLagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+1*nSide3+k)+1] = Lagrange(xlist, xi, j+1+k)  * dLagrange(xlist, eta, j2)     * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+2*nSide3+k)+1] = Lagrange(xlist, xi, j2)     * dLagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i);
          out_dshape[3*(nPts+3*nSide3+k)+1] = Lagrange(xlist, xi, j2-1-k) * dLagrange(xlist, eta, j)      * Lagrange(xlist, mu, i);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, xi, j)      * Lagrange(xlist, eta, j+1+k)  * dLagrange(xlist, mu, i);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, xi, j+1+k)  * Lagrange(xlist, eta, j2)     * dLagrange(xlist, mu, i);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, xi, j2)     * Lagrange(xlist, eta, j2-1-k) * dLagrange(xlist, mu, i);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j)      * dLagrange(xlist, mu, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2)  * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts)+1] = Lagrange(xlist, xi, nSide/2)  * dLagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, i);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, xi, nSide/2)  * Lagrange(xlist, eta, nSide/2)  * dLagrange(xlist, mu, i);
        nPts += 1;
      }

      // --- Front face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_dshape[3*(nPts+0)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, mu, j)  * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+1)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, mu, j)  * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+2)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+3)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);

        out_dshape[3*(nPts+0)+1] = Lagrange(xlist, xi, j)  * Lagrange(xlist, mu, j)  * dLagrange(xlist, eta, i);
        out_dshape[3*(nPts+1)+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j)  * dLagrange(xlist, eta, i);
        out_dshape[3*(nPts+2)+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * dLagrange(xlist, eta, i);
        out_dshape[3*(nPts+3)+1] = Lagrange(xlist, xi, j)  * Lagrange(xlist, mu, j2) * dLagrange(xlist, eta, i);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, xi, j)  * dLagrange(xlist, mu, j)  * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, xi, j2) * dLagrange(xlist, mu, j)  * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, xi, j2) * dLagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, xi, j)  * dLagrange(xlist, mu, j2) * Lagrange(xlist, eta, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = dLagrange(xlist, xi, j+1+k)  * Lagrange(xlist, mu, j)      * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+1*nSide3+k)+0] = dLagrange(xlist, xi, j2)     * Lagrange(xlist, mu, j+1+k)  * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+2*nSide3+k)+0] = dLagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j2)     * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+3*nSide3+k)+0] = dLagrange(xlist, xi, j)      * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i);

          out_dshape[3*(nPts+0*nSide3+k)+1] = Lagrange(xlist, xi, j+1+k)  * Lagrange(xlist, mu, j)      * dLagrange(xlist, eta, i);
          out_dshape[3*(nPts+1*nSide3+k)+1] = Lagrange(xlist, xi, j2)     * Lagrange(xlist, mu, j+1+k)  * dLagrange(xlist, eta, i);
          out_dshape[3*(nPts+2*nSide3+k)+1] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j2)     * dLagrange(xlist, eta, i);
          out_dshape[3*(nPts+3*nSide3+k)+1] = Lagrange(xlist, xi, j)      * Lagrange(xlist, mu, j2-1-k) * dLagrange(xlist, eta, i);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, xi, j+1+k)  * dLagrange(xlist, mu, j)      * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, xi, j2)     * dLagrange(xlist, mu, j+1+k)  * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, xi, j2-1-k) * dLagrange(xlist, mu, j2)     * Lagrange(xlist, eta, i);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, xi, j)      * dLagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, mu, nSide/2)  * Lagrange(xlist, eta, i);
        out_dshape[3*(nPts)+1] = Lagrange(xlist, xi, nSide/2)  * Lagrange(xlist, mu, nSide/2)  * dLagrange(xlist, eta, i);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, xi, nSide/2)  * dLagrange(xlist, mu, nSide/2) * Lagrange(xlist, eta, i);
        nPts += 1;
      }

      // --- Left face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_dshape[3*(nPts+0)+0] = Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, j)  * dLagrange(xlist, xi, i);
        out_dshape[3*(nPts+1)+0] = Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, j2) * dLagrange(xlist, xi, i);
        out_dshape[3*(nPts+2)+0] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * dLagrange(xlist, xi, i);
        out_dshape[3*(nPts+3)+0] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j)  * dLagrange(xlist, xi, i);

        out_dshape[3*(nPts+0)+1] = dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, j)  * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+1)+1] = dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+2)+1] = dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+3)+1] = dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, j)  * Lagrange(xlist, xi, i);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, j)  * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, j2) * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, j)  * Lagrange(xlist, xi, i);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = Lagrange(xlist, eta, j)      * Lagrange(xlist, mu, j+1+k)  * dLagrange(xlist, xi, i);
          out_dshape[3*(nPts+1*nSide3+k)+0] = Lagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, j2)     * dLagrange(xlist, xi, i);
          out_dshape[3*(nPts+2*nSide3+k)+0] = Lagrange(xlist, eta, j2)     * Lagrange(xlist, mu, j2-1-k) * dLagrange(xlist, xi, i);
          out_dshape[3*(nPts+3*nSide3+k)+0] = Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j)      * dLagrange(xlist, xi, i);

          out_dshape[3*(nPts+0*nSide3+k)+1] = dLagrange(xlist, eta, j)      * Lagrange(xlist, mu, j+1+k)  * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+1*nSide3+k)+1] = dLagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, j2)     * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+2*nSide3+k)+1] = dLagrange(xlist, eta, j2)     * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+3*nSide3+k)+1] = dLagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j)      * Lagrange(xlist, xi, i);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, eta, j)      * dLagrange(xlist, mu, j+1+k)  * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, eta, j+1+k)  * dLagrange(xlist, mu, j2)     * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, eta, j2)     * dLagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, eta, j2-1-k) * dLagrange(xlist, mu, j)      * Lagrange(xlist, xi, i);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = Lagrange(xlist, eta, nSide/2)  * Lagrange(xlist, mu, nSide/2)  * dLagrange(xlist, xi, i);
        out_dshape[3*(nPts)+1] = dLagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2)  * Lagrange(xlist, xi, i);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, eta, nSide/2)  * dLagrange(xlist, mu, nSide/2) * Lagrange(xlist, xi, i);
        nPts += 1;
      }

      // --- Right face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_dshape[3*(nPts+0)+0] = Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, j)  * dLagrange(xlist, xi, i2);
        out_dshape[3*(nPts+1)+0] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j)  * dLagrange(xlist, xi, i2);
        out_dshape[3*(nPts+2)+0] = Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * dLagrange(xlist, xi, i2);
        out_dshape[3*(nPts+3)+0] = Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, j2) * dLagrange(xlist, xi, i2);

        out_dshape[3*(nPts+0)+1] = dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, j)  * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+1)+1] = dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, j)  * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+2)+1] = dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+3)+1] = dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, j)  * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, j)  * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, j2) * Lagrange(xlist, xi, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = Lagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, j)      * dLagrange(xlist, xi, i2);
          out_dshape[3*(nPts+1*nSide3+k)+0] = Lagrange(xlist, eta, j2)     * Lagrange(xlist, mu, j+1+k)  * dLagrange(xlist, xi, i2);
          out_dshape[3*(nPts+2*nSide3+k)+0] = Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j2)     * dLagrange(xlist, xi, i2);
          out_dshape[3*(nPts+3*nSide3+k)+0] = Lagrange(xlist, eta, j)      * Lagrange(xlist, mu, j2-1-k) * dLagrange(xlist, xi, i2);

          out_dshape[3*(nPts+0*nSide3+k)+1] = dLagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, j)      * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+1*nSide3+k)+1] = dLagrange(xlist, eta, j2)     * Lagrange(xlist, mu, j+1+k)  * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+2*nSide3+k)+1] = dLagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, j2)     * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+3*nSide3+k)+1] = dLagrange(xlist, eta, j)      * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i2);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, eta, j+1+k)  * dLagrange(xlist, mu, j)      * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, eta, j2)     * dLagrange(xlist, mu, j+1+k)  * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, eta, j2-1-k) * dLagrange(xlist, mu, j2)     * Lagrange(xlist, xi, i2);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, eta, j)      * dLagrange(xlist, mu, j2-1-k) * Lagrange(xlist, xi, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = Lagrange(xlist, eta, nSide/2)  * Lagrange(xlist, mu, nSide/2)  * dLagrange(xlist, xi, i2);
        out_dshape[3*(nPts)+1] = dLagrange(xlist, eta, nSide/2)  * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, xi, i2);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, eta, nSide/2) * dLagrange(xlist, mu, nSide/2)  * Lagrange(xlist, xi, i2);
        nPts += 1;
      }

      // --- Back face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_dshape[3*(nPts+0)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, mu, j)  * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+1)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, mu, j)  * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+2)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+3)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);

        out_dshape[3*(nPts+0)+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j)  * dLagrange(xlist, eta, i2);
        out_dshape[3*(nPts+1)+1] = Lagrange(xlist, xi, j)  * Lagrange(xlist, mu, j)  * dLagrange(xlist, eta, i2);
        out_dshape[3*(nPts+2)+1] = Lagrange(xlist, xi, j)  * Lagrange(xlist, mu, j2) * dLagrange(xlist, eta, i2);
        out_dshape[3*(nPts+3)+1] = Lagrange(xlist, xi, j2) * Lagrange(xlist, mu, j2) * dLagrange(xlist, eta, i2);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, xi, j2) * dLagrange(xlist, mu, j)  * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, xi, j)  * dLagrange(xlist, mu, j)  * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, xi, j)  * dLagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, xi, j2) * dLagrange(xlist, mu, j2) * Lagrange(xlist, eta, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = dLagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j)      * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+1*nSide3+k)+0] = dLagrange(xlist, xi, j)      * Lagrange(xlist, mu, j+1+k)  * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+2*nSide3+k)+0] = dLagrange(xlist, xi, j+1+k)  * Lagrange(xlist, mu, j2)     * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+3*nSide3+k)+0] = dLagrange(xlist, xi, j2)     * Lagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i2);

          out_dshape[3*(nPts+0*nSide3+k)+1] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, mu, j)      * dLagrange(xlist, eta, i2);
          out_dshape[3*(nPts+1*nSide3+k)+1] = Lagrange(xlist, xi, j)      * Lagrange(xlist, mu, j+1+k)  * dLagrange(xlist, eta, i2);
          out_dshape[3*(nPts+2*nSide3+k)+1] = Lagrange(xlist, xi, j+1+k)  * Lagrange(xlist, mu, j2)     * dLagrange(xlist, eta, i2);
          out_dshape[3*(nPts+3*nSide3+k)+1] = Lagrange(xlist, xi, j2)     * Lagrange(xlist, mu, j2-1-k) * dLagrange(xlist, eta, i2);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, xi, j2-1-k) * dLagrange(xlist, mu, j)      * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, xi, j)      * dLagrange(xlist, mu, j+1+k)  * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, xi, j+1+k)  * dLagrange(xlist, mu, j2)     * Lagrange(xlist, eta, i2);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, xi, j2)     * dLagrange(xlist, mu, j2-1-k) * Lagrange(xlist, eta, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, mu, nSide/2) * Lagrange(xlist, eta, i2);
        out_dshape[3*(nPts)+1] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, mu, nSide/2) * dLagrange(xlist, eta, i2);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, xi, nSide/2) * dLagrange(xlist, mu, nSide/2) * Lagrange(xlist, eta, i2);
        nPts += 1;
      }

      // --- Top face ---
      for (int j0 = 0; j0 < nLevels2; j0++) {
        // Corners
        int j = j0 + i + 1;
        int j2 = i + 1 + (nSide2-1) - j0;
        out_dshape[3*(nPts+0)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+1)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, eta, j)  * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+2)+0] = dLagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+3)+0] = dLagrange(xlist, xi, j)  * Lagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);

        out_dshape[3*(nPts+0)+1] = Lagrange(xlist, xi, j)  * dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+1)+1] = Lagrange(xlist, xi, j2) * dLagrange(xlist, eta, j)  * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+2)+1] = Lagrange(xlist, xi, j2) * dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts+3)+1] = Lagrange(xlist, xi, j)  * dLagrange(xlist, eta, j2) * Lagrange(xlist, mu, i2);

        out_dshape[3*(nPts+0)+2] = Lagrange(xlist, xi, j)  * Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+1)+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j)  * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+2)+2] = Lagrange(xlist, xi, j2) * Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, i2);
        out_dshape[3*(nPts+3)+2] = Lagrange(xlist, xi, j)  * Lagrange(xlist, eta, j2) * dLagrange(xlist, mu, i2);
        nPts += 4;

        // Edges: Bottom, right, top, left
        int nSide3 = nSide2 - 2 * (j0+1);
        for (int k = 0; k < nSide3; k++) {
          out_dshape[3*(nPts+0*nSide3+k)+0] = dLagrange(xlist, xi, j+1+k)  * Lagrange(xlist, eta, j)      * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+1*nSide3+k)+0] = dLagrange(xlist, xi, j2)     * Lagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+2*nSide3+k)+0] = dLagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j2)     * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+3*nSide3+k)+0] = dLagrange(xlist, xi, j)      * Lagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i2);

          out_dshape[3*(nPts+0*nSide3+k)+1] = Lagrange(xlist, xi, j+1+k)  * dLagrange(xlist, eta, j)      * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+1*nSide3+k)+1] = Lagrange(xlist, xi, j2)     * dLagrange(xlist, eta, j+1+k)  * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+2*nSide3+k)+1] = Lagrange(xlist, xi, j2-1-k) * dLagrange(xlist, eta, j2)     * Lagrange(xlist, mu, i2);
          out_dshape[3*(nPts+3*nSide3+k)+1] = Lagrange(xlist, xi, j)      * dLagrange(xlist, eta, j2-1-k) * Lagrange(xlist, mu, i2);

          out_dshape[3*(nPts+0*nSide3+k)+2] = Lagrange(xlist, xi, j+1+k)  * Lagrange(xlist, eta, j)      * dLagrange(xlist, mu, i2);
          out_dshape[3*(nPts+1*nSide3+k)+2] = Lagrange(xlist, xi, j2)     * Lagrange(xlist, eta, j+1+k)  * dLagrange(xlist, mu, i2);
          out_dshape[3*(nPts+2*nSide3+k)+2] = Lagrange(xlist, xi, j2-1-k) * Lagrange(xlist, eta, j2)     * dLagrange(xlist, mu, i2);
          out_dshape[3*(nPts+3*nSide3+k)+2] = Lagrange(xlist, xi, j)      * Lagrange(xlist, eta, j2-1-k) * dLagrange(xlist, mu, i2);
        }
        nPts += 4*nSide3;
      }

      // Center node for even-ordered Lagrange quads (odd value of nSide)
      if (isOdd2) {
        out_dshape[3*(nPts)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts)+1] = Lagrange(xlist, xi, nSide/2) * dLagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, i2);
        out_dshape[3*(nPts)+2] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * dLagrange(xlist, mu, i2);
        nPts += 1;
      }
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd) {
      out_dshape[3*(nNodes-1)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2);
      out_dshape[3*(nNodes-1)+1] = Lagrange(xlist, xi, nSide/2) * dLagrange(xlist, eta, nSide/2) * Lagrange(xlist, mu, nSide/2);
      out_dshape[3*(nNodes-1)+2] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2) * dLagrange(xlist, mu, nSide/2);
    }
  }
}
