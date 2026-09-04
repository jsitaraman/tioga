// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)
// Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause

/*
 * Device port of the host donor containment test.
 *
 * MeshBlock::checkContainment() decides whether a query point lies inside a
 * given cell by computing the point's nodal interpolation weights and checking
 * that all of them fall in [-TOL, 1+TOL]. Every GPU search backend needs that
 * same test, on the same mesh, in the same precision, so it lives here rather
 * than in any one backend.
 *
 * Everything is double precision and mirrors math.c operation for operation,
 * so a backend that visits the same candidate cells as the host reaches the
 * same verdict on each of them.
 *
 * The mesh is described by the flat, already 0-based arrays that the backends
 * upload:
 *
 *   x            3*nnodes node coordinates
 *   conn         per-type connectivity, all types concatenated
 *   connOffset   [ntypes]   where each type starts inside conn
 *   cellStart    [ntypes+1] cumulative cell counts, so cell ids are global
 *   nvertArr     [ntypes]   vertices per cell of each type
 *   cellRes      [ncells]   cell resolution, BIGVALUE marks a fringe cell
 */

#ifndef SEARCHGPU_DEVICE_H
#define SEARCHGPU_DEVICE_H

#include "codetypes.h"

/* ================================================================== */
/* device port of the host containment math (math.c)                  */
/* ================================================================== */

/* 3x3 gaussian elimination with partial (zero-pivot only) swapping,
   bit-for-bit the same sequence of operations as solvec() in math.c */
__device__ __forceinline__
int d_solvec3(double a[3][3], double b[3])
{
  const double eps = 1e-8;
  for (int i = 0; i < 3; i++) {
    if (fabs(a[i][i]) < eps) {
      int flag = 1;
      for (int k = i + 1; k < 3 && flag; k++) {
        if (a[k][i] != 0.0) {
          flag = 0;
          for (int l = 0; l < 3; l++) {
            double t = a[k][l]; a[k][l] = a[i][l]; a[i][l] = t;
          }
          double t = b[k]; b[k] = b[i]; b[i] = t;
        }
      }
      if (flag) return 0;
    }
    for (int k = i + 1; k < 3; k++) {
      double fact = -a[k][i] / a[i][i];
      for (int j = 0; j < 3; j++) a[k][j] += fact * a[i][j];
      b[k] += fact * b[i];
    }
  }
  for (int i = 2; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i + 1; j < 3; j++) sum += a[i][j] * b[j];
    b[i] = (b[i] - sum) / a[i][i];
  }
  return 1;
}

/* trilinear-form Newton solve, mirrors newtonSolve() in math.c */
__device__ __forceinline__
void d_newtonSolve(const double f[8][3], double *u1, double *v1, double *w1)
{
  const int    itmax = 500;
  const double convergenceLimit = 1e-14;
  const double alph = 1.0;
  int    isolflag = 1;
  double u = 0.5, v = 0.5, w = 0.5;
  double rhs[3], lhs[3][3];
  int iter;

  for (iter = 0; iter < itmax; iter++) {
    double uv = u * v, vw = v * w, wu = w * u, uvw = u * v * w;

    for (int j = 0; j < 3; j++)
      rhs[j] = f[0][j] + f[1][j]*u + f[2][j]*v + f[3][j]*w
             + f[4][j]*uv + f[5][j]*vw + f[6][j]*wu + f[7][j]*uvw;

    double norm = rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2];
    if (sqrt(norm) <= convergenceLimit) break;

    for (int j = 0; j < 3; j++) {
      lhs[j][0] = f[1][j] + f[4][j]*v + f[6][j]*w + f[7][j]*vw;
      lhs[j][1] = f[2][j] + f[5][j]*w + f[4][j]*u + f[7][j]*wu;
      lhs[j][2] = f[3][j] + f[6][j]*u + f[5][j]*v + f[7][j]*uv;
    }

    isolflag = d_solvec3(lhs, rhs);
    if (isolflag == 0) break;

    u -= rhs[0]*alph;
    v -= rhs[1]*alph;
    w -= rhs[2]*alph;
  }
  if (iter == itmax || isolflag == 0) { u = 2.0; v = 0.0; w = 0.0; }
  *u1 = u; *v1 = v; *w1 = w;
}

/* mirrors computeNodalWeights() in math.c */
__device__ __forceinline__
void d_computeNodalWeights(const double xv[8][3], const double *xp,
                           double frac[8], int nvert)
{
  double f[8][3], u, v, w, omu, omv, omw, omuv;

  if (nvert == 4) {
    /* tetrahedron: solve the 3x3 barycentric system directly */
    double lhs[3][3], rhs[3];
    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) lhs[j][k] = xv[k][j] - xv[3][j];
      rhs[k] = xp[k] - xv[3][k];
    }
    if (d_solvec3(lhs, rhs)) {
      for (int k = 0; k < 3; k++) frac[k] = rhs[k];
      frac[3] = 1.0 - frac[0] - frac[1] - frac[2];
    } else {
      frac[0] = 1.0; frac[1] = frac[2] = frac[3] = 0.0;
    }
    return;
  }

  if (nvert == 5) {
    /* pyramid */
    for (int j = 0; j < 3; j++) {
      f[0][j] = xv[0][j] - xp[j];
      f[1][j] = xv[1][j] - xv[0][j];
      f[2][j] = xv[3][j] - xv[0][j];
      f[3][j] = xv[4][j] - xv[0][j];
      f[4][j] = xv[0][j] - xv[1][j] + xv[2][j] - xv[3][j];
      f[5][j] = xv[0][j] - xv[3][j];
      f[6][j] = xv[0][j] - xv[1][j];
      f[7][j] = -xv[0][j] + xv[1][j] - xv[2][j] + xv[3][j];
    }
    d_newtonSolve(f, &u, &v, &w);
    omu = 1.0 - u; omv = 1.0 - v; omw = 1.0 - w;
    frac[0] = omu*omv*omw;
    frac[1] = u*omv*omw;
    frac[2] = u*v*omw;
    frac[3] = omu*v*omw;
    frac[4] = w;
    return;
  }

  if (nvert == 6) {
    /* prism */
    for (int j = 0; j < 3; j++) {
      f[0][j] = xv[0][j] - xp[j];
      f[1][j] = xv[1][j] - xv[0][j];
      f[2][j] = xv[2][j] - xv[0][j];
      f[3][j] = xv[3][j] - xv[0][j];
      f[4][j] = 0.0;
      f[5][j] = xv[0][j] - xv[2][j] - xv[3][j] + xv[5][j];
      f[6][j] = xv[0][j] - xv[1][j] - xv[3][j] + xv[4][j];
      f[7][j] = 0.0;
    }
    d_newtonSolve(f, &u, &v, &w);
    omuv = 1.0 - u - v; omw = 1.0 - w;
    frac[0] = omuv*omw;
    frac[1] = u*omw;
    frac[2] = v*omw;
    frac[3] = omuv*w;
    frac[4] = u*w;
    frac[5] = v*w;
    return;
  }

  if (nvert == 8) {
    /* hexahedron */
    for (int j = 0; j < 3; j++) {
      f[0][j] = xv[0][j] - xp[j];
      f[1][j] = xv[1][j] - xv[0][j];
      f[2][j] = xv[3][j] - xv[0][j];
      f[3][j] = xv[4][j] - xv[0][j];
      f[4][j] = xv[0][j] - xv[1][j] + xv[2][j] - xv[3][j];
      f[5][j] = xv[0][j] - xv[3][j] + xv[7][j] - xv[4][j];
      f[6][j] = xv[0][j] - xv[1][j] + xv[5][j] - xv[4][j];
      f[7][j] = -xv[0][j] + xv[1][j] - xv[2][j] + xv[3][j]
              +  xv[4][j] - xv[5][j] + xv[6][j] - xv[7][j];
    }
    d_newtonSolve(f, &u, &v, &w);
    omu = 1.0 - u; omv = 1.0 - v; omw = 1.0 - w;
    frac[0] = omu*omv*omw;
    frac[1] = u*omv*omw;
    frac[2] = u*v*omw;
    frac[3] = omu*v*omw;
    frac[4] = omu*omv*w;
    frac[5] = u*omv*w;
    frac[6] = u*v*w;
    frac[7] = omu*v*w;
    return;
  }

  /* unsupported element: report "not contained" via an impossible weight */
  frac[0] = 2.0;
}

/* ================================================================== */
/* helpers                                                            */
/* ================================================================== */

/* cell id -> (type, index within type). ntypes is 1..4 in practice. */
__device__ __forceinline__
void d_cellType(int cell, const int *cellStart, int ntypes, int *n, int *i)
{
  int t = 0;
  for (; t < ntypes - 1; t++)
    if (cell < cellStart[t+1]) break;
  *n = t;
  *i = cell - cellStart[t];
}

__device__ __forceinline__
void d_gatherCell(double xv[8][3], int cell,
                  const double *x, const int *conn, const int *connOffset,
                  const int *cellStart, const int *nvertArr, int ntypes,
                  int *nvertOut)
{
  int n, i;
  d_cellType(cell, cellStart, ntypes, &n, &i);
  int nvert = nvertArr[n];
  const int *cv = conn + connOffset[n] + (size_t)nvert*i;
  for (int m = 0; m < nvert; m++) {
    int i3 = 3*cv[m];
    xv[m][0] = x[i3+0];
    xv[m][1] = x[i3+1];
    xv[m][2] = x[i3+2];
  }
  *nvertOut = nvert;
}

/* float conversions that round strictly outward, so a float AABB always
   contains the double AABB it came from */
__device__ __forceinline__ float f_down(double v)
{
  float f = (float)v;
  return ((double)f <= v) ? f : nextafterf(f, -INFINITY);
}
__device__ __forceinline__ float f_up(double v)
{
  float f = (float)v;
  return ((double)f >= v) ? f : nextafterf(f, INFINITY);
}

/* Exact containment, mirroring MeshBlock::checkContainment() for ihigh==0.
   Returns -1 (not in cell), 0 (accepted) or 1 (accepted, but the cell is a
   zero-resolution cell and the point sits on its boundary -- the host code
   keeps searching in that case). */
__device__ __forceinline__
int d_checkContainment(int cell, const double *xp,
                       const double *x, const int *conn, const int *connOffset,
                       const int *cellStart, const int *nvertArr, int ntypes,
                       const double *cellRes)
{
  double xv[8][3], frac[8];
  int nvert;
  d_gatherCell(xv, cell, x, conn, connOffset, cellStart, nvertArr, ntypes, &nvert);
  d_computeNodalWeights(xv, xp, frac, nvert);

  int flagged = 0;
  for (int m = 0; m < nvert; m++) {
    if ((frac[m]+TOL)*(frac[m]-1.0-TOL) > 0.0) return -1;
    if (fabs(frac[m]) < TOL && cellRes[cell] == BIGVALUE) flagged = 1;
  }
  return flagged;
}

#endif /* SEARCHGPU_DEVICE_H */
