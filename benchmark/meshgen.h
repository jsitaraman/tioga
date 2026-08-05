// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

/* Synthetic overset mesh block + query cloud generation, shared by the
   single-block (bench_search.C) and many-blocks-per-rank
   (bench_multiblock.C) benchmarks. Header-only, no state. */

#ifndef TIOGA_BENCH_MESHGEN_H
#define TIOGA_BENCH_MESHGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <vector>
#include <algorithm>

#include "codetypes.h"
#include "MeshBlock.h"


static double wtime(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC,&ts);
  return (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec;
}

/* deterministic, portable RNG so runs are reproducible across machines */
struct Rng
{
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed ? seed : 0x9e3779b97f4a7c15ull) {}
  uint64_t next(void)
  {
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
  }
  double uniform(void) { return (double)(next() >> 11) * (1.0/9007199254740992.0); }
};

/* ------------------------------------------------------------------ */
/* synthetic mesh                                                     */
/* ------------------------------------------------------------------ */

enum ElType { EL_HEX = 0, EL_PRISM, EL_TET, EL_MIXED };

static const char *elTypeName(int t)
{
  switch (t) {
  case EL_HEX:   return "hex";
  case EL_PRISM: return "prism";
  case EL_TET:   return "tet";
  default:       return "mixed";
  }
}

struct Mesh
{
  int      nnodes;
  int      ncells;
  int      ntypes;
  double  *x;         /* 3*nnodes, owned here                         */
  int     *iblank;    /* nnodes                                       */
  int      nv[2];
  int      nc[2];
  int     *conn[2];   /* per type, 1-based                            */
  int    **vconn;
  int      nobc;
  int     *obcnode;
  double   lo[3], hi[3];
};

/* Kuhn subdivision of a hex into 6 tets. Every cell uses the same paths from
   corner 0 to corner 6, so the diagonals agree across shared faces and the
   resulting tet mesh is conforming. */
static const int KUHN[6][4] = {
  {0,1,2,6}, {0,1,5,6}, {0,3,2,6}, {0,3,7,6}, {0,4,5,6}, {0,4,7,6}
};
/* the two prisms a hex splits into (triangle 0-1-2 under 4-5-6, etc.) */
static const int PRISM[2][6] = {
  {0,1,2,4,5,6}, {0,2,3,4,6,7}
};

/* origin: optional lower corner, so several blocks can be placed side by side
   in one domain. NULL means the unit cube. */
static void buildMesh(Mesh &m, int nx, int ny, int nz, double perturb,
                      int eltype, const double *origin = NULL)
{
  const double L = 1.0;
  const double ox = origin ? origin[0] : 0.0;
  const double oy = origin ? origin[1] : 0.0;
  const double oz = origin ? origin[2] : 0.0;
  const int    npx = nx+1, npy = ny+1, npz = nz+1;

  const int nhexcells = nx*ny*nz;
  m.nnodes = npx*npy*npz;
  m.x      = (double *)malloc(sizeof(double)*3*m.nnodes);
  m.iblank = (int *)malloc(sizeof(int)*m.nnodes);

  const double dx = L/nx, dy = L/ny, dz = L/nz;

  for (int k = 0; k < npz; k++)
    for (int j = 0; j < npy; j++)
      for (int i = 0; i < npx; i++) {
        int    n  = (k*npy + j)*npx + i;
        double X  = i*dx, Y = j*dy, Z = k*dz;
        /* smooth interior displacement; zero on the outer faces so the
           block stays a cube */
        double w  = sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
        m.x[3*n+0] = ox + X + perturb*dx*w*sin(3.0*M_PI*Y);
        m.x[3*n+1] = oy + Y + perturb*dy*w*sin(3.0*M_PI*Z);
        m.x[3*n+2] = oz + Z + perturb*dz*w*sin(3.0*M_PI*X);
        m.iblank[n] = 1;
      }

  m.lo[0] = ox; m.lo[1] = oy; m.lo[2] = oz;
  for (int j = 0; j < 3; j++) m.hi[j] = m.lo[j] + L;

  /* Every element is carved out of the same logical hex, so all four element
     type mixes describe the same geometry and are conforming.
     Hex node ordering expected by computeNodalWeights():
     bottom face 0-1-2-3 (ccw), top face 4-5-6-7 with 4 above 0 */
  std::vector<int> ec[2];
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
        int h[8];
        #define NODE(a,b,d) (((d)*npy + (b))*npx + (a) + BASE)
        h[0] = NODE(i  ,j  ,k  );
        h[1] = NODE(i+1,j  ,k  );
        h[2] = NODE(i+1,j+1,k  );
        h[3] = NODE(i  ,j+1,k  );
        h[4] = NODE(i  ,j  ,k+1);
        h[5] = NODE(i+1,j  ,k+1);
        h[6] = NODE(i+1,j+1,k+1);
        h[7] = NODE(i  ,j+1,k+1);
        #undef NODE

        int kind = eltype;
        if (eltype == EL_MIXED) kind = ((i+j+k) & 1) ? EL_PRISM : EL_HEX;

        if (kind == EL_HEX) {
          for (int m2 = 0; m2 < 8; m2++) ec[0].push_back(h[m2]);
        } else if (kind == EL_PRISM) {
          int t = (eltype == EL_MIXED) ? 1 : 0;
          for (int p = 0; p < 2; p++)
            for (int m2 = 0; m2 < 6; m2++) ec[t].push_back(h[PRISM[p][m2]]);
        } else {
          for (int p = 0; p < 6; p++)
            for (int m2 = 0; m2 < 4; m2++) ec[0].push_back(h[KUHN[p][m2]]);
        }
      }

  if (eltype == EL_MIXED) {
    m.ntypes = 2;
    m.nv[0] = 8; m.nv[1] = 6;
  } else {
    m.ntypes = 1;
    m.nv[0] = (eltype == EL_HEX) ? 8 : (eltype == EL_PRISM ? 6 : 4);
    m.nv[1] = 0;
  }
  m.ncells = 0;
  m.vconn = (int **)malloc(sizeof(int *)*m.ntypes);
  for (int t = 0; t < m.ntypes; t++) {
    m.nc[t] = (int)ec[t].size()/m.nv[t];
    m.conn[t] = (int *)malloc(sizeof(int)*(ec[t].size() ? ec[t].size() : 1));
    memcpy(m.conn[t], ec[t].data(), sizeof(int)*ec[t].size());
    m.vconn[t] = m.conn[t];
    m.ncells += m.nc[t];
  }
  (void)nhexcells;

  /* outer boundary nodes act as overset boundary nodes */
  std::vector<int> obc;
  for (int k = 0; k < npz; k++)
    for (int j = 0; j < npy; j++)
      for (int i = 0; i < npx; i++)
        if (i == 0 || i == npx-1 || j == 0 || j == npy-1 || k == 0 || k == npz-1)
          obc.push_back((k*npy + j)*npx + i + BASE);
  m.nobc = (int)obc.size();
  m.obcnode = (int *)malloc(sizeof(int)*(m.nobc > 0 ? m.nobc : 1));
  memcpy(m.obcnode, obc.data(), sizeof(int)*m.nobc);
}

static void freeMesh(Mesh &m)
{
  free(m.x); free(m.iblank); free(m.vconn); free(m.obcnode);
  for (int t = 0; t < m.ntypes; t++) free(m.conn[t]);
}

/* ------------------------------------------------------------------ */
/* query cloud                                                        */
/* ------------------------------------------------------------------ */

/* qfrac: edge length of the query box as a fraction of the block, centred
   in the block. qfrac=1 means the queries cover the whole block (nothing
   for the OBB pre-filter to reject); small qfrac is the regime the host
   pre-filter was written for. */
static void buildQueries(const Mesh &m, int nq, double qfrac, double dupFrac,
                         uint64_t seed, std::vector<double> &xs)
{
  Rng rng(seed);
  xs.resize(3*nq);
  double c[3], h[3];
  for (int j = 0; j < 3; j++) {
    c[j] = 0.5*(m.lo[j] + m.hi[j]);
    h[j] = 0.5*(m.hi[j] - m.lo[j])*qfrac;
  }
  int nuniq = 0;
  for (int i = 0; i < nq; i++) {
    if (nuniq > 0 && rng.uniform() < dupFrac) {
      int src = (int)(rng.uniform()*nuniq);
      if (src >= nuniq) src = nuniq-1;
      for (int j = 0; j < 3; j++) xs[3*i+j] = xs[3*src+j];
    } else {
      for (int j = 0; j < 3; j++)
        xs[3*i+j] = c[j] + (2.0*rng.uniform()-1.0)*h[j];
      nuniq = i+1;
    }
  }
}

/* ------------------------------------------------------------------ */
/* independent containment check, used to classify CPU/GPU mismatches  */
/* ------------------------------------------------------------------ */

extern "C" {
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}

/* Does cell `cell` (global cell numbering, types concatenated in order)
   contain xp, by the same criterion checkContainment() uses? */
static bool donorIsValid(const Mesh &m, int cell, const double *xp)
{
  if (cell < 0) return false;
  int t = 0, base = 0;
  for (; t < m.ntypes; t++) {
    if (cell < base + m.nc[t]) break;
    base += m.nc[t];
  }
  if (t == m.ntypes) return false;
  int nvert = m.nv[t];
  const int *cv = m.conn[t] + (size_t)nvert*(cell - base);
  double xv[8][3], frac[8];
  for (int i = 0; i < nvert; i++)
    for (int j = 0; j < 3; j++) xv[i][j] = m.x[3*(cv[i]-BASE)+j];
  computeNodalWeights(xv, const_cast<double *>(xp), frac, nvert);
  for (int i = 0; i < nvert; i++)
    if ((frac[i]+TOL)*(frac[i]-1.0-TOL) > 0) return false;
  return true;
}

#endif /* TIOGA_BENCH_MESHGEN_H */
