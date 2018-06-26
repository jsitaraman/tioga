#ifndef MATH_FUNCS_H
#define MATH_FUNCS_H
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

/** box product of vectors, a,b,c */
double scalarProduct(double a[3], double b[3], double c[3]);

/**
 * Subroutine to find volumes of arbitrary polyhedra with triangles and quads
 * as faces
 * Uses vertex coordinates and provided face connectivity
 *
 * Uses Gausss divergence theorem with quad faces treated as bilinear
 * Reference Appendix B of https://drum.umd.edu/dspace/handle/1903/310
 * Sitaraman, UMD, PhD thesis 2003
 *
 * Arguments:
 * xc - coordinates of vertices of the polyhedra
 * fconn - face connectivity
 * numverts-number of vertices in each face
 * nfaces - number of faces
 * nvert - number of vertices of the polyhedra
 */
double cellVolume(double xv[8][3], int numverts[6], int fconn[24],
                int nfaces, int nvert);

void solvec(double **a,double *b,int *iflag,int n);

void newtonSolve(double f[7][3],double *u1,double *v1,double *w1);

void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);

double computeCellVolume(double xv[8][3],int nvert);


#endif // MATH_FUNCS_H
