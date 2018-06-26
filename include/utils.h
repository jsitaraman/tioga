#ifndef UTILS_H
#define UTILS_H
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

/** Get etype index from global cell ID */
int get_cell_type(int* nc, int ntypes, int ic_in);

void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);

/** Check if a point is inside the provided hole map (return hole status) */
int checkHoleMap(double *x,int *nx,int *sam,double *extents);

/**
 fill a given hole map using iterative
 flood fill from outside the marked boundary.
 boundary is marked by "2"
*/
void fillHoleMap(int *holeMap, int ix[3],int isym);

/** Determine whether two OBB's intersect */
int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
           double vB[3][3],double xB[3],double dxB[3]);

void writebbox(OBB *obb,int bid);
void writePoints(double *x,int nsearch,int bid);

/**
 * Create a unique hash for list of coordinates with duplicates in
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes(double *x,double *rtag,int *itag,int nn);

#endif // UTILS_H
