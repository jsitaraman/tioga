/*!
 * \file superMesh.hpp
 * \brief Header file for SuperMesh class
 *
 * Creates a local supermesh for an element in one grid from
 * elements in another grid (See Farrell and Maddison, 2010)
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <array>
#include <map>
#include <vector>

#include "points.hpp"

struct tetra
{
  std::array<point,4> nodes; //! Positions of nodes in tet
  int donorID;               //! Donor-grid cell ID to which tetra belongs
};

struct triangle
{
  std::array<point,3> nodes; //! Positions of nodes in tri
  int donorID;               //! Donor-grid cell ID to which tetra belongs
};

class SuperMesh
{
public:

  /* ---- Default Functions ---- */

  SuperMesh(std::vector<point> &_target, std::vector<point> &_donors, int _nDonors, int _nDims);
  SuperMesh(std::vector<point> &_target, std::vector<point> &_donors, int _nDonors, int _nDims, int _rank, int _ID);

  SuperMesh();
  ~SuperMesh();

  /* ---- Member Variables ---- */

  std::vector<point> target;  //! Target cell's node positions for which to create local supermesh
  std::vector<point> donors;  //! Node positions of cells from donor grid which overlap target cell

  int nSimps;     //! Total number of simplices comprising the supermesh
  int nv_simp;    //! Number of vertices in simplex (3 or 4)
  int order;      //! Order of quadrature rule to use
  int nDonors;    //! Number of donor cells
  int nDims;      //! Dimension of the problem: 2D or 3D

  std::vector<point> faceXCs; //! Centroids of target cell faces (for clipping)
  std::vector<Vec3> normals;  //! Outward face normals for target cell (for clipping)

  std::vector<tetra> tets;    //! Tetrahedrons comprising the supermesh [for 3D]
  std::vector<triangle> tris; //! triangles comprising the supermesh [for 2D]

  /* ---- Member Functions ---- */

  void setup(std::vector<point> &_target, std::vector<point> &_donors, int _nDonors, int nDims);

  //! Using given grids and target cell, build the local supermesh
  void buildSuperMesh(void);

  //! Print the simplices of the SuperMesh to a CSV file
  void printSuperMesh(int rank, int ID);

  std::vector<point> getSuperMeshPoints(void);

  int ID, rank;

private:

  void setupMaps(void);

  void buildSuperMeshTri(void);
  void buildSuperMeshTet(void);
};

/* --- Extra Helper Functions --- */

//! Subdivide the given hexahedron into 5 tetrahedrons
std::vector<tetra> splitHexIntoTets(const point* hexNodes);

//! Subdivide the given quadrilateral into 2 triangles
std::vector<triangle> splitQuadIntoTris(const point* quadNodes);

//! Use the given face and outward normal to clip the given tet and return the new set of tets
std::vector<tetra> clipTet(tetra &tet, const point &xc, Vec3 &norm);

//! Use the given face and outward normal to clip the given triangle and return the new set of tris
std::vector<triangle> clipTri(triangle &tri, const point &xc, Vec3 &norm);

double getAreaTri(std::array<point,3> &nodes);

double getVolumeTet(std::array<point,4> &nodes);
