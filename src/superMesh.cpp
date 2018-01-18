/*!
 * \file superMesh.cpp
 * \brief SuperMesh class definition
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "superMesh.hpp"

#include <fstream>
#include <set>
#include <iomanip>

#include "points.hpp"
#include "funcs.hpp"

static std::map<int,std::array<int,3>> flipTet1;
static std::map<std::array<int,2>,std::array<int,4>> flipTet2;
static std::map<int,std::array<int,2>> flipTri1;

/*
 * Tetrahedron Node-Ordering Conventions:
 *
 *            2                                     2
 *          ,/|`\                                 ,/|`\
 *        ,/  |  `\                             ,/  |  `\
 *      ,/    '.   `\                         ,6    '.   `5
 *    ,/       |     `\                     ,/       8     `\
 *  ,/         |       `\                 ,/         |       `\
 * 0-----------'.--------1               0--------4--'.--------1
 *  `\.         |      ,/                 `\.         |      ,/
 *     `\.      |    ,/                      `\.      |    ,9
 *        `\.   '. ,/                           `7.   '. ,/
 *           `\. |/                                `\. |/
 *              `3                                    `3
 */

SuperMesh::SuperMesh()
{
  setupMaps();
}

SuperMesh::SuperMesh(std::vector<point> &_target, std::vector<point>& _donors, int _nDonors, int _nDims)
{
  setupMaps();

  setup(_target,_donors,_nDonors,_nDims);
}

SuperMesh::SuperMesh(std::vector<point> &_target, std::vector<point>& _donors, int _nDonors, int _nDims, int _rank, int _ID)
{
  rank = _rank;
  ID = _ID;

  setupMaps();

  setup(_target,_donors,_nDonors,_nDims);
}

void SuperMesh::setupMaps(void)
{
  flipTri1[0] = {{1,2}};  flipTri1[1] = {{2,0}}; flipTri1[2] = {{0,1}};

  flipTet2[{{0,1}}] = {{0,1,2,3}};  flipTet2[{{0,2}}] = {{0,2,3,1}};
  flipTet2[{{0,3}}] = {{0,3,1,2}};  flipTet2[{{1,2}}] = {{1,2,0,3}};
  flipTet2[{{1,3}}] = {{1,3,2,0}};  flipTet2[{{2,3}}] = {{2,3,0,1}};


  flipTet1[0] = {{1,3,2}};  flipTet1[1] = {{0,2,3}};
  flipTet1[2] = {{0,3,1}};  flipTet1[3] = {{0,1,2}};
}

SuperMesh::~SuperMesh()
{

}

void SuperMesh::setup(std::vector<point> &_target, std::vector<point> &_donors, int _nDonors, int _nDims)
{
  target = _target;
  donors = _donors;
  nDonors = _nDonors;
  nDims = _nDims;

  nv_simp = nDims+1;

  buildSuperMesh();
}

void SuperMesh::buildSuperMesh(void)
{
  if (nDims==2)
    buildSuperMeshTri();
  else
    buildSuperMeshTet();
}

void SuperMesh::buildSuperMeshTri(void)
{
  // Step 1: Split the donor hexahedrons into triangles to prepare for clipping
  tris.resize(0);
  for (int i=0; i<nDonors; i++) {
    auto tmpTris = splitQuadIntoTris(&donors[4*i]);
    tris.insert(tris.end(),tmpTris.begin(),tmpTris.end());
  }

  // Step 2: Get the clipping planes from the target-cell faces
  normals.resize(4);

  point xc;
  for (auto &pt:target)
    xc += pt;
  xc /= target.size();

  std::vector<point> facePts = {target[0],target[1]};
  faceXCs.push_back(getCentroid(facePts));
  normals[0] = getEdgeNormal(facePts,xc);

  facePts = {target[1],target[2]};
  faceXCs.push_back(getCentroid(facePts));
  normals[1] = getEdgeNormal(facePts,xc);

  facePts = {target[2],target[3]};
  faceXCs.push_back(getCentroid(facePts));
  normals[2] = getEdgeNormal(facePts,xc);

  facePts = {target[3],target[0]};
  faceXCs.push_back(getCentroid(facePts));
  normals[3] = getEdgeNormal(facePts,xc);

  // Step 3: Use the faces to clip the tets
  for (uint i=0; i<4; i++) {
    std::vector<triangle> newTris;
    for (uint j=0; j<tris.size(); j++) {
      triangle tri = tris[j];
      auto tmpTris = clipTri(tri, faceXCs[i], normals[i]);

      newTris.insert(newTris.end(),tmpTris.begin(),tmpTris.end());
    }
    tris = newTris;
  }
}

void SuperMesh::buildSuperMeshTet(void)
{
  // Step 1: Split the donor hexahedrons into tets to prepare for clipping
  tets.resize(0);
  for (int i=0; i<nDonors; i++) {
    auto tmpTets = splitHexIntoTets(&donors[8*i]);
    tets.insert(tets.end(),tmpTets.begin(),tmpTets.end());
  }

  // Step 2: Get the clipping planes from the target-cell faces
  normals.resize(6); // Consider slightly-more-accurate approach of splitting faces into tris instead (12 faces)
  faceXCs.resize(6);

  point xc;
  for (auto &pt:target)
    xc += pt;
  xc /= target.size();

  std::vector<point> facePts = {target[0],target[1],target[2],target[3]};
  faceXCs[0] = getCentroid(facePts);
  normals[0] = getFaceNormalQuad(facePts,xc);

  facePts = {target[4],target[5],target[6],target[7]};
  faceXCs[1] = getCentroid(facePts);
  normals[1] = getFaceNormalQuad(facePts,xc);

  facePts = {target[0],target[3],target[7],target[4]};
  faceXCs[2] = getCentroid(facePts);
  normals[2] = getFaceNormalQuad(facePts,xc);

  facePts = {target[2],target[1],target[5],target[6]};
  faceXCs[3] = getCentroid(facePts);
  normals[3] = getFaceNormalQuad(facePts,xc);

  facePts = {target[0],target[1],target[5],target[4]};
  faceXCs[4] = getCentroid(facePts);
  normals[4] = getFaceNormalQuad(facePts,xc);

  facePts = {target[3],target[2],target[6],target[7]};
  faceXCs[5] = getCentroid(facePts);
  normals[5] = getFaceNormalQuad(facePts,xc);

  // Step 3: Use the faces to clip the tets
  for (uint i=0; i<6; i++) {
    std::vector<tetra> newTets;
    for (uint j=0; j<tets.size(); j++) {
      tetra tet = tets[j];
      auto tmpTets = clipTet(tet, faceXCs[i], normals[i]);
      newTets.insert(newTets.end(),tmpTets.begin(),tmpTets.end());
    }
    tets = newTets;
  }
}

void SuperMesh::printSuperMesh(int rank, int ID)
{
  std::string filename = "mesh_" + std::to_string((long long)rank) + "_" + std::to_string((long long)ID) + ".csv";
  std::ofstream mesh(filename.c_str());

  if (nDims == 2) {
    mesh << "Triangle,x,y" << std::endl;
    for (int i=0; i<tris.size(); i++) {
      point pt = tris[i].nodes[0];
      mesh << i << "," << pt.x << "," << pt.y << std::endl;
      pt = tris[i].nodes[1];
      mesh << i << "," << pt.x << "," << pt.y << std::endl;
      pt = tris[i].nodes[2];
      mesh << i << "," << pt.x << "," << pt.y << std::endl;
    }
  } else {
    mesh << "Tet,x,y,z" << std::endl;
    for (int i=0; i<tets.size(); i++) {
      point pt = tets[i].nodes[0];
      std::cout.setf(std::ios::fixed);
      std::cout.precision(8);
      mesh << i << ", " << std::right << std::setw(14) << pt.x << "," << std::right << std::setw(14) << pt.y << "," << std::right << std::setw(14) << pt.z << std::endl;
      pt = tets[i].nodes[1];
      mesh << i << ", " << std::right << std::setw(14) << pt.x << "," << std::right << std::setw(14) << pt.y << "," << std::right << std::setw(14) << pt.z << std::endl;
      pt = tets[i].nodes[2];
      mesh << i << ", " << std::right << std::setw(14) << pt.x << "," << std::right << std::setw(14) << pt.y << "," << std::right << std::setw(14) << pt.z << std::endl;
      pt = tets[i].nodes[3];
      mesh << i << ", " << std::right << std::setw(14) << pt.x << "," << std::right << std::setw(14) << pt.y << "," << std::right << std::setw(14) << pt.z << std::endl;
    }
  }
  mesh.close();
}

std::vector<point> SuperMesh::getSuperMeshPoints(void)
{
  std::vector<point> out_pts;

  if (nDims == 3)
  {
    out_pts.resize(tets.size()*4);

    for (int i = 0; i < tets.size(); i++)
      for (int j = 0; j < 4; j++)
        out_pts[4*i+j] = tets[i].nodes[j];
  }
  else
  {
    out_pts.resize(tris.size()*3);

    for (int i = 0; i < tris.size(); i++)
      for (int j = 0; j < 3; j++)
        out_pts[3*i+j] = tris[i].nodes[j];
  }

  return out_pts;
}

std::vector<tetra> splitHexIntoTets(const point* hexNodes)
{
  std::vector<tetra> newTets(5);

  short ind[5][4] = {{0,1,4,3},{2,1,6,3},{5,1,6,4},{7,3,4,6},{1,3,6,4}};

  for (short i=0; i<5; i++)
    for (short j=0; j<4; j++)
      newTets[i].nodes[j] = hexNodes[ind[i][j]];

  return newTets;
}

std::vector<tetra> clipTet(tetra &tet, const point& xc, Vec3 &norm)
{
  /* --- WARNING: Assuming a linear, planar face --- */

  std::vector<tetra> outTets;

  std::set<int> deadPts;
  for (int i=0; i<4; i++) {
    // Check each point of tetra to see which must be removed
    Vec3 dx = tet.nodes[i] - xc;
    double dot = dx*norm;
    if (dot > 0) // Point lies on cut-side of clipping plane
      deadPts.insert(i);
  }

  /*
   * Perform the clipping and subdivide the new volume into new tets
   * Only 3 cases in which the clipping can occur
   * New points are created at the intersections of the original tet's edges
   * with the clipping plane: http://geomalgorithms.com/a05-_intersect-1.html
   */
  switch (deadPts.size()) {
    case 0: {
      // No intersection
      outTets.push_back(tet);
      break;
    }

    case 1: {
      // Remove 1 point to get a prism; split prism into 3 new tets
      int kill = *(deadPts.begin());  // The point to remove

      // Get the new points by intersecting the tet's edges with the clipping plane
      // Have to be careful about orientation of final tet
      std::array<int,3> ePts = flipTet1[kill];

      // Find the intersection points
      std::array<point,3> newPts;
      for (int i=0; i<3; i++) {
        Vec3 ab = tet.nodes[ePts[i]] - tet.nodes[kill];
        Vec3 ac = xc - tet.nodes[kill];
        newPts[i] = ab*((norm*ac)/(norm*ab)) + tet.nodes[kill];
      }

      outTets.resize(3);
      outTets[0].nodes = {{tet.nodes[ePts[0]], tet.nodes[ePts[1]], newPts[0], tet.nodes[ePts[2]]}};
      outTets[1].nodes = {{tet.nodes[ePts[2]], newPts[0], newPts[2], newPts[1]}};
      outTets[2].nodes = {{tet.nodes[ePts[1]], tet.nodes[ePts[2]], newPts[1],newPts[0]}};
      break;
    }

    case 2: {
      // Tet cut in half through 4 edges; split into 3 new tets
      // Get the points we're keeping and 'killing'
      std::array<int,2> kill, keep;
      int m=0, n=0;
      for (int i=0; i<4; i++) {
        if (deadPts.count(i)) {
          kill[m] = i;
          m++;
        }
        else {
          keep[n] = i;
          n++;
        }
      }

      /* Re-orient tet (shuffle nodes) based on kept nodes so that
       * clipping becomes standardized; 'base case' is keeping {0,1}
       * One possible case for each edge being removed */
      std::array<int,4> ind = flipTet2[keep];

      // Intersect the plane with the edges to get the new points
      std::array<point,4> newPts;
      point a,b;
      Vec3 ab,ac;

      // Edge 0-3
      a = tet.nodes[ind[0]];
      b = tet.nodes[ind[3]];
      ab = b - a;
      ac = xc - a;
      newPts[0] = ab*((norm*ac)/(norm*ab)) + a;

      // Edge 1-3
      a = tet.nodes[ind[1]];
      b = tet.nodes[ind[3]];
      ab = b - a;
      ac = xc - a;
      newPts[1] = ab*((norm*ac)/(norm*ab)) + a;

      // Edge 1-2
      a = tet.nodes[ind[1]];
      b = tet.nodes[ind[2]];
      ab = b - a;
      ac = xc - a;
      newPts[2] = ab*((norm*ac)/(norm*ab)) + a;

      // Edge 0-2
      a = tet.nodes[ind[0]];
      b = tet.nodes[ind[2]];
      ab = b - a;
      ac = xc - a;
      newPts[3] = ab*((norm*ac)/(norm*ab)) + a;

      // Setup the new tets
      outTets.resize(3);
      outTets[0].nodes = {{tet.nodes[ind[1]],newPts[0],newPts[3],tet.nodes[ind[0]]}};
      outTets[1].nodes = {{newPts[0],newPts[3],newPts[1],tet.nodes[ind[1]]}};
      outTets[2].nodes = {{newPts[1],newPts[3],newPts[2],tet.nodes[ind[1]]}};
      break;
    }

    case 3: {
      // The opposite of case 1; new tet is one corner of original tet
      int keep = -1;
      for (int i=0; i<4; i++) {
        if (!deadPts.count(i)) {
          keep = i;
          break;
        }
      }

      // Get the new points by intersecting the tet's edges with the clipping plane
      // Have to be careful about orientation of final tet, so map to a 'standard' orientation

      std::array<int,3> ePts = flipTet1[keep];

      // Setup outgoing tet; node 3 is the 'kept' node
      // Find the intersection points
      outTets.resize(1);
      outTets[0].nodes[3] = tet.nodes[keep];
      for (int i=0; i<3; i++) {
        Vec3 ab = tet.nodes[ePts[i]] - tet.nodes[keep];
        Vec3 ac = xc - tet.nodes[keep];
        outTets[0].nodes[i] = ab*((norm*ac)/(norm*ab)) + tet.nodes[keep];
      }
      break;
    }

    case 4: {
      // Entire tet is beyond clipping face
      break;
    }
  }

  return outTets;
}


std::vector<triangle> clipTri(triangle &tri, const point& xc, Vec3 &norm)
{
  /* --- WARNING: Assuming a linear edge --- */

  std::vector<triangle> outTris;

  std::set<int> deadPts;
  for (int i=0; i<3; i++) {
    // Check each point of triangle to see which must be removed
    Vec3 dx = tri.nodes[i] - xc;
    double dot = dx*norm;
    if (dot > 0) // Point lies on cut-side of clipping plane
      deadPts.insert(i);
  }

  /*
   * Perform the clipping and subdivide the new volume into new tris
   */
  switch (deadPts.size()) {
    case 0: {
      // No intersection.
      outTris.push_back(tri);
      break;
    }
    case 1: {
      // Removing one corner of tri
      int kill = *(deadPts.begin());  // The point to remove

      // Get the points being kept; map to a 'standard' triangle
      std::array<int,2> ePts = flipTri1[kill];

      // Find the cutting-plane intersection points
      Vec3 ab = tri.nodes[ePts[0]] - tri.nodes[kill];
      Vec3 ac = xc - tri.nodes[kill];
      point newPt1 = ab*(norm*ac)/(norm*ab) + tri.nodes[kill];

      ab = tri.nodes[ePts[1]] - tri.nodes[kill];
      point newPt2 = ab*(norm*ac)/(norm*ab) + tri.nodes[kill];

      outTris.resize(2);
      outTris[0].nodes = {{tri.nodes[ePts[0]],tri.nodes[ePts[1]],newPt1}};
      outTris[1].nodes = {{tri.nodes[ePts[1]],newPt2,newPt1}};
      break;
    }
    case 2: {
      // Keeping one corner of tri
      int keep = -1;
      for (int i=0; i<3; i++) {
        if (!deadPts.count(i)) {
          keep = i;
          break;
        }
      }

      std::array<int,2> ePts = flipTri1[keep];

      // Setup outgoing tri; node 2 is the 'kept' node
      // Find the intersection points
      outTris.resize(1);
      outTris[0].nodes[2] = tri.nodes[keep];
      for (int i=0; i<2; i++) {
        Vec3 ab = tri.nodes[ePts[i]] - tri.nodes[keep];
        Vec3 ac = xc - tri.nodes[keep];
        outTris[0].nodes[i] = ab*((norm*ac)/(norm*ab)) + tri.nodes[keep];
      }
      break;
    }
    case 3: {
      // Entire tri is beyond clipping face
      break;
    }
  }

  return outTris;
}


std::vector<triangle> splitQuadIntoTris(const point* quadNodes)
{
  std::vector<triangle> newTris(2);

  newTris[0].nodes = {{quadNodes[0],quadNodes[1],quadNodes[3]}};
  newTris[1].nodes = {{quadNodes[1],quadNodes[2],quadNodes[3]}};

  return newTris;
}

double getAreaTri(std::array<point,3> &nodes)
{
  Vec3 dx1 = nodes[1] - nodes[0];
  Vec3 dx2 = nodes[2] - nodes[0];
  Vec3 cross = dx1.cross(dx2);

  return 1./2.*std::abs(cross.z);
}

double getVolumeTet(std::array<point,4> &nodes)
{
  Vec3 dx1 = nodes[1] - nodes[0];
  Vec3 dx2 = nodes[2] - nodes[0];
  Vec3 dx3 = nodes[3] - nodes[0];
  Vec3 cross = dx1.cross(dx2);

  return 1./6.*std::abs(dx3*cross);
}
