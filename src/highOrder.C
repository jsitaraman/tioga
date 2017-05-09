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
#include "MeshBlock.h"

#include <omp.h>

#include <algorithm>
#include <set>
#include <vector>

#ifdef _GPU
#include "cuda_funcs.h"
#include "highOrder_kernels.h"
#endif

#define ROW 0
#define COLUMN 1
#define NFRAC 1331

#define FRINGE -1
#define HOLE 0
#define NORMAL 1

#define HOLE_CUT 0
#define FIELD_CUT 1

extern "C" 
{
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}

using namespace tg_funcs;

void MeshBlock::extraConn(void)
{
  // Cell-to-cell connectivity
  int nface = nDims*2; // Quad or hex only, same as rest of artBnd stuff
  c2c.assign(ncells*nface,-1);

  for (int ic = 0; ic < ncells; ic++)
  {
    for (int j = 0; j < nface; j++)
    {
      int ff = c2f[nface*ic+j];
      int ic1 = f2c[2*ff];
      int ic2 = f2c[2*ff+1];

      if (ic1 == ic)
        c2c[nface*ic+j] = ic2;
      else
        c2c[nface*ic+j] = ic1;
    }
  }

  // List of all solid wall-boundary faces
  std::set<int> cut_faces;
  for (int ff = 0; ff < nfaces; ff++)
  {
    // All nodes of face should be wall nodes or not, so just look at node 0
    int count = 0;
    for (int j = 0; j < nfv[0]; j++)
    {
      int iv = fconn[0][ff*nfv[0] + j];
      for (int k = 0; k < nwbc; k++)
        if (iv == wbcnode[k])
          count++;
    }

    if (count == nfv[0])
      cut_faces.insert(ff);
  }

  nCutHole = cut_faces.size();
  cutFacesW.reserve(nCutHole);
  for (auto &ff : cut_faces)
    cutFacesW.push_back(ff);

  // List of all overset-boundary faces
  cut_faces.clear();
  for (int ff = 0; ff < nfaces; ff++)
  {
    // All nodes of face should be wall nodes or not, so just look at node 0
    int iv = fconn[0][ff*nfv[0] + 0];
    for (int k = 0; k < nobc; k++)
      if (iv == obcnode[k])
      {
        cut_faces.insert(ff);
        break;
      }
  }

  nCutFringe = cut_faces.size();
  cutFacesO.reserve(nCutFringe);
  for (auto &ff : cut_faces)
    cutFacesO.push_back(ff);
}

void MeshBlock::setGridVelocity(double *grid_vel)
{
  vg = grid_vel;

  // Setup temporary buffers for performing iteration-level hole cutting
  ibc_2.assign(ncells, NORMAL);
  x2.resize(3*nnodes);
}

void MeshBlock::calcNextGrid(double dt)
{
  for (int i = 0; i < nnodes; i++)
    for (int d = 0; d < 3; d++)
      x2[3*i+d] = x[3*i+d] + vg[3*i+d] * dt;

  xtmp = x;
  x = x2.data();

  ibc_tmp = iblank_cell;
  iblank_cell = ibc_2.data();
}

void MeshBlock::swapPointers(void)
{
  xtmp = x;
//  x = x2.data();

  ibc_0.assign(iblank_cell,iblank_cell+ncells);

  ibc_tmp = iblank_cell;
  iblank_cell = ibc_2.data();
}

void MeshBlock::resetCurrentGrid(void)
{
  x = xtmp;
  iblank_cell = ibc_tmp;
}

int MeshBlock::getIterIblanks(void)
{
  unblanks.clear();

  // ibc_0 is the current blanking, iblank_cell is the new blanking,
  // and ibc_2 is the estimate for the end of the time step
  for (int ic = 0; ic < ncells; ic++)
  {
    if ((ibc_0[ic] != NORMAL) && (ibc_2[ic] == NORMAL || iblank_cell[ic] == NORMAL))
    {
      unblanks.insert(ic);
      iblank_cell[ic] = FRINGE;
    }
  }

  return unblanks.size();
}

void MeshBlock::clearUnblanks(void)
{
  for (auto & ic : unblanks)
    iblank_cell[ic] = NORMAL;

  nreceptorCells = 0;

  unblanks.clear();
}

// void MeshBlock::directCut(int nGroups, int* groupIDs, int* cutType, int* nGf,
//     int** cutFaces, int gridID, int nGrids, MPI_Comm &scomm)
// {
  /*
   * Variables to add (incl. callbacks for):
   * int nGroups     (# of cutting groups present on THIS RANK)
   * int* groupIDs   (Cutting-group IDs on this rank [default to grid ID])
   * int* cutType    (group type - hole or field cutting)
   * int* nGf        (# of faces per group)
   * int** cutFaces  (list of faces for each cut group)
   * int* f2v        (if not already)
   *
   * std::vector<std::vector<double>> cutBbox  (bounding box of each group)
   */
void MeshBlock::getCutGroupBoxes(std::vector<double> &cutBox, std::vector<std::vector<double>> &faceBox, int nGroups_glob)
{
  cutBox.resize(6*nGroups_glob);
  //faceBox.resize(nGroups);
  faceBox.resize(nGroups_glob);

  // Generate bounding boxes for each cutting group and face
  for (int g = 0; g < nGroups; g++)
  {
    int G = groupIDs[g];

    // Initialize this group's bounding box
    for (int d = 0; d < nDims; d++)
    {
      cutBox[6*G+d]   = BIGVALUE;
      cutBox[6*G+d+3] = -BIGVALUE;
    }

    // Loop over all faces of this cutting group
    faceBox[G].resize(nGf[g]);
    for (int i = 0; i < nGf[g]; i++)
    {
      // Initialize this face's bounding box
      for (int d = 0; d < nDims; d++)
      {
        faceBox[G][6*i+d]   = BIGVALUE;
        faceBox[G][6*i+d+3] = -BIGVALUE;
      }

      // Loop over vertices of this face
      int ff = cutFaces[g][i];
      for (int n = 0; n < nfv[0]; n++)
      {
        int iv = fconn[0][ff*nfv[0] + n];
        for (int d = 0; d < nDims; d++)
        {
          cutBox[6*G+d]   = std::min(cutBox[6*G+d], x[3*iv+d]);
          cutBox[6*G+d+3] = std::max(cutBox[6*G+d+3], x[3*iv+d]);
          faceBox[g][6*i+d]   = std::min(faceBox[g][6*i+d], x[3*iv+d]);
          faceBox[g][6*i+d+3] = std::max(faceBox[g][6*i+d+3], x[3*iv+d]);
        }
      }
    }
  }
}

//void MeshBlock::getCutGroupFaces(std::vector<std::vector<double>> &faceNodes, int nGroups_glob)
//{
//  int nvert = nfv[0];

//  // Generate bounding boxes for each cutting group and face
//  for (int g = 0; g < nGroups; g++)
//  {
//    int G = groupIDs[g];

//    // Loop over all faces of this cutting group
//    faceNodes[G].resize(3*nvert*nGf[g]);
//    for (int i = 0; i < nGf[g]; i++)
//    {
//      // Loop over vertices of this face
//      int ff = cutFaces[g][i];
//      for (int n = 0; n < nvert; n++)
//      {
//        int iv = fconn[0][ff*nvert + n];
//        for (int d = 0; d < nDims; d++)
//          faceNodes[G][3*(nvert*i+n)+d] = x[3*iv+d];
//      }
//    }
//  }
//}

int MeshBlock::getCuttingFaces(std::vector<double> &faceNodesW, std::vector<double> &faceNodesO)
{
  int nvert = nfv[0];

  // Generate bounding boxes for each cutting group and face
  faceNodesW.resize(nDims*nvert*nCutHole);
  for (int i = 0; i < nCutHole; i++)
  {
    // Loop over vertices of this face
    int ff = cutFacesW[i];
    for (int n = 0; n < nvert; n++)
    {
      int iv = fconn[0][ff*nvert + n];
      for (int d = 0; d < nDims; d++)
        faceNodesW[nDims*(nvert*i+n)+d] = x[3*iv+d];
    }
  }

  faceNodesO.resize(nDims*nvert*nCutFringe);
  for (int i = 0; i < nCutFringe; i++)
  {
    // Loop over vertices of this face
    int ff = cutFacesO[i];
    for (int n = 0; n < nvert; n++)
    {
      int iv = fconn[0][ff*nvert + n];
      for (int d = 0; d < nDims; d++)
        faceNodesO[nDims*(nvert*i+n)+d] = x[3*iv+d];
    }
  }

  return nvert;
}


void MeshBlock::getDirectCutCells(std::vector<std::unordered_set<int>> &cellList, std::vector<double> &cutBox_global, int nGroups_glob)
{
  /* Use ADT Search to find all cells [on this rank] which intersect with each
   * cutting group [which is not on this rank]
   * NOTE: Only 'background' grids (defined as grids with no cutting groups)
   * will be cut by 'field'-type groups */
  for (int G = 0; G < nGroups_glob; G++)
  {
    if (myGroups.count(G)) continue;
    if (nGroups > 0 && cutType_glob[G] == FIELD_CUT) continue;
    adt->searchADT_box(elementList.data(),cellList[G],&cutBox_global[6*G]);
  }
}

//void MeshBlock::directCut(std::vector<double> &faceBox, int nCutFaces, std::vector<int> &cutFlag)
//{
//  cutFlag.assign(ncells,0);

//  std::unordered_set<int> cellList;
//  for (int ff = 0; ff < nCutFaces; ff++)
//  {
//    cellList.clear();

//    // Find all cells that the cutting face might pass through
//    adt->searchADT_box(elementList, cellList, &faceBox[6*ff]);

//    // Check each cell to determine if the face intersects it
//    for (auto &ic : cellList)
//    {
//      //getRefLocNewton() // Galbratih says to use Nelder-Mead though...
//      /*if (intersectionCheck())
//      {
//        cutFlag[ic] = 1;
//      }*/
//    }
//  }
//}

//#define FILL_LOOP
#define FILL_QUEUE

void MeshBlock::directCut(std::vector<double> &cutFaces, int nCut, int nvertf,
    CutMap &cutMap, int cutType)
{
  /* Peform hole cutting with solid or overset boundary faces
   *
   * Cells intersecting cut faces are immediately flagged, and normal/hole
   * status is filled by advancing outwards from flagged cells
   *
   * If two (or more) cutting faces are (close to) the same distance away from a
   * cell, the normals of the faces are averged and used to determine a more
   * accurate cut flag
   *
   * See Galbraith's thesis, Section 4.2
   */
  if (nCut == 0)
  {
    cutMap.flag.assign(ncells,DC_NORMAL);
    return;
  }
  else
    cutMap.flag.assign(ncells,DC_UNASSIGNED);

  std::unordered_set<int> cellList, checked_cells;
  std::set<int> paintQueue;

  double bbox[6];
  std::vector<double> xv(nv[0]*nDims);
  int stride = nDims*nvertf;
  int nvert = nv[0];
  int nface = nDims*2;
  double tol = 1e-6;

  for (int ff = 0; ff < nCut; ff++)
  {
    checked_cells.clear();
    cellList.clear();

    if (rrot)
      getBoundingBox(&cutFaces[ff*stride], nvertf, nDims, bbox, Rmat.data());
    else
      getBoundingBox(&cutFaces[ff*stride], nvertf, nDims, bbox);

    // Find all cells that the cutting face might pass through
    adt->searchADT_box(elementList.data(), cellList, bbox);

    // Check each cell to determine if the face intersects it
    for (auto ic : cellList)
    {
      if (cutMap.flag[ic] != DC_CUT && !checked_cells.count(ic)) // If not already found to be cut
      {
        checked_cells.insert(ic);

        // Load up the cell nodes into an array
        for (int i = 0; i < nvert; i++)
          for (int d = 0; d < nDims; d++)
            xv[nDims*i+d] = x[nDims*vconn[0][ic*nvert+i]+d];


        /// DEBUGGING
//        double xc_check[3] = {.069, .027, 0.0};
        double xc_check[3] = {.14, -.18, 0.0};
        double xc[3] = {0,0,0};
        double dr = 0;
        getCentroid(xv.data(), nvert, 3, xc);
        for (int d = 0; d < nDims; d++)
          dr += (xc[d]-xc_check[d])*(xc[d]-xc_check[d]);
        if (std::sqrt(dr) < .05)
          printf("Cell ID to use for debugging: %d\n",ic); // 552


        // Find distance from face to cell
        Vec3 vec = intersectionCheck2(&cutFaces[ff*stride], nvertf, xv.data(), nvert, nDims);
//        double dist1 = intersectionCheck2(&cutFaces[ff*stride], nvertf, xv.data(), nvert, nDims);
        double dist1 = vec.norm();
        if (ic == 552 && (ff == 128 || ff == 129)) // Faces 128, 129
        {
          for (int i = 0; i < nvertf; i++)
            printf("Face %d: Node %d: %f %f %f\n",ff,i,cutFaces[ff*stride+i*nDims+0],
                cutFaces[ff*stride+i*nDims+1],cutFaces[ff*stride+i*nDims+2]);
          printf("\n");
          printf("IC %d: ff %d: vec.norm() = %.4e\n",ic,ff,dist1);
        }

        if (vec.norm() == 0.) // They intersect
        {
          cutMap.flag[ic] = DC_CUT;
          paintQueue.insert(ic);
        }
      }

      if (cutMap.flag[ic] == DC_CUT)
      {
        // Check all neighboring elements and assign status
        for (int j = 0; j < nface; j++)
        {
          int ic2 = c2c[nface*ic+j];
          if (checked_cells.count(ic2)) continue;
          checked_cells.insert(ic2);

          if (ic2 >= 0 && cutMap.flag[ic2] != DC_CUT)
          {
            for (int i = 0; i < nvert; i++)
              for (int d = 0; d < nDims; d++)
                xv[nDims*i+d] = x[nDims*vconn[0][ic2*nvert+i]+d];

            /// DEBUGGING
            double xc_check[3] = {.14, -.18, 0.0};
            double xc[3] = {0,0,0};
            double dr = 0;
            getCentroid(xv.data(), nvert, 3, xc);
            for (int d = 0; d < nDims; d++)
              dr += (xc[d]-xc_check[d])*(xc[d]-xc_check[d]);
            if (std::sqrt(dr) < .05)
              printf("Cell ID-2 to use for debugging: %d\n",ic); // 552

            Vec3 vec = intersectionCheck2(&cutFaces[ff*stride], nvertf, xv.data(), nvert, nDims);
            double dist = vec.norm();
//            double ptc[3], ptf[3];
//            getCentroid(&cutFaces[ff*stride], nvertf, 3, ptf);
//            getCentroid(xv.data(), nvert, 3, ptc);
//            vec = point(ptc) - point(ptf);
            if (ic2 == 552 && (ff == 128 || ff == 129)) // Faces 128, 129
            {
              for (int i = 0; i < nvertf; i++)
                printf("Face %d: Node %d: %f %f %f\n",ff,i,cutFaces[ff*stride+i*nDims+0],
                    cutFaces[ff*stride+i*nDims+1],cutFaces[ff*stride+i*nDims+2]);
              printf("\n");
              printf("IC2 %d: ff %d: vec.norm() = %.4e\n",ic2,ff,dist);
            }

            if (dist == 0)
            {
              cutMap.flag[ic2] = DC_CUT;
            }
            else if (cutMap.flag[ic2] == DC_UNASSIGNED || dist < (cutMap.dist[ic2] - tol))
            {
              // Unflagged cell, or have a closer face to use
              cutMap.dist[ic2] = dist;
              cutMap.norm[ic2] = faceNormal(&cutFaces[ff*stride], nDims);
              cutMap.nMin[ic2] = 1;

              if (cutMap.norm[ic2] * vec < 0) /// TODO: decide on standard orientation
                cutMap.flag[ic2] = DC_HOLE; // outwards normal = inside hole
              else
                cutMap.flag[ic2] = DC_NORMAL;
            }
            else if (std::abs(dist-cutMap.dist[ic2]) < tol)
            {
              // Same dist. to two faces; avg. their normals to decide
              Vec3 norm = faceNormal(&cutFaces[ff*stride], nDims);

              int N = cutMap.nMin[ic2];
              for (int d = 0; d < nDims; d++)
                cutMap.norm[ic2][d] = (N*cutMap.norm[ic2][d] + norm[d]) / (N + 1.);
              cutMap.nMin[ic2]++;

              if (cutMap.norm[ic2] * vec < 0) /// TODO: decide on standard orientation
                cutMap.flag[ic2] = DC_HOLE;
              else
                cutMap.flag[ic2] = DC_NORMAL;
            }

            paintQueue.insert(ic2);
          }

        }
      }
    }

  }

  // Now paint-fill the remainder of the grid based upon the cutting boundary

  // ----- Linear loop version -----
#ifdef FILL_LOOP
  int nUnassigned = 0;
  for (auto &flag : cutMap.flag)
    if (flag == DC_UNASSIGNED)
      nUnassigned++;

  std::vector<int> mark(ncells);

  while (nUnassigned > 0)
  {
    for (int ic = 0; ic < ncells; ic++)
    {
      if (cutMap.flag[ic] == DC_UNASSIGNED or mark[ic])
        continue;

      if (cutMap.flag[ic] == DC_CUT)
      {
        if (cutType == 1) // Solid-boundary surface
          cutMap.flag[ic] = DC_HOLE;
        else              // Overset-boundary surface for background grid
          cutMap.flag[ic] = DC_NORMAL;
      }

      for (int j = 0; j < nface; j++)
      {
        int ic2 = c2c[nface*ic+j];
        if (ic2 >= 0 and cutMap.flag[ic2] == DC_UNASSIGNED)
        {
          cutMap.flag[ic2] = cutMap.flag[ic];
          mark[ic2] = 0;
          nUnassigned--;
        }
      }

      mark[ic] = 1;
    }
  }
#endif

  // ----- Queue version -----
#ifdef FILL_QUEUE
//  for (int ic = 0; ic < ncells; ic++)
//  {
//    if (cutMap.flag[ic] == DC_CUT)
//      cutMap.flag[ic] = DC_HOLE;
//    else
//      cutMap.flag[ic] = DC_NORMAL;
//  } /// DEBUGGING

  int nstack = 0;
  std::vector<int> icstack(ncells,-1);
  for (auto ic : paintQueue)
    icstack[nstack++] = ic; /// TODO: remove set-based version?

  while (nstack > 0)
  {
    nstack--;
    int ic = icstack[nstack];

    if (ic < 0) continue;

    if (cutMap.flag[ic] == DC_CUT)
    {
      if (cutType == 1) // Solid-boundary surface
        cutMap.flag[ic] = DC_HOLE;
      else              // Overset-boundary surface for background grid
        cutMap.flag[ic] = DC_NORMAL;
    }

    for (int j = 0; j < nface; j++)
    {
      int ic2 = c2c[nface*ic+j];
      if (ic2 >= 0 and cutMap.flag[ic2] == DC_UNASSIGNED)
      {
        cutMap.flag[ic2] = cutMap.flag[ic];
        icstack[nstack++] = ic2;
      }
    }
  }
#endif
}

void MeshBlock::unifyCutFlags(std::vector<CutMap> &cutMap)
{
  int nhole = 0; /// DEBUGGING
  for (int ic = 0; ic < ncells; ic++)
  {
    iblank_cell[ic] = NORMAL;
    for (int g = 0; g < cutMap.size(); g++)
    {
      if (cutMap[g].flag.size() < ncells) continue;

      if (cutMap[g].flag[ic] == DC_HOLE)
      {
        iblank_cell[ic] = HOLE;
        nhole++;
        break;
      }
    }
  }
  printf("%d: Final # of hole cells = %d\n",myid,nhole);
}

int get_cell_type(int* nc, int ntypes, int ic_in)
{
  // Return the cell type index and cell index within that type
  int count = 0;
  for (int n = 0; n < ntypes; n++)
  {
    if (count + nc[n] > ic_in)
      return n;

    count += nc[n];
  }
}

int get_cell_index(int* nc, int ntypes, int ic_in, int &ic_out)
{
  // Return the cell type index and cell index within that type
  int count = 0;
  for (int n = 0; n < ntypes; n++)
  {
    if (count + nc[n] > ic_in)
    {
      ic_out = ic_in - count;
      return n;
    }
    count += nc[n];
  }
}


void MeshBlock::getCellIblanks(const MPI_Comm meshComm)
{
  if (!iartbnd)
  {
    if (!iblank_cell) iblank_cell = (int *)malloc(sizeof(int)*ncells);
  }

  int icell = 0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      int flag = 1;
      iblank_cell[icell] = NORMAL;
      int ncount = 0;
      for (int m = 0; m < nvert; m++)
      {
        int inode = vconn[n][nvert*i+m]-BASE;
        if (iblank[inode] == HOLE)
        {
          iblank_cell[icell] = HOLE;
          flag = 0;
          break;
        }

        if (iblank[inode] == FRINGE) ncount++;
      }

      if (flag && ncount == nvert)
        iblank_cell[icell] = HOLE; // FRINGE

      icell++;
    }
  }

//  if (myid == 0) /// DEBUGGING
//  {
//    for (int i = 0; i < ncells; i++)
//      iblank_cell[i] = NORMAL;
//  }

//  /// DEBUGGING / TEMPORARY FIX
//  std::vector<int> ncf(1,6);
//  std::vector<std::vector<int>> c2c(1);
//  c2c[0].assign(ncells*6, -1);

//  for (int ic = 0; ic < ncells; ic++)
//  {
//    for (int j = 0; j < 6; j++)
//    {
//      int ff = c2f[6*ic+j];
//      int ic1 = f2c[2*ff];
//      int ic2 = f2c[2*ff+1];

//      if (ic1 == ic)
//        c2c[0][6*ic+j] = ic2;
//      else
//        c2c[0][6*ic+j] = ic1;
//    }
//  }

//  // New method to try: Expand hole regions by one layer where possible
//  // For all cells neighboring a hole cell:
//  //   If all surrounding cells are hole or fringe, set to hole cell
//  std::vector<int> iblank_tmp(ncells,0);

//  for (int ic = 0; ic < ncells; ic++)
//  {
//    int n = get_cell_type(nc,ntypes,ic);
//    if (iblank_cell[ic] == HOLE)
//    {
//      // Set all surrounding cells to 'HOLE' if no neighbors are 'NORMAL'
//      for (int j = 0; j < ncf[n]; j++)
//      {
//        int ic1 = c2c[n][ic*ncf[n]+j];
//        if (ic1 < 0 || iblank_cell[ic1] != FRINGE || iblank_tmp[ic1] == 1) // Already set to hole, or not possible to set as hole
//          continue;

//        bool set_hole = true;
//        int n1 = get_cell_type(nc,ntypes,ic1);
//        for (int k = 0; k < ncf[n]; k++)
//        {
//          int ic2 = c2c[n1][ic1*ncf[n1]+k];
//          if (ic2 < 0 || iblank_cell[ic2] == NORMAL)
//          {
//            set_hole = false;
//            break;
//          }
//        }

//        if (set_hole)
//          iblank_tmp[ic1] = 1;
//      }
//    }
//  }

//  for (int ic = 0; ic < ncells; ic++)
//  {
//    if (iblank_tmp[ic])
//      iblank_cell[ic] = HOLE;
//    else if (iblank_cell[ic] == FRINGE)
//      iblank_cell[ic] = NORMAL;
//  }

//  int rank;
//  MPI_Comm_rank(meshComm, &rank);
}

/**
void MeshBlock::getCellIblanks(const MPI_Comm meshComm)
{
  /// HACK FOR GOLF BALL CASE: Golf ball radius = .0625, outer radius = .14
  //double cut_rad2 = .1*.1;
  double cut_rad2 = .75*.75;
  //double cut_rad2 = 1*1;
  if (meshtag == 1) /// MUST BE BACKGROUND GRID
  {
    int icell = 0;
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for (int ic = 0; ic < nc[n]; ic++)
      {
        iblank_cell[icell] = NORMAL;

        double xc[3] = {0,0,0};
        for (int nd = 0; nd < 8; nd++)
        {
          int ind = 3*vconn[n][nvert*ic+nd];
//          for (int d = 0; d < 3; d++)
//            xc[d] += x[ind+d];

          double rad2 = 0;
          for (int d = 0; d < 3; d++)
            rad2 += x[ind+d]*x[ind+d];

          if (rad2 < cut_rad2)
          {
            iblank_cell[icell] = HOLE;
            break;
          }
        }
//        for (int d = 0; d < 3; d++)
//          xc[d] /= 8.;

//        double rad2 = 0;
//        for (int d = 0; d < 3; d++)
//          rad2 += xc[d]*xc[d];
////        rad2 = sqrt(rad2);

//        if (rad2 < cut_rad2)
//          iblank_cell[icell] = HOLE;

        icell++;
      }
    }
  }
  else
  {
    for (int ic = 0; ic < ncells; ic++)
      iblank_cell[ic] = NORMAL;
  }
}
*/

//void MeshBlock::getCellIblanks(const MPI_Comm meshComm)
//{
//  /// HACK FOR TAYLOR-GREEN VORTEX TEST CASE: INNER BOX (HALF-)LENGTH .25*PI
//  double safety_fac = .8;
//  double L = .25*PI * safety_fac;

//  if (meshtag == 1) /// MUST BE BACKGROUND GRID
//  {
//    int icell = 0;
//    for (int n = 0; n < ntypes; n++)
//    {
//      int nvert = nv[n];
//      for (int ic = 0; ic < nc[n]; ic++)
//      {
//        iblank_cell[icell] = NORMAL;

//        bool inside = true;
//        for (int nd = 0; nd < 8 && inside; nd++)
//        {
//          int ind = 3*vconn[n][nvert*ic+nd];
//          for (int d = 0; d < 3; d++)
//          {
//            if (std::abs(x[ind+d] - offset[d]) > L)
//            {
//              inside = false;
//              break;
//            }
//          }
//        }

//        if (inside)
//          iblank_cell[icell] = HOLE;

//        icell++;
//      }
//    }
//  }
//  else
//  {
//    for (int ic = 0; ic < ncells; ic++)
//      iblank_cell[ic] = NORMAL;
//  }
//}

void MeshBlock::calcFaceIblanks(const MPI_Comm &meshComm)
{
  nreceptorFaces = 0;

  std::set<int> artBndFaces;

  for (int ff = 0; ff < nfaces; ff++)
  {
    iblank_face[ff] = NORMAL;

    int ic1 = f2c[2*ff];   // By convention, always >= 0
    int ic2 = f2c[2*ff+1]; // By convention, == -1 if a boundary/MPI face

    if (ic2 < 0)
    {
      // Boundary or MPI face
      iblank_face[ff] = abs(iblank_cell[ic1]);
    }
    else
    {
      // Internal face
      int isum = abs(iblank_cell[ic1]) + abs(iblank_cell[ic2]);
      if (isum == HOLE+NORMAL)  // Artificial Boundary
      {
        iblank_face[ff] = FRINGE;
        artBndFaces.insert(ff);
      }
      else if (isum == HOLE+HOLE)  // Blanked face
      {
        iblank_face[ff] = HOLE;
      }
    }
  }

  // Now, ensure consistency of MPI face blanking across processes

  int meshRank, nProcMesh;
  MPI_Comm_rank(meshComm, &meshRank);
  MPI_Comm_size(meshComm, &nProcMesh);

  // Get the number of mpiFaces on each processor (for later communication)
  std::vector<int> nMpiFaces_proc(nProcMesh);
  MPI_Allgather(&nMpiFaces,1,MPI_INT,nMpiFaces_proc.data(),1,MPI_INT,meshComm);
  int maxNMpi = *std::max_element(nMpiFaces_proc.begin(), nMpiFaces_proc.end());

  std::vector<int> mpiIblank(nMpiFaces);
  std::vector<int> mpiIblank_proc(nProcMesh*maxNMpi);
  std::vector<int> mpiFid_proc(nProcMesh*maxNMpi);

  std::vector<int> recvCnts(nProcMesh);
  std::vector<int> recvDisp(nProcMesh);
  for (int i=0; i<nProcMesh; i++) {
    recvCnts[i] = nMpiFaces_proc[i];
    recvDisp[i] = i*maxNMpi;
  }

  for (int i = 0; i < nMpiFaces; i++)
    mpiIblank[i] = iblank_face[mpiFaces[i]];

  // Get iblank data for all mpi faces
  MPI_Allgatherv(mpiFaces, nMpiFaces, MPI_INT, mpiFid_proc.data(), recvCnts.data(), recvDisp.data(), MPI_INT, meshComm);
  MPI_Allgatherv(mpiIblank.data(), nMpiFaces, MPI_INT, mpiIblank_proc.data(), recvCnts.data(), recvDisp.data(), MPI_INT, meshComm);

  for (int F = 0; F < nMpiFaces; F++) {
    int ff = mpiFaces[F];
    int p = mpiProcR[F];
    for (int i = 0; i < nMpiFaces_proc[p]; i++) {
      if (mpiFid_proc[p*maxNMpi+i] == mpiFidR[F]) {
        int isum = mpiIblank[F] + mpiIblank_proc[p*maxNMpi+i];
        if (isum != NORMAL+NORMAL) {
          // Not a normal face; figure out if hole or fringe
          if (isum == HOLE+HOLE || iblank_cell[f2c[2*ff]] == HOLE)
            iblank_face[ff] = HOLE;
          else
          {
            iblank_face[ff] = FRINGE;
            artBndFaces.insert(ff);
          }
        }
      }
    }
  }

  // Lastly, set any explicit [solver-defined] overset faces (if not blanked)
  for (int i = 0; i < nOverFaces; i++)
  {
    int ff = overFaces[i];
    if (iblank_face[ff] != HOLE)
    {
      iblank_face[ff] = FRINGE;
      artBndFaces.insert(ff);
    }
  }

  // Setup final Artificial Boundary face list
  free(ftag);
  ftag = (int*)malloc(sizeof(int)*artBndFaces.size());

  for (auto &ff: artBndFaces) ftag[nreceptorFaces++] = ff;
}

void MeshBlock::setArtificialBoundaries(void)
{
  std::set<int> overFaces;

  for (int ff = 0; ff < nfaces; ff++)
  {
    if (iblank_face[ff] == FRINGE)
      overFaces.insert(ff);
  }
  nreceptorFaces = overFaces.size();

  // Setup final storage of A.B. face indices
  free(ftag);
  ftag = (int*)malloc(sizeof(int)*nreceptorFaces);

  int ind = 0;
  for (auto &ff: overFaces) ftag[ind++] = ff;
}

void MeshBlock::clearOrphans(int *itmp)
{
  if (iartbnd)
  {
    int fpt = 0;
    for (int i = 0; i < nreceptorFaces; i++)
    {
      bool reject = false;
      for (int j = 0; j < pointsPerFace[i]; j++)
      {
        if (itmp[fpt] == 0) reject = true;
        fpt++;
      }
      if (reject)
      {
        // Change both connected cells to hole cells, if not already
        /// TODO: this requires creation of new AB faces and additional flux
        /// point matching afterwards
        int ic1 = f2c[2*ftag[i]];
        int ic2 = f2c[2*ftag[i]+1];
        iblank_face[ftag[i]] = HOLE;
        if (ic1 >= 0)
          iblank_cell[ic1] = HOLE;
        if (ic2 >= 0)
          iblank_cell[ic2] = HOLE;
      }
    }
  }
  else if (ihigh)
    {
      int m=0;
      for (int i = 0; i < nreceptorCells;i++)
	{
	  bool reject = false;
	  for (int j = 0; j < pointsPerCell[i];j++)
	    {
	      if (itmp[m]==0) reject = true;
	      m++;
	    }
	  if (reject) 
	    {
	      iblank_cell[ctag[i]-1]=0; // changed to hole for off-body
	    }
	}
    }
  else
    {
      int m=0;
      for (int i = 0; i < nnodes;i++)
	{
	  if (picked[i]) 
	    {
	      bool reject = false;
	      if (itmp[m]==0) reject = true;
	      if (reject) iblank[i]=1; // changed to field for near-body
                                       // perhaps not the right thing to do
              m++;
	    }
	}
    }
}


void MeshBlock::getInternalNodes(void)
{
  nreceptorCells = 0;

  if (ihigh)
  {
    // Get a list of positions of all internal DOFs (solution points)
    // for all fringe cells
    free(ctag);
    ctag = (int *)malloc(sizeof(int)*ncells);

    // Gather a list of cell IDs for all receptor (fringe) cells
    for (int i = 0; i < ncells; i++)
      if (iblank_cell[i] == FRINGE)
        ctag[nreceptorCells++] = i+BASE;

    free(pointsPerCell);
    pointsPerCell = (int *)malloc(sizeof(int)*nreceptorCells);

    // Get total number of internal nodes (solution points) for our fringe cells
    maxPointsPerCell=0;
    ntotalPoints=0;

    for (int i = 0; i < nreceptorCells; i++)
    {
      get_nodes_per_cell(&(ctag[i]),&(pointsPerCell[i]));
      ntotalPoints += pointsPerCell[i];
      maxPointsPerCell = max(maxPointsPerCell,pointsPerCell[i]);
    }

    rxyz.resize(ntotalPoints*3);

    int m = 0;
    for (int i = 0; i < nreceptorCells; i++)
    {
      get_receptor_nodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
      m += (3*pointsPerCell[i]);
    }
  }
  else
  {
    ntotalPoints=0;
    free(picked);
    picked = (int *) malloc(sizeof(int)*nnodes);

    for (int i = 0; i < nnodes; i++) {
      picked[i] = 0;
      if (iblank[i] == FRINGE) {
        picked[i] = 1;
        ntotalPoints++;
      }
    }

    rxyz.resize(ntotalPoints*3);

    int m = 0;
    for (int i = 0; i < nnodes; i++)
    {
      if (picked[i])
      {
        int i3 = 3*i;
        for (int j = 0; j < 3; j++)
          rxyz[m++] = x[i3+j];
      }
    }
  }
}

void MeshBlock::getFringeNodes(bool unblanking)
{
  nreceptorFaces = 0;
  nreceptorCells = 0;
  ntotalPoints = 0;

  nFacePoints = 0;
  nCellPoints = 0;

  if (iartbnd && !unblanking)
  {
    // Gather all Artificial Boundary face flux points into fringe-point list
    free(ftag);
    ftag = (int*)malloc(sizeof(int)*nfaces);

    // Gather a list of cell IDs for all receptor (fringe) cells
    for (int i = 0; i < nfaces; i++)
      if (iblank_face[i] == FRINGE)
        ftag[nreceptorFaces++] = i+BASE;

    free(pointsPerFace);
    pointsPerFace = (int *)malloc(sizeof(int)*nreceptorFaces);

    // Get total number of face nodes (flux points) for our AB faces
    maxPointsPerFace = 0;

    for (int i = 0; i < nreceptorFaces; i++)
    {
      get_nodes_per_face(&(ftag[i]),&(pointsPerFace[i]));
      nFacePoints += pointsPerFace[i];
      maxPointsPerFace = max(maxPointsPerFace,pointsPerFace[i]);
    }

    rxyz.resize(nFacePoints*3);
    ntotalPoints = nFacePoints;

#ifdef _GPU
    get_face_nodes_gpu(ftag,nreceptorFaces,pointsPerFace,rxyz.data());
#else
    // Find the position of each flux point using callback function
    int m = 0;
    for (int i = 0; i < nreceptorFaces; i++)
    {
      get_face_nodes(&(ftag[i]),&(pointsPerFace[i]),&(rxyz[m]));
      m += (3*pointsPerFace[i]);
    }
#endif
  }

  if (ihigh)
  {
    /* Add in interior nodes from fringe elements (i.e. unblank cells, or
     * non-AB fringe cells) */
    if (iartbnd)
    {
      free(ctag);
      ctag = (int*)malloc(sizeof(int)*unblanks.size());

      // Gather a list of cell IDs for all receptor (fringe) cells
      for (auto &ic : unblanks)
        ctag[nreceptorCells++] = ic;
    }
    else
    {
      free(ctag);
      ctag = (int *)malloc(sizeof(int)*ncells);

      // Gather a list of cell IDs for all receptor (fringe) cells
      for (int i = 0; i < ncells; i++)
        if (iblank_cell[i] == FRINGE)
          ctag[nreceptorCells++] = i+BASE;
    }

    free(pointsPerCell);
    pointsPerCell = (int *)malloc(sizeof(int)*nreceptorCells);

    // Get total number of internal nodes (solution points) for our fringe cells
    maxPointsPerCell = 0;

    for (int i = 0; i < nreceptorCells; i++)
    {
      get_nodes_per_cell(&(ctag[i]),&(pointsPerCell[i]));
      nCellPoints += pointsPerCell[i];
      maxPointsPerCell = max(maxPointsPerCell,pointsPerCell[i]);
    }

    ntotalPoints = nCellPoints + nFacePoints;
    rxyz.resize(ntotalPoints*3);

#ifdef _GPU
    get_cell_nodes_gpu(ctag,nreceptorCells,pointsPerCell,&(rxyz[3*nFacePoints]));
#else
    int m = nFacePoints*3;
    for (int i = 0; i < nreceptorCells; i++)
    {
      get_receptor_nodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
      m += (3*pointsPerCell[i]);
    }
#endif
  }
  else
  {
    // Gather all fringe nodes into fringe-point list

    free(picked);
    picked = (int *) malloc(sizeof(int)*nnodes);

    for (int i = 0; i < nnodes; i++) {
      picked[i] = 0;
      if (iblank[i] == FRINGE) {
        picked[i] = 1;
        ntotalPoints++;
      }
    }

    rxyz.resize(ntotalPoints*3);

    int m = 0;
    for (int i = 0; i < nnodes; i++)
      if (picked[i])
        for (int j = 0; j < 3; j++)
          rxyz[m++] = x[3*i+j];
  }
}

void MeshBlock::getExtraQueryPoints(OBB *obc, int &nints, int *&intData,
                                    int &nreals, double *&realData)
{
  nints = 0;
  nreals = 0;

  // When high-order solvers being used, receptor points stored in 'rxyz'
  // [Recall: points loaded using the AB face list 'ftag' + callback funcs,
  // or using the list of fringe cells 'ctag' + callback funcs]
  std::vector<int> inode(ntotalPoints, -1);

  // Determine which of our receptor points overlap with given grid's OBB
  for (int i = 0; i < ntotalPoints; i++)
  {
    double xd[3] = {0,0,0};
    int i3 = 3*i;
    // Dot product of distance to OBB center with OBB axes
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        xd[j] += (rxyz[i3+k]-obc->xc[k])*obc->vec[j][k];

    if (fabs(xd[0]) <= obc->dxc[0] &&
        fabs(xd[1]) <= obc->dxc[1] &&
        fabs(xd[2]) <= obc->dxc[2])
    {
      inode[nints] = i;
      nints++;
      nreals += 3;
    }
  }

  intData = (int *)malloc(sizeof(int)*nints);
  realData = (double *)malloc(sizeof(double)*nreals);

  int m = 0;
  for (int i = 0; i < nints; i++)
  {
    int i3 = 3*inode[i];
    intData[i] = inode[i];
    realData[m++] = rxyz[i3];
    realData[m++] = rxyz[i3+1];
    realData[m++] = rxyz[i3+2];
  }
}

void MeshBlock::processPointDonorsGPU(void)
{
#ifdef _GPU
  interp2ListSize = 0;
  ninterp2 = 0;
  for (int i = 0; i < nsearch; i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]] == NORMAL)
      interp2ListSize++;

  interpList2.resize(interp2ListSize);

  if (interp2ListSize == 0) return;

  std::vector<double> rst2(rst.data(), rst.data()+rst.size());
  std::vector<int> donorId2(donorId.data(), donorId.data()+donorId.size());

  rst.resize(interp2ListSize*nDims);
  donorId.resize(interp2ListSize);

  int nWeightsTotal = 0;
  for (int i = 0; i < nsearch; i++)
  {
    if (donorId2[i] < 0 || iblank_cell[donorId2[i]] != NORMAL) continue;

    interpList2[ninterp2].nweights = get_n_weights(donorId2[i]);
    interpList2[ninterp2].receptorInfo[0] = isearch[2*i];   // procID
    interpList2[ninterp2].receptorInfo[1] = isearch[2*i+1]; // pointID
    interpList2[ninterp2].donorID = donorId2[i];
    nWeightsTotal += interpList2[ninterp2].nweights;

    for (int d = 0; d < 3; d++)
      rst[3*ninterp2+d] = rst2[3*i+d];

    donorId[ninterp2] = donorId2[i];

    ninterp2++;
  }

  // Get the interpolation weights for each of the interpolation points
  donors_d.assign(donorId.data(), donorId.size(), &stream_handle);
  mb_d.rst.assign(rst.data(), rst.size(), &stream_handle);

  weights_d.resize(nWeightsTotal);
  donor_frac_gpu(donors_d.data(), ninterp2, mb_d.rst.data(), weights_d.data());
#endif
}

void MeshBlock::processPointDonors(void)
{
  int buffsize = NFRAC;
  double* frac = (double *)malloc(sizeof(double)*buffsize);

  ninterp2 = 0;
  for (int i = 0; i < nsearch; i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]] == NORMAL)
      ninterp2++;

  interp2ListSize = ninterp2;
  interpList2.resize(interp2ListSize);

  int m = 0;
  for (int i = 0; i < nsearch; i++)
  {
    if (donorId[i] > -1 && iblank_cell[donorId[i]] == NORMAL)
    {
      if (ihigh)
      {
        int icell = donorId[i]+BASE;
        interpList2[m].inode.resize(1);
        interpList2[m].nweights = 0;
        interpList2[m].donorID = icell;

        // Use High-Order callback function to get interpolation weights
        donor_frac(&(icell), &(xsearch[3*i]), &(interpList2[m].nweights),
                   interpList2[m].inode.data(), frac, &(rst[3*i]), &buffsize);

        interpList2[m].weights.resize(interpList2[m].nweights);
        for(int j = 0; j < interpList2[m].nweights; j++)
          interpList2[m].weights[j] = frac[j];

        interpList2[m].receptorInfo[0] = isearch[2*i];   // procID
        interpList2[m].receptorInfo[1] = isearch[2*i+1]; // pointID

        m++;
      }
      else
      {
        int icell = donorId[i];
        int isum = 0;
        int n = -1;
        for (n = 0; n < ntypes; n++)
        {
          isum += nc[n];
          if (icell < isum)
          {
            icell = icell - (isum-nc[n]);
            break;
          }
        }

        int nvert = nv[n];
        interpList2[m].inode.resize(nvert);
        interpList2[m].nweights = nvert;
        interpList2[m].weights.resize(interpList2[m].nweights);

        double xv[8][3];
        for (int ivert = 0; ivert < nvert; ivert++)
        {
          interpList2[m].inode[ivert] = vconn[n][nvert*icell+ivert]-BASE;
          int i3 = 3*interpList2[m].inode[ivert];
          for (int j = 0; j < 3; j++)
            xv[ivert][j] = x[i3+j];
        }

        double xp[3];
        double frac2[8];

        xp[0] = xsearch[3*i];
        xp[1] = xsearch[3*i+1];
        xp[2] = xsearch[3*i+2];
        computeNodalWeights(xv,xp,frac2,nvert);

        for (int j = 0; j < nvert; j++)
          interpList2[m].weights[j] = frac2[j];         // interpolation weights

        interpList2[m].receptorInfo[0] = isearch[2*i];   // procID
        interpList2[m].receptorInfo[1] = isearch[2*i+1]; // pointID
        m++;
      }
    }
  }
  free(frac);
}

#ifdef _GPU
//void MeshBlock::setupBuffersGPU(int nsend, std::vector<int> &intData, std::vector<PACKET> &sndPack)
void MeshBlock::setupBuffersGPU(int nsend, std::vector<int> &intData, std::vector<VPACKET> &sndPack)
{
  if (ninterp2 == 0)
  {
    nSpts = 0;
    buf_disp.assign(nsend, 0);
    buf_inds.resize(0);
    intData.resize(0);

    sndPack.resize(nsend);
    for (int p = 0; p < nsend; p++)
    {
      sndPack[p].nints = 0;
      sndPack[p].nreals = 0;
      sndPack[p].intData.resize(0); ///
    }

    return;
  }

  nSpts = interpList2[0].nweights;

  std::vector<int> n_ints(nsend);
  intData.resize(ninterp2*2);
  for (int i = 0; i < ninterp2; i++)
  {
    int p = interpList2[i].receptorInfo[0];
    intData[2*i] = p;
    intData[2*i+1] = n_ints[p];
    n_ints[p]++;
  }

  sndPack.resize(nsend);
  for (int p = 0; p < nsend; p++)
  {
    sndPack[p].nints = n_ints[p];
//    if (sndPack[p].intData)
//      delete [] sndPack[p].intData;
//    sndPack[p].intData = new int[n_ints[p]];
    sndPack[p].intData.resize(n_ints[p]); ///
  }

  buf_disp.assign(nsend, 0);
  for (int p = 1; p < nsend; p++)
    buf_disp[p] = buf_disp[p-1] + n_ints[p-1];

  // Store final position within MPI buffer to place data
  buf_inds.resize(ninterp2);
  for (int i = 0; i < ninterp2; i++)
  {
    int p = intData[2*i];
    int ind = intData[2*i+1];
    buf_inds[i] = buf_disp[p] + ind;
    sndPack[p].intData[ind] = interpList2[i].receptorInfo[1];
  }
  buf_inds_d.assign(buf_inds.data(), buf_inds.size(), &stream_handle);
}

void MeshBlock::set_stream_handle(cudaStream_t stream, cudaEvent_t event)
{
  stream_handle = stream;
  event_handle = event;
}
#endif

void MeshBlock::getInterpolatedSolutionAtPoints(int *nints, int *nreals,
                int **intData, double **realData, double *q, int nvar,
                int interpType)
{
  (*nints) = ninterp2;
  (*nreals) = ninterp2*nvar;
  if ((*nints) == 0) return;

  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));

  std::vector<double> qq(nvar,0);
  int icount = 0;
  int dcount = 0;

  if (iartbnd)
  {
    for (int i = 0; i < ninterp2; i++)
    {
      qq.assign(nvar,0);
      int ic = interpList2[i].donorID;
      for (int spt = 0; spt < interpList2[i].nweights; spt++)
      {
        double weight = interpList2[i].weights[spt];
        for (int k = 0; k < nvar; k++) {
          double val = get_q_spt(ic,spt,k);
          qq[k] += val*weight;
        }
      }
      (*intData)[icount++] = interpList2[i].receptorInfo[0];
      (*intData)[icount++] = interpList2[i].receptorInfo[1];
      for (int k = 0; k < nvar; k++)
        (*realData)[dcount++] = qq[k];
    }
  }
  else if (ihigh)
  {
    if (interpType==ROW)
    {
      for (int i = 0; i < ninterp2; i++)
      {
        qq.assign(nvar,0);
        int inode = interpList2[i].inode[0]-BASE;
        for (int m = 0; m < interpList2[i].nweights; m++)
        {
          double weight = interpList2[i].weights[m];
          for (int k = 0; k < nvar; k++)
            qq[k] += q[inode+m*nvar+k]*weight;
        }
        (*intData)[icount++] = interpList2[i].receptorInfo[0];
        (*intData)[icount++] = interpList2[i].receptorInfo[1];
        for (int k = 0; k < nvar; k++)
          (*realData)[dcount++] = qq[k];
      }
    }
  }
  else
  {
    if (interpType==ROW)
    {
      for (int i = 0; i < ninterp2; i++)
      {
        for (int k = 0; k < nvar;k++) qq[k]=0;
        for (int m = 0; m < interpList2[i].nweights;m++)
        {
          int inode = interpList2[i].inode[m];
          double weight = interpList2[i].weights[m];
          if (weight < 0 || weight > 1.0) {
            traced(weight);
            printf("warning: weights are not convex\n");
          }
          for (int k = 0; k < nvar;k++)
            qq[k] += q[inode*nvar+k]*weight;
        }
        (*intData)[icount++] = interpList2[i].receptorInfo[0];
        (*intData)[icount++] = interpList2[i].receptorInfo[1];
        for (int k = 0; k < nvar;k++)
          (*realData)[dcount++] = qq[k];
      }
    }
    else if (interpType==COLUMN)
    {
      for (int i = 0; i < ninterp2; i++)
      {
        for (int k = 0; k < nvar;k++) qq[k]=0;
        for (int m = 0; m < interpList2[i].nweights; m++)
        {
          int inode = interpList2[i].inode[m];
          double weight = interpList2[i].weights[m];
          for (int k = 0; k < nvar; k++)
            qq[k] += q[k*nnodes+inode]*weight;
        }
        (*intData)[icount++] = interpList2[i].receptorInfo[0];
        (*intData)[icount++] = interpList2[i].receptorInfo[1];
        for (int k = 0; k < nvar; k++)
          (*realData)[dcount++] = qq[k];
      }
    }
  }
  //
  // no column-wise storage for high-order data
  //
}

void MeshBlock::getInterpolatedGradientAtPoints(int &nints, int &nreals,
                int *&intData, double *&realData, double *q, int nvar)
{
  nints = ninterp2;
  nreals = ninterp2*3*nvar;
  if (nints == 0) return;

  intData = (int *)malloc(sizeof(int)*2*nints);
  realData = (double *)malloc(sizeof(double)*nreals);

  std::vector<double> qq(3*nvar,0);
  int icount = 0;
  int dcount = 0;

  MPI_Pcontrol(1,"get_interp_grad");
  for (int i = 0; i < ninterp2; i++)
  {
    qq.assign(3*nvar,0);
    int ic = interpList2[i].donorID;
    for (int spt = 0; spt < interpList2[i].nweights; spt++)
    {
      double weight = interpList2[i].weights[spt];
      for (int dim = 0; dim < 3; dim++)
      {
        for (int k = 0; k < nvar; k++)
        {
          double val = get_grad_spt(ic,spt,dim,k);
          qq[dim*nvar+k] += val*weight;
        }
      }
    }
    intData[icount++] = interpList2[i].receptorInfo[0];
    intData[icount++] = interpList2[i].receptorInfo[1];
    for (int dim = 0; dim < 3; dim++)
      for (int k = 0; k < nvar; k++)
        realData[dcount++] = qq[dim*nvar+k];
  }
  MPI_Pcontrol(-1,"get_interp_grad");
}

void MeshBlock::getInterpolatedSolutionArtBnd(int &nints, int &nreals,
                std::vector<int> &intData, std::vector<double> &realData, int nvar)
{  
  nints = ninterp2;
  nreals = ninterp2*nvar;
  if (nints == 0) return;

  intData.resize(2*nints);
  realData.resize(nreals);

  int estride, sstride, vstride;
  double *q_spts = get_q_spts(estride, sstride, vstride);

  PUSH_NVTX_RANGE("getInterpU", 4);
  MPI_Pcontrol(1,"get_interpolated_U");
#pragma omp parallel for
  for (int i = 0; i < ninterp2; i++)
  {
    for (int k = 0; k < nvar; k++)
      realData[nvar*i+k] = 0.;

    double* eptr = q_spts + interpList2[i].donorID*estride;
    for (int spt = 0; spt < interpList2[i].nweights; spt++)
    {
      double weight = interpList2[i].weights[spt];
      for (int k = 0; k < nvar; k++)
      {
        realData[nvar*i+k] += weight * eptr[spt*sstride + k*vstride];
      }
    }
    intData[2*i] = interpList2[i].receptorInfo[0];
    intData[2*i+1] = interpList2[i].receptorInfo[1];
  }
  MPI_Pcontrol(-1,"get_interpolated_U");
  POP_NVTX_RANGE;
}

void MeshBlock::getInterpolatedGradientArtBnd(int &nints, int &nreals,
                std::vector<int> &intData, std::vector<double> &realData, int nvar)
{
  nints = ninterp2;
  nreals = ninterp2*nvar*nDims;
  if (nints == 0) return;

  intData.resize(2*nints);
  realData.resize(nreals);

  int estride, sstride, vstride, dstride;
  double *dq_spts = get_dq_spts(estride, sstride, vstride, dstride);

  PUSH_NVTX_RANGE("getInterpGrad", 5);
  MPI_Pcontrol(1, "get_interpolated_grad");
#pragma omp parallel for
  for (int i = 0; i < ninterp2; i++)
  {
    for (int dim = 0; dim < nDims; dim++)
      for (int k = 0; k < nvar; k++)
        realData[nvar*(nDims*i+dim)+k] = 0;

    double *eptr = dq_spts + interpList2[i].donorID*estride;
    for (int spt = 0; spt < interpList2[i].nweights; spt++)
    {
      double weight = interpList2[i].weights[spt];
      for (int dim = 0; dim < nDims; dim++)
        for (int k = 0; k < nvar; k++)
          realData[nvar*(nDims*i+dim)+k] += weight * eptr[spt*sstride + k*vstride + dim*dstride];
    }
    intData[2*i] = interpList2[i].receptorInfo[0];
    intData[2*i+1] = interpList2[i].receptorInfo[1];
  }
  MPI_Pcontrol(-1, "get_interpolated_grad");
  POP_NVTX_RANGE;
}

#ifdef _GPU
void MeshBlock::interpSolution_gpu(double *q_out_d, int nvar)
{
  if (ninterp2 == 0) return;

  int estride, sstride, vstride;
  double *q_d = get_q_spts_d(estride, sstride, vstride);

  // Perform the interpolation
  interp_u_wrapper(q_d, q_out_d, donors_d.data(), weights_d.data(),
      buf_inds_d.data(), ninterp2, nSpts, nvar, estride, sstride, vstride,
      stream_handle);
}

void MeshBlock::interpGradient_gpu(double *dq_out_d, int nvar)
{
  if (ninterp2 == 0) return;

  int estride, sstride, vstride, dstride;
  double *dq_d = get_dq_spts_d(estride, sstride, vstride, dstride);

  // Perform the interpolation
  interp_du_wrapper(dq_d, dq_out_d, donors_d.data(), weights_d.data(),
      buf_inds_d.data(), ninterp2, nSpts, nvar, nDims, estride, sstride,
      vstride, dstride, stream_handle);
}
#endif

void MeshBlock::updatePointData(double *q,double *qtmp,int nvar,int interptype)  
{
  int i,j,k,n,m;
  double *qout;
  int index_out;
  int npts;
  //
  if (ihigh) 
    {
      npts=NFRAC;
      qout=(double *)malloc(sizeof(double)*nvar*npts);	
      //
      m=0;
      for (int i = 0; i < nreceptorCells;i++)
	{
    if (iblank_cell[ctag[i]-1]==-1)
	    {
	      convert_to_modal(&(ctag[i]),&(pointsPerCell[i]),&(qtmp[m]),&npts,
			       &index_out,qout);
	      index_out-=BASE;
	      k=0;
	      for (int j = 0; j < npts;j++)
		for(n=0;n<nvar;n++)
		  {
		    q[index_out+j*nvar+n]=qout[k];
		    k++;
		  }
	    }
	  m+=(pointsPerCell[i]*nvar);
	}
      free(qout);
    }
  else
    {
      m=0;
      for (int i = 0; i < nnodes;i++)
	{
	  if (picked[i]) 
	    {
	      if (iblank[i]==-1) updateSolnData(i,&(qtmp[m]),q,nvar,interptype);
	      m+=nvar;
	    }
	}
    }
}

void MeshBlock::updateFringePointData(double *qtmp, int nvar)
{
  if (!ihigh) FatalError("updateFringePointData not applicable to non-high order solvers");

#ifdef _GPU
  if (nreceptorFaces > 0)
    face_data_to_device(ftag, nreceptorFaces, 0, qtmp);

  if (nreceptorCells > 0)
    cell_data_to_device(ctag, nreceptorCells, 0, qtmp+nvar*nFacePoints);
#else
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//  int fpt_start = 0;
#pragma omp parallel for
  for(int i = 0; i < nreceptorFaces; i++)
  {
    if (iblank_face[ftag[i]-BASE] == FRINGE)
    {
      int fpt_start = pointsPerFace[0]*nvar*i;
      for (int j = 0; j < pointsPerFace[i]; j++)
        for (int n = 0; n < nvar; n++)
          get_q_fpt(ftag[i], j, n) = qtmp[fpt_start+j*nvar+n];
    }
//    fpt_start += (pointsPerFace[i]*nvar);
  }
#endif
}

void MeshBlock::updateFringePointGradient(double *dqtmp, int nvar)
{
  if (!ihigh) FatalError("updateFringePointData not applicable to non-high order solvers");

  PUSH_NVTX_RANGE("tg_update_fringeGrad", 2);
#ifdef _GPU
  if (nreceptorFaces > 0)
    face_data_to_device(ftag, nreceptorFaces, 1, dqtmp);
#else
//  int fpt_start = 0;
#pragma omp parallel for
  for (int i = 0; i < nreceptorFaces; i++)
  {
    if (iblank_face[ftag[i]-BASE] == FRINGE)
    {
      int fpt_start = pointsPerFace[0]*nvar*nDims*i;
      for (int j = 0; j < pointsPerFace[i]; j++)
        for (int dim = 0; dim < 3; dim++)
          for (int n = 0; n < nvar; n++)
            get_grad_fpt(ftag[i], j, dim, n) = dqtmp[fpt_start+nvar*(j*3+dim)+n];
    }
//    fpt_start += (pointsPerFace[i]*nvar*3);
  }
#endif
  POP_NVTX_RANGE;
}

void MeshBlock::getDonorDataGPU(int dataFlag)
{
  MPI_Pcontrol(1, "getDonorDataGPU");

  // Get a sorted list of all donor cells on this rank
  std::vector<int> donorIDs(ninterp2);
  for (int i = 0; i < ninterp2; i++)
    donorIDs[i] = interpList2[i].donorID;

  std::sort(donorIDs.begin(), donorIDs.end());
  donorIDs.erase(std::unique(donorIDs.begin(),donorIDs.end()), donorIDs.end());

  data_from_device(donorIDs.data(), donorIDs.size(), dataFlag);

  MPI_Pcontrol(-1, "getDonorDataGPU");
}

void MeshBlock::sendFringeDataGPU(int gradFlag)
{
  MPI_Pcontrol(1, "sendFringeDataGPU");
  face_data_to_device(ftag, nreceptorFaces, gradFlag, NULL);
  MPI_Pcontrol(-1, "sendFringeDataGPU");
}
