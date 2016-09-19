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

#include <algorithm>
#include <set>

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

/*
void MeshBlock::directCut(int nGroups, int* groupIDs, int* cutType, int* nGf,
    int** cutFaces, int gridID, int nGrids, MPI_Comm &scomm)
{
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
   *
  int nDims = 3; /// TODO: add to MB.h ...

  /// TODO: pre-process groups on this rank vs. groups on ALL ranks/grids
  /// TODO: Move up to Tioga class
  int maxID = 0;
  for (int g = 0; g < nGroups; g++)
    maxID = max(maxID, groupIDs[g]);

  int nGroups_glob = maxID;
  MPI_Allreduce(&maxId, &nGroups_glob, 1, MPI_INT, MPI_MAX, scomm);

  // Full list of group types, and list of group IDs for each grid
  std::vector<int> cutType_glob(nGroups_glob, -1);
  std::vector<int> gridGroups(nGrids*nGroups_glob);
  std::unordered_set<int> myGroups;

  for (int g = 0; g < nGroups; g++)
  {
    int G = groupIDs[g];
    cutType_glob[G] = cutType[g];
    gridGroups[gridID*nGroups_glob+G] = 1;
    myGroups.insert(G);
  }

  MPI_Allreduce(&cutType_glob[0], MPI_IN_PLACE, nGroups_glob, MPI_INT, MPI_MAX, scomm);
  MPI_Allreduce(&gridGroups[0], MPI_IN_PLACE, nGrids*nGroups_glob, MPI_INT, MPI_MAX, scomm);

  std::vector<double> cutBox(6*nGroups_glob);
  std::vector<std::vector<double>> faceBox(nGroups);


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
    for (int i = 0; i < nGf[g]; i++)
    {
      // Initialize this face's bounding box
      for (int d = 0; d < nDims; d++)
      {
        faceBox[g][6*i+d]   = BIGVALUE;
        faceBox[g][6*i+d+3] = -BIGVALUE;
      }

      // Loop over vertices of this face
      int ff = cutFaces[g][i];
      for (int n = 0; n < nfv[0]; n++)
      {
        int iv = f2v[ff*nfv[0] + n];
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

  /// TODO: re-locate up to Tioga class?
  // Get the global bounding box info across all the partitions for all meshes
  std::vector<double> cutBox_global(6*nGroups_glob);
  for (int G = 0; G < nGroups_glob; G++)
  {
    MPI_Allreduce(&cutBox[6*G],  &cutBox_global[6*G],  3,MPI_DOUBLE,MPI_MIN,scomm);
    MPI_Allreduce(&cutBox[6*G+3],&cutBox_global[6*G+3],3,MPI_DOUBLE,MPI_MAX,scomm);
  }

  /* Use ADT Search to find all cells [on this rank] which intersect with each
   * cutting group [which is not on this rank]
   * NOTE: Only 'background' grids (defined as grids with no cutting groups)
   * will be cut by 'field'-type groups *
  std::vector<std::unordered_set<int>> cellList(nGroups_glob);
  for (int G = 0; G < nGroups_glob; G++)
  {
    if (myGroups.count(G)) continue;
    if (nGroups > 0 && cutType_glob[G] == FIELD_CUT) continue;
    adt->searchADT_box(elementList,cellList[G],&cutBox_global[6*G]);
  }

  // use 'PC' object to send/recv faceBox data to correct ranks
  // Ranks for which cellList[G] != 0 must send data
  int nsend = 0;
  int nrecv = 0;
  std::vector<int> sendMap, sendInts, recvMap, recvGroups;
  for (int p = 0; p < numprocs; p++)
  {
    int grid = gridIds[p];
    for (int G = 0; G < nGroups_glob; G++)
    {
      if (gridGroups[grid*nGroups_glob+G]) // && cellList[G].size() > 0)
      {
        nrecv++;
        nsend++;
        recvMap.push_back(p);
        recvGroups.push_back(G);
        sendMap.push_back(p);
        sendInts.push_back(cellList[G].size());
      }
    }
  }
}
*/

void MeshBlock::getCellIblanks(void)
{
  int icell = 0;
  if (!iartbnd)
  {
    if (!iblank_cell) iblank_cell = (int *)malloc(sizeof(int)*ncells);
  }

  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      int flag = 1;
      iblank_cell[icell] = NORMAL;
      int ncount = 0;
      for (int m = 0; m < nvert && flag; m++)
      {
        int inode = vconn[n][nvert*i+m]-BASE;
        if (iblank[inode] == HOLE)
        {
          iblank_cell[icell] = HOLE;
          flag = 0;
        }
        ncount = ncount + (iblank[inode] == FRINGE);
      }

//      if (flag && ncount == nvert)
//        iblank_cell[icell] = FRINGE;

      icell++;
    }
  }
}

void MeshBlock::calcFaceIblanks(const MPI_Comm &meshComm)
{
  nreceptorFaces = 0;

  // First, correct iblank_cell to contain only normal or hole cells (no fringe)
  for (int ic = 0; ic < ncells; ic++)
  {
    if (iblank_cell[ic] == FRINGE)
      iblank_cell[ic] = HOLE;
  }

  std::set<int> artBndFaces;

  for (int ff = 0; ff < nfaces; ff++)
  {
    iblank_face[ff] = NORMAL;

    int ic1 = f2c[2*ff];   // By convention, always >= 0
    int ic2 = f2c[2*ff+1]; // By convention, == -1 if a boundary/MPI face

    if (ic2 < 0)
    {
      // Boundary or MPI face
      iblank_face[ff] = iblank_cell[ic1];
    }
    else
    {
      // Internal face
      int isum = iblank_cell[ic1] + iblank_cell[ic2];
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

  nreceptorFaces = artBndFaces.size();
//printf("Rank %d: Grid %d: # Overset Faces = %d\n",myid,meshtag,nreceptorFaces);
  // Setup final Artificial Boundary face list
  free(ftag);
  ftag = (int*)malloc(sizeof(int)*nreceptorFaces);

  int ind = 0;
  for (auto &ff: artBndFaces) ftag[ind++] = ff;
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

    free(rxyz);
    //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
    rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);

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

    free(rxyz);
    rxyz = (double *)malloc(sizeof(double)*ntotalPoints*3);

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

void MeshBlock::getBoundaryNodes(void)
{
  nreceptorFaces = 0;

  if (iartbnd)
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
    ntotalPoints = 0;

    for (int i = 0; i < nreceptorFaces; i++)
    {
      get_nodes_per_face(&(ftag[i]),&(pointsPerFace[i]));
      ntotalPoints += pointsPerFace[i];
      maxPointsPerFace = max(maxPointsPerFace,pointsPerFace[i]);
    }

    free(rxyz);
    rxyz = (double *)malloc(sizeof(double)*ntotalPoints*3);

    // Find the position of each flux point using callback function
    int m = 0;
    for (int i = 0; i < nreceptorFaces; i++)
    {
      get_face_nodes(&(ftag[i]),&(pointsPerFace[i]),&(rxyz[m]));
      m += (3*pointsPerFace[i]);
    }
  }
  else
  {
    // Gather all fringe nodes into fringe-point list
    ntotalPoints = 0;

    free(picked);
    picked = (int *) malloc(sizeof(int)*nnodes);

    for (int i = 0; i < nnodes; i++) {
      picked[i] = 0;
      if (iblank[i] == FRINGE) {
        picked[i] = 1;
        ntotalPoints++;
      }
    }

    free(rxyz);
    rxyz = (double *)malloc(sizeof(double)*ntotalPoints*3);

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

void MeshBlock::processPointDonors(void)
{
  int buffsize = NFRAC;
  double* frac = (double *)malloc(sizeof(double)*buffsize);
  interp2ListSize = ninterp2;
  ninterp2 = 0;

  for (int i = 0; i < nsearch; i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]] == NORMAL)
      ninterp2++;

  if (interpList2)
  {
    for (int i = 0; i < interp2ListSize; i++)
    {
      free(interpList2[i].inode);
      free(interpList2[i].weights);
    }
    free(interpList2);
  }

  interp2ListSize = ninterp2;
  interpList2 = (INTERPLIST *)malloc(sizeof(INTERPLIST)*interp2ListSize);

  for (int i = 0; i < interp2ListSize; i++)
  {
    interpList2[i].inode=NULL;
    interpList2[i].weights=NULL;
  }

  int m = 0;
  for (int i = 0; i < nsearch; i++)
  {
    if (donorId[i] > -1 && iblank_cell[donorId[i]] == NORMAL)
    {
      if (ihigh)
      {
        int icell = donorId[i]+BASE;
        interpList2[m].inode = (int *) malloc(sizeof(int));
        interpList2[m].nweights = 0;
        interpList2[m].donorID = icell;

        // Use High-Order callback function to get interpolation weights
        donor_frac(&(icell), &(xsearch[3*i]), &(interpList2[m].nweights),
                   interpList2[m].inode, frac, &(rst[3*i]), &buffsize);

        interpList2[m].weights = (double *)malloc(sizeof(double)*interpList2[m].nweights);
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
        interpList2[m].inode = (int *) malloc(sizeof(int)*nvert);
        interpList2[m].nweights = nvert;
        interpList2[m].weights = (double *)malloc(sizeof(double)*interpList2[m].nweights);

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
      for (int i = 0; i < ninterp2;i++)
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
}

void MeshBlock::getInterpolatedSolutionArtBnd(int &nints, int &nreals,
                std::vector<int> &intData, std::vector<double> &realData, int nvar)
{
  nints = ninterp2;
  nreals = ninterp2*nvar;
  if (nints == 0) return;

  intData.resize(2*nints);
  realData.resize(nreals);

  std::vector<double> qq(nvar,0);
  int icount = 0;
  int dcount = 0;

  for (int i = 0; i < ninterp2; i++)
  {
    qq.assign(nvar,0);
    int ic = interpList2[i].donorID;
    for (int spt = 0; spt < interpList2[i].nweights; spt++)
    {
      double weight = interpList2[i].weights[spt];
      for (int k = 0; k < nvar; k++)
        qq[k] += weight * get_q_spt(ic,spt,k);
    }
    intData[icount++] = interpList2[i].receptorInfo[0];
    intData[icount++] = interpList2[i].receptorInfo[1];
    for (int k = 0; k < nvar; k++)
      realData[dcount++] = qq[k];
  }
}
	
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

void MeshBlock::updateFluxPointData(double *qtmp, int nvar)
{
  if (!ihigh) FatalError("updateFluxPointData not applicable to non-high order solvers");

  int fpt_start = 0;
  for(int i = 0; i < nreceptorFaces; i++)
  {
    if (iblank_face[ftag[i]-BASE] == -1)
    {
      for (int j = 0; j < pointsPerFace[i]; j++)
        for (int n = 0; n < nvar; n++)
          get_q_fpt(ftag[i], j, n) = qtmp[fpt_start+j*nvar+n];
    }
    fpt_start += (pointsPerFace[i]*nvar);
  }
}

void MeshBlock::updateFluxPointGradient(double *dqtmp, int nvar)
{
  if (!ihigh) FatalError("updateFluxPointData not applicable to non-high order solvers");

  int fpt_start = 0;
  for(int i = 0; i < nreceptorFaces; i++)
  {
    if (iblank_face[ftag[i]-BASE] == -1)
    {
      for (int j = 0; j < pointsPerFace[i]; j++)
        for (int dim = 0; dim < 3; dim++)
          for (int n = 0; n < nvar; n++)
            get_grad_fpt(ftag[i], j, dim, n) = dqtmp[fpt_start+nvar*(j*3+dim)+n];
    }
    fpt_start += (pointsPerFace[i]*nvar*3);
  }
}

void MeshBlock::getDonorDataGPU(int dataFlag)
{
  // Get a sorted list of all donor cells on this rank
  std::set<int> donors;
  for (int i = 0; i < ninterp2; i++)
    donors.insert(interpList2[i].donorID);

  std::vector<int> donorIDs;
  donorIDs.reserve(donors.size());
  for (auto &ic:donors) donorIDs.push_back(ic);

  data_from_device(donorIDs.data(), donorIDs.size(), dataFlag);
}

void MeshBlock::sendFringeDataGPU(int gradFlag)
{
  data_to_device(ftag, nreceptorFaces, gradFlag);
}
