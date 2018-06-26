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
#include "utils.h"
#include "math_funcs.h"
#include "linklist.h"

#define verbose false

using namespace tg_funcs;

void MeshBlock::getDonorPacket(PACKET *sndPack, int nsend)
{
  int* icount = (int *)malloc(sizeof(int)*nsend);
  int* dcount = (int *)malloc(sizeof(int)*nsend);

  // count numbers to send first
  for (int i = 0; i < nsearch; i++)
  {
    if (donorId[i] > -1)
    {
      int k = isearch[2*i];
      sndPack[k].nints  += 3;
      sndPack[k].nreals ++; // Donor res
    }
  }

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].intData = (int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData = (double *)malloc(sizeof(double)*sndPack[k].nreals);
  }

  for (int i = 0; i < nsend; i++)
    icount[i] = dcount[i] = 0;

  for (int i = 0; i < nsearch; i++)
  {
    if (donorId[i] > -1)
    {
      int k = isearch[2*i];
      sndPack[k].intData[icount[k]++]  = meshtag;             // mesh tag
      sndPack[k].intData[icount[k]++]  = isearch[2*i+1];      // point id from receptor grid
      sndPack[k].intData[icount[k]++]  = i;                   // point id on the donor side
      sndPack[k].realData[dcount[k]++] = cellRes[donorId[i]]; // donor resolution
    }
  }

  free(icount);
  free(dcount);
}
void MeshBlock::initializeDonorList(void)
{
  if (donorList)
  {
    for (int i = 0; i < donorListLength; i++)
      deallocateLinkList(donorList[i]);

    free(donorList);
  }

  donorListLength = nnodes;
  donorList = (DONORLIST **)malloc(sizeof(DONORLIST *)*donorListLength);
  for (int i = 0; i < donorListLength; i++)
    donorList[i] = NULL;
}

void MeshBlock::insertAndSort(int pointid,int senderid,int meshtagdonor, int remoteid,
            double donorRes)
{
  DONORLIST *temp1;
  temp1=(DONORLIST *)malloc(sizeof(DONORLIST));
  temp1->donorData[0] = senderid;
  temp1->donorData[1] = meshtagdonor;
  temp1->donorData[2] = remoteid;
  temp1->donorRes     = donorRes;
  insertInList(&donorList[pointid],temp1);
}

void MeshBlock::processDonors(HOLEMAP *holemap, int nmesh, int **donorRecords,double **receptorResolution,
			      int *nrecords)
{
  // First mark all hole points
  std::vector<int> iflag(nmesh);

  for (int i = 0; i < nnodes; i++)
  {
    iblank[i] = NORMAL;
    if (verbose) tracei(i);
    if (donorList[i]==NULL)
    {
      if (verbose) printf("No donor found for %d\n", i);
      for (int j = 0; j < nmesh; j++)
        if (j != (meshtag-BASE) && holemap[j].existWall)
        {
          if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents))
          {
            iblank[i] = HOLE;
            break;
          }
        }
    }
    else
    {
      DONORLIST *temp = donorList[i];

      for (int j = 0; j < nmesh; j++) iflag[j] = 0;

      while (temp != NULL)
      {
        int meshtagdonor = temp->donorData[1]-BASE;
        iflag[meshtagdonor] = 1;
        if (verbose) {
          tracei(meshtagdonor);
          traced(temp->donorRes);
        }
        temp = temp->next;
      }
      for (int j = 0; j < nmesh; j++)
      {
        if (j != (meshtag-BASE) && holemap[j].existWall)
        {
          if (!iflag[j])
            if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam,holemap[j].extents))
            {
              iblank[i] = HOLE;
              break;
            }
        }
      }
    }
  }

  for (int i = 0; i < nwbc; i++)
  {
   if (iblank[wbcnode[i]-BASE] == HOLE) {
     printf("--------------------------------------------------------------------\n");
     printf("Alarm from process %d : wall node is being tagged as a hole %d %p\n",myid,wbcnode[i]-BASE,
        donorList[wbcnode[i]-BASE]);
     int ii=wbcnode[i]-BASE;
     printf("xloc=%e %e %e\n",x[3*ii],x[3*ii+1],x[3*ii+2]);
     printf("Computations will continue, but may suffer from accuracy problems\n");
     printf("Please recheck positions of your grids\n");
     printf("--------------------------------------------------------------------\n");
    }
  }

  // Mark mandatory fringes as neighbors (up to nfringe depth) of hole points
  int *mtag = (int *)malloc(sizeof(int)*nnodes);
  int *mtag1 = (int *)malloc(sizeof(int)*nnodes);
 
  for (int i = 0; i < nnodes; i++)
  {
    mtag[i] = mtag1[i] = HOLE;
    if (iblank[i] == HOLE)
      mtag[i] = mtag1[i] = NORMAL;
  }

  for (int iter = 0; iter < nfringe; iter++)
  {
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for (int i = 0; i < nc[n]; i++)
      {
        for (int m = 0; m < nvert; m++)
        {
          if (mtag[(vconn[n][nvert*i+m]-BASE)] == NORMAL)
          {
            for (int mm = 0; mm < nvert; mm++)
              if (m != mm && mtag[vconn[n][nvert*i+mm]-BASE] != NORMAL)
                mtag1[vconn[n][nvert*i+mm]-BASE] = NORMAL;
          }
        }
      }
    }

    for (int i = 0; i < nnodes; i++) mtag[i] = mtag1[i];
  }

  for (int i = 0; i < nnodes; i++)
    if (mtag1[i] && iblank[i])
      nodeRes[i]=BIGVALUE;

  free(mtag);
  free(mtag1);

  // Now find fringes
  *nrecords=0;
  for (int i = 0; i < nnodes; i++)
  {
    //if (myid==553 && i==29670) verbose=1;
    if (verbose) {
       tracei(i);
       tracei(iblank[i]);
    }
    if (donorList[i] != NULL && iblank[i] != HOLE)
    {
      DONORLIST *temp = donorList[i];
      if (verbose) traced(nodeRes[i]);
      while (temp != NULL)
      {
        if (verbose) traced(temp->donorRes);
        if (temp->donorRes < nodeRes[i])
        {
          iblank[i] = -temp->donorData[1] - 1; // meshtag range [0,nmesh-1]; subtracting 1 to avoid '-0' (want actual negative value)
          if (verbose) {
          tracei(iblank[i]);
          tracei(temp->donorData[0]);
          tracei(temp->donorData[1]);
          tracei(temp->donorData[2]);}
          (*nrecords)++;
          break;
        }
        temp = temp->next;
      }
    }
  }

  // set the records to send back to the donor
  // process
  (*donorRecords) = (int *)malloc(sizeof(int)*2*(*nrecords));
  (*receptorResolution) = (double *)malloc(sizeof(double)*(*nrecords));

  int m = 0;
  int k = 0;
  for (int i = 0; i < nnodes; i++)
  {
    if (iblank[i] < 0)
    {
      DONORLIST *temp = donorList[i];
      while(temp != NULL)
      {
        if (temp->donorRes < nodeRes[i])
        {
          break;
        }
      }
      //(*receptorResolution)[k++] = nodeRes[i];
      //if (temp->donorRes < 0) nodeRes[i] = BIGVALUE;
      (*receptorResolution)[k++] = (resolutionScale > 1.0) ? -nodeRes[i] : nodeRes[i];
      (*donorRecords)[m++] = temp->donorData[0];
      (*donorRecords)[m++] = temp->donorData[2];
      if (verbose) {
        tracei(iblank[i]);
        tracei(m);
        tracei((*donorRecords)[m-1]);
        tracei((*donorRecords)[m-2]);
        traced((*receptorResolution)[k-1]);
      }
    }
  }
}

void MeshBlock::initializeInterpList(int ninterp_input)
{
  ninterp = ninterp_input;
  interpListSize = ninterp_input;

  interpList.resize(interpListSize);

  deallocateLinkList2(cancelList);
  cancelList = NULL;
  ncancel = 0;

  free(interp2donor);
  interp2donor = (int *)malloc(sizeof(int)*nsearch);
  for(int i = 0; i < nsearch; i++)
    interp2donor[i]=-1;
    
}
		
void MeshBlock::findInterpData(int &recid,int irecord,double receptorRes2)
{
  double xv[8][3];
  std::vector<double> xv2;
  std::vector<double> frac;
  std::vector<int> inode;
  INTEGERLIST *clist;

  double receptorRes = fabs(receptorRes2);
  int procid = isearch[2*irecord];
  int pointid = isearch[2*irecord+1];
  int meshtagrecv = tagsearch[irecord];
  int i3 = 3*irecord;

  double xp[3];
  xp[0] = xsearch[i3];
  xp[1] = xsearch[i3+1];
  xp[2] = xsearch[i3+2];

  int isum = 0;
  int n = 0;
  int idonor = 0;
  for (n = 0; n < ntypes; n++)
  {
    isum += nc[n];
    if (donorId[irecord] < isum)
    {
      idonor = donorId[irecord]-(isum-nc[n]);
      break;
    }
  }

  // Get the physical node positions of the donor cell
  int acceptFlag = 1;
  int nvert = nv[n];

  if (nvert > 8)
    nvert = nNodesToFirstOrder(ncf[n], nv[n]);

  inode.resize(nv[n]);
  if (nv[n] > 8)
    xv2.resize(nv[n]*3);

  for (int m = 0; m < nv[n]; m++)
  {
    inode[m]=vconn[n][nv[n]*idonor+m]-BASE;
    i3 = 3*inode[m];
    //if (iblank[inode[m]] != NORMAL)
    if (iblank[inode[m]] <= 0 && receptorRes2 > 0.0) // || nodeRes[inode[m]] == BIGVALUE)
    {
      if (nodeRes[inode[m]]==BIGVALUE) acceptFlag=0;
      if (abs(iblank[inode[m]])==meshtagrecv+1) acceptFlag=0;
      if (iblank[inode[m]]==0) acceptFlag=0;
    }

    if (nvert > 8) printf("ERROR! Shouldn't have nvert > 8 now.  nvert = %d\n",nvert);
//      for (int j = 0; j < 3; j++)
//        xv2[m*3+j] = x[i3+j];
//    else
    if (m < nvert)
      for (int j = 0; j < 3; j++)
        xv[m][j] = x[i3+j];
  }

  if (verbose) tracei(acceptFlag);
  if (acceptFlag == 0 && receptorRes2 != BIGVALUE) return;

  if (receptorRes2 == BIGVALUE) /// && !iartbnd) // Node tagged as 'mandatory fringe' receptor
  {
    clist=cancelList;

    // go to the end of the list
    if (clist !=NULL) while(clist->next !=NULL) clist=clist->next;

    for (int m = 0; m < nv[n]; m++)
    {
      inode[m] = vconn[n][nv[n]*idonor+m]-BASE;
      if (verbose) tracei(inode[m]);
      if (verbose) traced(nodeRes[inode[m]]);
      if (verbose) {
        tracei(procid);
        tracei(pointid);
        traced(receptorRes);
        tracei(irecord);
        tracei(donorId[irecord]);
      }
      if (iblank[inode[m]] < 0 && nodeRes[inode[m]] != BIGVALUE)
      {
        iblank[inode[m]]=NORMAL;
        if (clist == NULL)
        {
          clist=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
          clist->inode=inode[m];
          clist->next=NULL;
          cancelList=clist;
        }
        else
        {
          clist->next=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
          clist->next->inode=inode[m];
          clist->next->next=NULL;
          clist=clist->next;
        }
        ncancel++;
      }
    }
  }

  /* Compute the (linear) shape-function interpolation weights of the receptor
   * node within the donor cell */
  frac.resize(nvert);

//  if (nvert <= 8)
//  {
    computeNodalWeights(xv,xp,frac.data(),nvert); /// TODO: support nvert > 8
//  }
//  else
//  {
//    double xref[3];
//    getRefLocNewton(xv2.data(), &xp[0], &xref[0], nvert, 3);  /// TODO: mixed grids; mayabe use callback...?
//    shape_hex(point(&xref[0]),frac,nvert);
//  }

  interp2donor[irecord] = recid;
  interpList[recid].cancel = 0;
  interpList[recid].nweights = nvert;
  interpList[recid].receptorInfo[0] = procid;
  interpList[recid].receptorInfo[1] = pointid;
  if (verbose) {
    tracei(interpList[recid].receptorInfo[0]);
    tracei(interpList[recid].receptorInfo[1]);
  }
  interpList[recid].inode.resize(nvert);
  interpList[recid].weights.resize(nvert);
  for (int m = 0; m < nvert; m++)
  {
    interpList[recid].inode[m]   = inode[m];
    interpList[recid].weights[m] = frac[m];
  }
  recid++;
}

void MeshBlock::set_ninterp(int ninterp_input)
{
  ninterp=ninterp_input;
}
  
void MeshBlock::getCancellationData(int &nrecords,int*& intData)
{
  nrecords = ncancel;
  if (ncancel > 0)
  {
    intData = (int *)malloc(sizeof(int)*nrecords*2);
    int i = 0;
    for (INTEGERLIST* clist = cancelList; clist != NULL; clist = clist->next)
    {
      int inode = clist->inode;
      intData[i++] = donorList[inode]->donorData[0];
      intData[i++] = donorList[inode]->donorData[2];
    }
  }
}

void MeshBlock::cancelDonor(int irecord)
{
  int iptr = interp2donor[irecord];
  if (iptr > -1) interpList[iptr].cancel = 1;
}

void MeshBlock::getInterpData(int& nrecords, int*& donorData)
{
  nrecords = 0;
  for (int i = 0; i < ninterp; i++)
    if (!interpList[i].cancel) nrecords++;

  donorData = (int *)malloc(sizeof(int)*2*nrecords);

  for (int i = 0, k = 0; i < ninterp; i++)
  {
    if (!interpList[i].cancel)
    {
       donorData[k++] = interpList[i].receptorInfo[0];
       donorData[k++] = interpList[i].receptorInfo[1];
    }
  }
}

void MeshBlock::clearIblanks(void)
{
  int i;
  for(i=0;i<nnodes;i++)
     if (iblank[i] < 0) iblank[i]=1;
}

void MeshBlock::getStats(int mstats[2])
{
  int i;
  mstats[0]=mstats[1]=0;
  for (i=0;i<nnodes;i++)
    {
      if (iblank[i]==0) mstats[0]++;
      if (iblank[i] < 0) mstats[1]++;
    }
}

void MeshBlock::setIblanks(int inode)
{
  iblank[inode] = FRINGE;
}
