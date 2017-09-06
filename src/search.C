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
#include "tioga.h" /// DEBUGGING (profiling - timer)
#include "MeshBlock.h"
extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  void writebbox(OBB *obb,int bid);
  void writePoints(double *x,int nsearch,int bid);
}

void MeshBlock::search(void)
{
  if (nsearch == 0) {
    donorCount=0;
    return;
  }

  if (rrot && adt)
  {
    //rebuildADT();
    /// TODO: worry about updating ADT based on changing OBB of search points
    adt->setTransform(Rmat, offset, nDims);
#ifdef _GPU
    adt_d.setTransform(Rmat, offset, nDims);
    mb_d.setTransform(Rmat, offset, nDims);
    mb_d.updateSearchPoints(nsearch,isearch.data(),xsearch.data());
#endif
  }
  else
  {
    double xd[3];
    double dxc[3];
    double xmin[3];
    double xmax[3];

    // form the bounding box of the
    // query points
    OBB obq;
    findOBB(xsearch.data(),obq.xc,obq.dxc,obq.vec,nsearch);

    // find all the cells that may have intersections with
    // the OBB
    std::vector<int> icell(ncells, -1);

    ncells_adt = 0;
    int ic = 0;
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for(int i = 0; i < nc[n]; i++)
      {
        // find each cell that has
        // overlap with the bounding box
        xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
        xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
        for (int m = 0; m < nvert; m++)
        {
          int i3 = 3*(vconn[n][nvert*i+m]-BASE);
          for (int j = 0; j < 3; j++)
          {
            xd[j] = 0;
            for (int k = 0; k < 3; k++)
              xd[j] += (x[i3+k]-obq.xc[k])*obq.vec[j][k];
            xmin[j] = min(xmin[j],xd[j]);
            xmax[j] = max(xmax[j],xd[j]);
          }
          for (int j = 0; j < 3; j++)
          {
            xd[j] = (xmax[j]+xmin[j])*0.5;
            dxc[j] = (xmax[j]-xmin[j])*0.5;
          }
        }
        /// HACK - commenting out this if statement makes this loop basically useless (ncells_adt will == ncells)
        if (fabs(xd[0]) <= (dxc[0]+obq.dxc[0]) &&
            fabs(xd[1]) <= (dxc[1]+obq.dxc[1]) &&
            fabs(xd[2]) <= (dxc[2]+obq.dxc[2]))
        {
          // create a list of all the cells that have bounding box intersection
          // with the QP bounding box
          icell[ncells_adt++] = ic;
        }
        ic++;
      }
    }

    // now find the axis aligned bounding box of each cell in the list to
    // build the ADT
    elementBbox.resize(ncells_adt*6);
    elementList.resize(ncells_adt);

    for (int adt_ic = 0; adt_ic < ncells_adt; adt_ic++)
    {
      int loc_ic = -1;
      int mb_ic = icell[adt_ic];
      int isum = 0;
      int n;
      for (n = 0; n < ntypes; n++)
      {
        isum += nc[n];
        if (mb_ic < isum)
        {
          loc_ic = mb_ic-(isum-nc[n]);
          break;
        }
      }

      int nvert = nv[n];
      xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
      xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
      for (int m = 0; m < nvert; m++)
      {
        int i3 = 3*(vconn[n][nvert*loc_ic+m]-BASE);
        for (int j = 0; j < 3; j++)
        {
          xmin[j] = min(xmin[j],x[i3+j]);
          xmax[j] = max(xmax[j],x[i3+j]);
        }
      }

      elementBbox[6*adt_ic+0] = xmin[0];
      elementBbox[6*adt_ic+1] = xmin[1];
      elementBbox[6*adt_ic+2] = xmin[2];
      elementBbox[6*adt_ic+3] = xmax[0];
      elementBbox[6*adt_ic+4] = xmax[1];
      elementBbox[6*adt_ic+5] = xmax[2];

      elementList[adt_ic] = mb_ic;
    }

    // build the ADT now
    if (adt)
    {
      adt->clearData();
    }
    else
    {
      adt = new ADT[1];
    }

    int ndim = 6;
    adt->buildADT(ndim,ncells_adt,elementBbox.data());

#ifdef _GPU
    mb_d.dataToDevice(nDims,nnodes,ncells,ncells_adt,nsearch,nv,nc,elementList.data(),
                      elementBbox.data(),isearch.data(),xsearch.data(),myid);

    adt_d.copyADT(adt);
#endif
  }

  donorId.resize(nsearch);

//  Timer dtime("Device ADT Time: ");
//  Timer htime("Host ADT Time: ");
#ifdef _GPU
  //cudaDeviceSynchronize();
  //dtime.startTimer();
  searchADT(adt_d,mb_d);
  //cudaDeviceSynchronize();
  //dtime.stopTimer();

  rst.assign(mb_d.rst.data(), mb_d.rst.size());
  donorId.assign(mb_d.donorId.data(), mb_d.donorId.size());

  for (int i = 0; i < nsearch; i++)
  {
    if (donorId[i] > -1)
    {
      haveDonors = true;
      break;
    }
  }
#else
//  htime.startTimer();
//#pragma omp parallel for // CONTAINMENT CHECK MUST BE THREAD-SAFE
  for (int i = 0; i < nsearch; i++)
    adt->searchADT(this, &donorId[i], &xsearch[3*i], &rst[3*i]);

  donorCount = 0;
  for (int i = 0; i < nsearch; i++)
    if (donorId[i] > -1)
      donorCount++;
//  htime.stopTimer();
#endif
//  printf("%d: nsearch %d\n",myid,nsearch);
//  dtime.showTime(5);
//  htime.showTime(5);
}
