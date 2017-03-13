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
    /// TODO: worry about updating ADT based on changing OBB of search points
    adt->setTransform(Rmat.data(), offset, nDims);
#ifdef _GPU
    adt_d.setTransform(Rmat.data(), offset, nDims);
    mb_d.setTransform(Rmat.data(), offset, nDims);
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

    if (rrot)
    {
      /// HACK - adding a 20% safety factor to ADT bounds in case of moving grids
      for (int i = 0; i < nDims; i++)
        obq.dxc[i] *= 1.5;
    }

    // find all the cells that may have intersections with
    // the OBB
    std::vector<int> icell(ncells, -1);

    int iptr = -1;
    ncells_adt = 0;
    int p = 0;
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
        if (fabs(xd[0]) <= (dxc[0]+obq.dxc[0]) &&
            fabs(xd[1]) <= (dxc[1]+obq.dxc[1]) &&
            fabs(xd[2]) <= (dxc[2]+obq.dxc[2]))
        {
          // create a LIFO stack
          // with all the cells that
          // have bounding box intersection with
          // the QP bounding box
          icell[p] = iptr;
          iptr = p;
          ncells_adt++;
        }
        p++;
      }
    }

    // now find the axis aligned bounding box
    // of each cell in the LIFO stack to build the
    // ADT
    elementBbox.resize(ncells_adt*6);
    elementList.resize(ncells_adt);

    int k = iptr;
    int l = 0;
    p = 0;
    while (k != -1)
    {
      int ic = -1;
      int cellindex = k;
      int isum = 0;
      int n;
      for (n = 0; n < ntypes; n++)
      {
        isum += nc[n];
        if (cellindex < isum)
        {
          ic = cellindex-(isum-nc[n]);
          break;
        }
      }

      int nvert = nv[n];
      xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
      xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
      for (int m = 0; m < nvert; m++)
      {
        int i3 = 3*(vconn[n][nvert*ic+m]-BASE);
        for (int j = 0; j < 3; j++)
        {
          xmin[j] = min(xmin[j],x[i3+j]);
          xmax[j] = max(xmax[j],x[i3+j]);
        }
      }

      elementBbox[l++] = xmin[0];
      elementBbox[l++] = xmin[1];
      elementBbox[l++] = xmin[2];
      elementBbox[l++] = xmax[0];
      elementBbox[l++] = xmax[1];
      elementBbox[l++] = xmax[2];

      elementList[p++] = k;

      k = icell[k];
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
printf("Hey! You shouldn't be here for this case!\n");
#ifdef _GPU
    mb_d.dataToDevice(nDims,nnodes,ncells,ncells_adt,nsearch,nv,nc,elementList.data(),
                      elementBbox.data(),isearch.data(),xsearch.data());

    adt_d.copyADT(adt);
#endif
  }

  donorId.resize(nsearch);

  Timer dtime("Device ADT Time: ");
  Timer htime("Host ADT Time: ");
#ifdef _GPU

  //cudaDeviceSynchronize();
  //dtime.startTimer();
  searchADT(adt_d,mb_d);
  //cudaDeviceSynchronize();
  //dtime.stopTimer();

  rst.assign(mb_d.rst.data(), mb_d.rst.size());
  donorId.assign(mb_d.donorId.data(), mb_d.donorId.size());
#else
  htime.startTimer();
  donorCount = 0;
  ipoint = 0;
  for (int i = 0; i < nsearch; i++)
  {
    adt->searchADT(this, &(donorId[i]), &(xsearch[3*i]));

    if (donorId[i] > -1)
      donorCount++;

    ipoint += 3;
  }
  htime.stopTimer();

  printf("%d: nsearch %d\n",myid,nsearch);
  dtime.showTime(5);
  htime.showTime(5);
#endif
}
