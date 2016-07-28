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
#include "tioga.h"

extern "C"{
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
			double vB[3][3],double xB[3],double dxB[3]);			   
}

void tioga::exchangeBoxes(void)
{
  int* alltags=(int *)malloc(sizeof(int)*numprocs);
  MPI_Allgather(&mytag, 1, MPI_INT, alltags,1,MPI_INT,scomm);

  // count number of other processors to communicate to
  // in overset grid scenario, usually you do not communicate
  // to mesh blocks that carry the same mesh tag (i.e. you don't
  // talk to your sister partitions)
  int nsend = 0;
  int nrecv = 0;
  for (int i = 0; i < numprocs; i++) if (alltags[i] != mytag) nsend++;

  // In general we communicate forward
  // and backward, separate lists are maintained for
  // flexibility
  nrecv = nsend;
  int* sndMap=(int *)malloc(sizeof(int)*nsend);
  int* rcvMap=(int *)malloc(sizeof(int)*nrecv);

  for (int i = 0, m = 0; i < numprocs; i++)
    if (alltags[i]!=mytag)
    {
      sndMap[m]=rcvMap[m]=i;
      m++;
    }

  pc->setMap(nsend,nrecv,sndMap,rcvMap);

  PACKET *sndPack = (PACKET *)malloc(sizeof(PACKET)*nsend);
  PACKET *rcvPack = (PACKET *)malloc(sizeof(PACKET)*nrecv);

  for (int k = 0; k < nsend; k++)
  {
    sndPack[k].nints = 0;
    sndPack[k].nreals = 15;
    sndPack[k].realData = (double *)malloc(sizeof(double)*sndPack[k].nreals);
    int m = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        sndPack[k].realData[m++]=mb->obb->vec[i][j];
    for (int i = 0; i < 3; i++)
      sndPack[k].realData[m++]=mb->obb->xc[i];
    for (int i = 0; i < 3; i++)
      sndPack[k].realData[m++]=mb->obb->dxc[i];
  }

  pc->sendRecvPackets(sndPack,rcvPack);

  free(obblist);
  obblist = (OBB *)malloc(sizeof(OBB)*nrecv);

  for (int k = 0; k < nrecv; k++)
  {
    int m = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        obblist[k].vec[i][j]=rcvPack[k].realData[m++];
    for (int i = 0; i < 3; i++)
      obblist[k].xc[i]=rcvPack[k].realData[m++];
    for (int i = 0; i < 3; i++)
      obblist[k].dxc[i]=rcvPack[k].realData[m++];
  }

  int m = 0;
  for (int k = 0; k < nrecv; k++)
  {
    if ( obbIntersectCheck(mb->obb->vec,mb->obb->xc,mb->obb->dxc,
                           obblist[k].vec,obblist[k].xc,obblist[k].dxc) ||
         obbIntersectCheck(obblist[k].vec,obblist[k].xc,obblist[k].dxc,
                           mb->obb->vec,mb->obb->xc,mb->obb->dxc))
    {
      rcvMap[m]=sndMap[k];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          obblist[m].vec[i][j]=obblist[k].vec[i][j];
      for (int i = 0; i < 3; i++)
        obblist[m].xc[i]=obblist[k].xc[i];
      for (int i = 0; i < 3; i++)
        obblist[m].dxc[i]=obblist[k].dxc[i];
      m++;
    }
  }
  nsend = nrecv = m;

  for (int i = 0; i < nsend; i++) sndMap[i] = rcvMap[i];

  // clear packets before nsend and nrecv are modified in pc->setMap
  pc->clearPackets(sndPack,rcvPack);

  pc->setMap(nsend,nrecv,sndMap,rcvMap);  

  // free local memory
  free(alltags);
  free(sndMap);
  free(rcvMap);
  free(sndPack);
  free(rcvPack);    
}
