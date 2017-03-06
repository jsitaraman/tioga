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
  int i,j,k,m,k1;
  double *obbdata,*obbdatalocal;
  int *sndMap;
  int *rcvMap;
  int *displs;
  int nsend;
  int nrecv;
  int overlap_present;
  PACKET *sndPack,*rcvPack;
  //
  obbdata=(double *)malloc(sizeof(double)*16*numprocs);
  obbdatalocal=(double *)malloc(sizeof(double)*16);
  m=0;
  obbdatalocal[m++]=(double)mytag;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      obbdatalocal[m++]=mb->obb->vec[i][j];
  for(i=0;i<3;i++)
    obbdatalocal[m++]=mb->obb->xc[i];
  for(i=0;i<3;i++)
    obbdatalocal[m++]=mb->obb->dxc[i];
  //
  MPI_Allgather(obbdatalocal, 16, MPI_DOUBLE, obbdata,16,MPI_DOUBLE,scomm);
  nrecv=0;
  for(k=0;k<numprocs;k++) 	
    if ((int)obbdata[16*k] != mytag) nrecv++;
  //
  if (obblist) free(obblist);
  obblist=(OBB *) malloc(sizeof(OBB)*nrecv);
  sndMap=(int *)malloc(sizeof(int)*nrecv);
  rcvMap=(int *)malloc(sizeof(int)*nrecv);
  //
  k1=0;
  for(k=0;k<numprocs;k++)
    {
      if ((int)obbdata[16*k] != mytag)
	{	 
	  m=1;      
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      {
		obblist[k1].vec[i][j]=obbdata[16*k+m];
		m++;
	      }
	  for(i=0;i<3;i++)
	    {
	      obblist[k1].xc[i]=obbdata[16*k+m];
	      m++;
	    }	      
	  for(i=0;i<3;i++)
	    {
	      obblist[k1].dxc[i]=obbdata[16*k+m];
	      m++;
	    }
          sndMap[k1]=k;
	  k1++;
	}
    }
  //
  m=0;
  for(k=0;k<nrecv;k++)
    {
      if ( obbIntersectCheck(mb->obb->vec,mb->obb->xc,mb->obb->dxc,
			     obblist[k].vec,obblist[k].xc,obblist[k].dxc) ||
	   obbIntersectCheck(obblist[k].vec,obblist[k].xc,obblist[k].dxc,
			     mb->obb->vec,mb->obb->xc,mb->obb->dxc)) 
	{
	  if ((int)obbdata[16*sndMap[k]] < 0 || mytag < 0) 
	    {
	      mb->check_intersect_p4est(&sndMap[k],&overlap_present);	      
	    }
	  else
	    {
	      overlap_present=1;
	    }
	  if (overlap_present==1)
	    {
	      rcvMap[m]=sndMap[k];
	      for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		  obblist[m].vec[i][j]=obblist[k].vec[i][j];
	      for(i=0;i<3;i++)
		obblist[m].xc[i]=obblist[k].xc[i];
	      for(i=0;i<3;i++)
		obblist[m].dxc[i]=obblist[k].dxc[i];
	      m++;
	    }
	}
    }
  nsend=nrecv=m;
  for(i=0;i<nsend;i++) sndMap[i]=rcvMap[i];
  //printf("%d %d %d\n",myid,nsend,nrecv);
  //
  pc->setMap(nsend,nrecv,sndMap,rcvMap);  
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  pc->initPackets(sndPack,rcvPack);
  for(k=0;k<nsend;k++) {
   sndPack[k].nints=0;
   sndPack[k].nreals=6;
   sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
   mb->getReducedOBB(&obblist[k],sndPack[k].realData);
  }
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  for(k=0;k<nrecv;k++) {
   for(j=0;j<3;j++) obblist[k].xc[j]=rcvPack[k].realData[j];
   for(j=0;j<3;j++) obblist[k].dxc[j]=rcvPack[k].realData[j+3];
   }  
  pc->clearPackets(sndPack,rcvPack);
  //
  // Free local memory
  //
  free(sndMap);
  free(rcvMap);
  free(sndPack);
  free(rcvPack);
  free(obbdata);
  free(obbdatalocal);
}
