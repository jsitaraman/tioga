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
extern "C" 
{ 
  void fillHoleMap(int *holeMap, int ix[3],int isym);
}

/**
 * Create hole maps (structured cartesian maps of solid surfaces) for all grids
 * Using simplified 'vision space bins' concept from Sitaraman et al, JCP 2010
 * This routine is not efficient since it does mutiple global reduce ops
 * have to change it at a later date when there is more time to develop code
 */
void tioga::getHoleMap(void)
{
  // Get the bounding box of all wall boundary nodes on this rank
  double wbox[6];
  int existWall;
  int meshtag;

  if (holeMap)
  {
    for (int i = 0; i < nmesh; i++)
      if (holeMap[i].existWall)
      {
        free(holeMap[i].sam);
        free(holeMap[i].samLocal);
      }
    delete [] holeMap;
  }

  // Grid/Mesh ID, whether rank has wall nodes, and bbox of wall nodes if so
  mb->getWallBounds(&meshtag, &existWall, wbox);

  // Get 'nmesh' : number of grids (separate mesh tags) in system
  MPI_Allreduce(&meshtag, &nmesh, 1, MPI_INT, MPI_MAX, scomm);
  nmesh += 1;

  holeMap = new HOLEMAP[nmesh];

  std::vector<int> existHoleLocal(nmesh);
  std::vector<int> existHole(nmesh);

  existHoleLocal[meshtag] = existWall;

  MPI_Allreduce(existHoleLocal.data(),existHole.data(),nmesh,MPI_INT,MPI_MAX,scomm);

  for (int i = 0; i < nmesh; i++) holeMap[i].existWall = existHole[i];

  std::vector<double> bboxLocal(6*nmesh);
  std::vector<double> bboxGlobal(6*nmesh);

  for (int i = 0; i < 3*nmesh; i++) bboxLocal[i]          =  BIGVALUE;
  for (int i = 0; i < 3*nmesh; i++) bboxLocal[i+3*nmesh]  = -BIGVALUE;
  for (int i = 0; i < 3*nmesh; i++) bboxGlobal[i]         =  BIGVALUE;
  for (int i = 0; i < 3*nmesh; i++) bboxGlobal[i+3*nmesh] = -BIGVALUE;

  for (int i = 0;  i < 3; i++)
  {
    bboxLocal[3*meshtag+i] = wbox[i];
    bboxLocal[3*meshtag+i+3*nmesh] = wbox[i+3];
  }

  // Get the global bounding box info across all the partitions for all meshes
  MPI_Allreduce(&bboxLocal[0],        &bboxGlobal[0],        3*nmesh,MPI_DOUBLE,MPI_MIN,scomm);
  MPI_Allreduce(&(bboxLocal[3*nmesh]),&(bboxGlobal[3*nmesh]),3*nmesh,MPI_DOUBLE,MPI_MAX,scomm);

  // Find the bounding box for each mesh from the globally reduced data
  for (int i = 0; i < nmesh; i++)
  {
    if (holeMap[i].existWall)
    {
      double ds[3];
      for (int j = 0; j < 3; j++)
      {
        holeMap[i].extents[j] = bboxGlobal[3*i+j];
        holeMap[i].extents[j+3] = bboxGlobal[3*i+j+3*nmesh];
        ds[j] = holeMap[i].extents[j+3]-holeMap[i].extents[j];
      }
      double dsmax = max(max(ds[0], ds[1]), ds[2]);
      double dsmin = min(min(ds[0], ds[1]), ds[2]);
      double dx_avg = .5*(dsmax + dsmin);
      int NX_HOLEMAP = 75; // Originally 64
      double dsbox = dx_avg/NX_HOLEMAP; // Accounting somewhat for high aspect ratios

      for (int j = 0; j < 3; j++)
      {
        // Extend bounding box by 2/NX_HOLEMAP of max dimension in each direction
        holeMap[i].extents[j]   -= (2*dsbox);
        holeMap[i].extents[j+3] += (2*dsbox);
        // nx should end up equal to 68 for maximum dimension
        holeMap[i].nx[j] = floor(max((ds[j]+4*dsbox)/dsbox, 1.));
      }

      int bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1] * holeMap[i].nx[2];
      holeMap[i].sam = (int *)malloc(sizeof(int)*bufferSize);
      holeMap[i].samLocal = (int *)malloc(sizeof(int)*bufferSize);

      for (int j = 0; j < bufferSize; j++)
        holeMap[i].sam[j] = holeMap[i].samLocal[j] = 0;
    }
  }

  // mark the wall boundary cells in the holeMap
  if (holeMap[meshtag].existWall)
  {
    mb->markWallBoundary(holeMap[meshtag].samLocal,holeMap[meshtag].nx,holeMap[meshtag].extents);
  }

  // allreduce the holeMap of each mesh
  for (int i = 0; i < nmesh; i++)
  {
    if (holeMap[i].existWall)
    {
      int bufferSize = holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
      MPI_Allreduce(holeMap[i].samLocal,holeMap[i].sam,bufferSize,MPI_INT,MPI_MAX,scomm);
    }
  }

  // now fill the holeMap
  for (int i = 0; i < nmesh; i++)
    if (holeMap[i].existWall) fillHoleMap(holeMap[i].sam,holeMap[i].nx,isym);

  // output the hole map
  //this->outputHoleMap();
}

/**
 * Output the hole map to a tecplot compatible file
*/
void tioga::outputHoleMap(void)
{
  int i,k;
  int nnodes,ncells;
  int ns1,ns2;
  int ii,jj,kk,m;
  FILE *fp;
  double ds[3];
  char intstring[7];
  char fname[80];

  for(i=0;i<nmesh;i++)
    if (holeMap[i].existWall)
       {
	 sprintf(intstring,"%d",100000+i+100*myid);
	 sprintf(fname,"holeMap%s.dat",&(intstring[1]));
	 fp=fopen(fname,"w");
	 fprintf(fp,"TITLE =\"Tioga output\"\n");
	 fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
	 nnodes=(holeMap[i].nx[0]+1)*(holeMap[i].nx[1]+1)*(holeMap[i].nx[2]+1);
	 ncells=(holeMap[i].nx[0])*(holeMap[i].nx[1])*(holeMap[i].nx[2]);	 
	 fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,ncells);
	 fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");
	 for(k=0;k<3;k++) ds[k]=(holeMap[i].extents[k+3]-holeMap[i].extents[k])/(holeMap[i].nx[k]);
	 //
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",ii*ds[0]);
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",jj*ds[1]);
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",kk*ds[2]);
	 m=0;
	 for(kk=0;kk<holeMap[i].nx[2];kk++)
	   for(jj=0;jj<holeMap[i].nx[1];jj++)
	     for(ii=0;ii<holeMap[i].nx[0];ii++)
	       {
		 fprintf(fp,"%f\n",(double)holeMap[i].sam[m]);
		 m++;
	       }
	 
	 m=0;
         ns1=holeMap[i].nx[0]+1;
	 ns2=(holeMap[i].nx[1]+1)*ns1;
	 for(kk=0;kk<holeMap[i].nx[2];kk++)
	   for(jj=0;jj<holeMap[i].nx[1];jj++)
	     for(ii=0;ii<holeMap[i].nx[0];ii++)
	       {
		 m=kk*ns2+jj*ns1+ii+1;
		 fprintf(fp,"%d %d %d %d %d %d %d %d\n",m,m+1,m+1+ns1,m+ns1,
			 m+ns2,m+1+ns2,m+ns2+ns1+1,m+ns1+ns2);
	       }
       }
 fclose(fp);
}
	 
