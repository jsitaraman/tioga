// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include <limits>
#include <vector>
#include <array>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
extern "C" 
{ 
  void fillHoleMap(int *holeMap, int ix[3],int isym);

};

/**
 * Create hole maps for all grids
 * this routine is not efficient
 * since it does mutiple global reduce ops
 * have to change it at a later date when
 * there is more time to develop code
 */
void tioga::getHoleMap(void)
{
  int i,j,k,m;
  int ii,jj,kk;
  // double wbox[6];
  std::vector<std::array<double,6>> wbox(nblocks);
  //int existWall;
  std::vector<int> existWall(nblocks);
  int meshtag,maxtag, mtagtmp;
  int *existHoleLocal;
  int *existHole;
  double *bboxLocal;
  double *bboxGlobal;
  double ds[3],dsmax,dsbox;
  int bufferSize;
  FILE *fp;
  char fname[80];
  char intstring[7];
 //
 // get the local bounding box
 //
  meshtag = -BIGINT; //std::numeric_limits<int>::lowest();
  for (int i=0; i<nblocks; i++) {
    auto& mb = mblocks[i];
    mb->getWallBounds(&mtagtmp,&existWall[i],wbox[i].data());
    if (mtagtmp > meshtag) meshtag = mtagtmp;
  }
  MPI_Allreduce(&meshtag,&maxtag,1,MPI_INT,MPI_MAX,scomm);
 //
 if (holeMap) 
   {
     for(i=0;i<nmesh;i++)
       if (holeMap[i].existWall) TIOGA_FREE(holeMap[i].sam);
     delete [] holeMap;
   }
 holeMap=new HOLEMAP[maxtag];
 //
 existHoleLocal=(int *)malloc(sizeof(int)*maxtag);
 existHole=(int *)malloc(sizeof(int)*maxtag);
 //
 for(i=0;i<maxtag;i++) existHole[i]=existHoleLocal[i]=0;
 //
 for (int i=0; i<nblocks; i++) {
   existHoleLocal[mtags[i]-1]=existWall[i];
 }
 //
 MPI_Allreduce(existHoleLocal,existHole,maxtag,MPI_INT,MPI_MAX,scomm);
 //
 for(i=0;i<maxtag;i++) holeMap[i].existWall=existHole[i];
 //
 bboxLocal=(double *) malloc(sizeof(double)*6*maxtag);
 bboxGlobal=(double *) malloc(sizeof(double)*6*maxtag);
 //
 for(i=0;i<3*maxtag;i++) bboxLocal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxLocal[i+3*maxtag]=-BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i+3*maxtag]=-BIGVALUE;

 //
 for (int n=0; n<nblocks; n++) {
   meshtag = mtags[n];
   for (i = 0; i < 3; i++) {
     bboxLocal[3 * (meshtag - 1) + i] = wbox[n][i];
     bboxLocal[3 * (meshtag - 1) + i + 3 * maxtag] = wbox[n][i + 3];
   }
 }
 //
 // get the global bounding box info across all the
 // partitions for all meshes
 //
 MPI_Allreduce(bboxLocal, bboxGlobal, 3 * maxtag, MPI_DOUBLE, MPI_MIN, scomm);
 MPI_Allreduce(
   &(bboxLocal[3 * maxtag]), &(bboxGlobal[3 * maxtag]), 3 * maxtag, MPI_DOUBLE,
   MPI_MAX, scomm);
 //
 // find the bounding box for each mesh
 // from the globally reduced data
 //
 for (i = 0; i < maxtag; i++) {
   if (holeMap[i].existWall) {
     for (j = 0; j < 3; j++) {
       holeMap[i].extents[j] = bboxGlobal[3 * i + j];
       holeMap[i].extents[j + 3] = bboxGlobal[3 * i + j + 3 * maxtag];
       ds[j] = holeMap[i].extents[j + 3] - holeMap[i].extents[j];
     }
     dsmax = TIOGA_MAX(ds[0], ds[1]);
     dsmax = TIOGA_MAX(dsmax, ds[2]);
     dsbox = dsmax / HOLEMAPSIZE;

     for (j = 0; j < 3; j++) {
       holeMap[i].extents[j] -= (2 * dsbox);
       holeMap[i].extents[j + 3] += (2 * dsbox);
       holeMap[i].nx[j] = floor(
         TIOGA_MAX((holeMap[i].extents[j + 3] - holeMap[i].extents[j]) / dsbox, 1));
     }
     bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1] * holeMap[i].nx[2];
     holeMap[i].sam = (int*)malloc(sizeof(int) * bufferSize);
     holeMap[i].samLocal = (int*)malloc(sizeof(int) * bufferSize);
     for (j = 0; j < bufferSize; j++)
       holeMap[i].sam[j] = holeMap[i].samLocal[j] = 0;
   }
 }
 //
 // mark the wall boundary cells in the holeMap
 //
 for (int ib=0;ib<nblocks;ib++) {
   auto& mb = mblocks[ib];
   meshtag = mb->getMeshTag();
   if (holeMap[meshtag - 1].existWall) {
    mb->markWallBoundary(
      holeMap[meshtag - 1].samLocal, holeMap[meshtag - 1].nx,
      holeMap[meshtag - 1].extents);
   }
 }
 //
 // allreduce the holeMap of each mesh
 //
 for(i=0;i<maxtag;i++)
   {
    if (holeMap[i].existWall) 
     {
      bufferSize=holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
      MPI_Allreduce(holeMap[i].samLocal,holeMap[i].sam,bufferSize,MPI_INT,MPI_MAX,scomm);
     }
   }
 //
 for(i=0;i<maxtag;i++)
   if (holeMap[i].existWall) TIOGA_FREE(holeMap[i].samLocal);
 //
 // set the global number of meshes to maxtag
 //
 nmesh=maxtag;
 //
 // now fill the holeMap
 //
 for(i=0;i<maxtag;i++)
   if (holeMap[i].existWall) fillHoleMap(holeMap[i].sam,holeMap[i].nx,isym);
 //
 // output the hole map
 //
 //this->outputHoleMap();
 //
 // free local memory
 //
 TIOGA_FREE(existHoleLocal);
 TIOGA_FREE(existHole);
 TIOGA_FREE(bboxLocal);
 TIOGA_FREE(bboxGlobal);
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
	 
