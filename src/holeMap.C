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
  int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
			double vB[3][3],double xB[3],double dxB[3]);			   

};
# define NSAM 64
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
  double wbox[6];
  int existWall;
  int meshtag,maxtag;
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
  int itag,iflag,nb,anyval,recvmeshtag;
  int nsend,nrecv;
  int *proc2meshtagmap,*meshtag2procmap,*meshtagcft;
  int *commcount,*commcountGlobal;
  int *procListRecv,*procListSend,*blockmap;
  OBB wobb;
  PACKET *sndPack,*rcvPack;
  MPI_Request *mpirequest;
  MPI_Status *mpistatus;
  MPI_Status status2,status3;
  //
  proc2meshtagmap=(int *)malloc(sizeof(int)*pc->numprocs);
  commcount=(int *)malloc(sizeof(int)*pc->numprocs);
  meshtag2procmap=(int *)malloc(sizeof(int)*pc->numprocs);
  //
  mb->getWallBounds(&meshtag,&existWall,wbox);
  MPI_Allgather(&meshtag,1,MPI_INT,proc2meshtagmap,1,MPI_INT,scomm);
  maxtag=0;
  for(i=0;i<pc->numprocs;i++) maxtag=max(maxtag,proc2meshtagmap[i]);
  //
  meshtagcft=(int *)malloc(sizeof(int)*(maxtag+1));
  commcount=(int *)malloc(sizeof(int)*(maxtag+1));
  blockmap=(int *)malloc(sizeof(int)*maxtag);
  //
  // invert the map now
  //
  for(i=0;i<maxtag+1;i++) {
    meshtagcft[i]=0;
    commcount[i]=0;
  }
  for(i=0;i<pc->numprocs;i++)
      { 
      meshtagcft[proc2meshtagmap[i]]++;
      meshtag2procmap[i]=0;      
     }
  for(i=1;i<maxtag+1;i++)
      meshtagcft[i]=meshtagcft[i]+meshtagcft[i-1];
  for(i=0;i<pc->numprocs;i++)
    {
      itag=proc2meshtagmap[i]-1;
      meshtag2procmap[meshtagcft[itag]+commcount[itag]]=i;
      commcount[itag]++;
    }
  free(commcount);
  //
  commcount=(int *)malloc(sizeof(int)*pc->numprocs);
  commcountGlobal=(int *)malloc(sizeof(int)*pc->numprocs);
  procListRecv=(int *)malloc(sizeof(int)*pc->numprocs);
  //
  for(i=0;i<pc->numprocs;i++) {
    commcount[i]=0;
    commcountGlobal[i]=0;
  }
  //
  //
  // get the local bounding box
  //  
 if (holeMap) 
   {
     for(i=0;i<nmesh;i++)
       if (holeMap[i].existWall)
	 {
	   if (holeMap[i].sam) free(holeMap[i].sam);
	 }
     delete [] holeMap;
   }
 holeMap=new HOLEMAP[maxtag];
 for(i=0;i<maxtag;i++)
   holeMap[i].sam=NULL;
 //
 existHoleLocal=(int *)malloc(sizeof(int)*maxtag);
 existHole=(int *)malloc(sizeof(int)*maxtag);
 //
 for(i=0;i<maxtag;i++) existHole[i]=existHoleLocal[i]=0;
 //
 existHoleLocal[meshtag-1]=existWall;
 //
 MPI_Allreduce(existHoleLocal,existHole,maxtag,MPI_INT,MPI_MAX,scomm);
 //
 for(i=0;i<maxtag;i++) { holeMap[i].existWall=existHole[i];}
 //
 bboxLocal=(double *) malloc(sizeof(double)*6*maxtag);
 bboxGlobal=(double *) malloc(sizeof(double)*6*maxtag);
 //
 for(i=0;i<3*maxtag;i++) bboxLocal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxLocal[i+3*maxtag]=-BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i+3*maxtag]=-BIGVALUE;
 //
 for(i=0;i<3;i++)
   {
     bboxLocal[3*(meshtag-1)+i]=wbox[i];
     bboxLocal[3*(meshtag-1)+i+3*maxtag]=wbox[i+3];
   }
 //
 // get the global bounding box info across all the
 // partitions for all meshes
 //
 MPI_Allreduce(bboxLocal,bboxGlobal,3*maxtag,MPI_DOUBLE,MPI_MIN,scomm);
 MPI_Allreduce(&(bboxLocal[3*maxtag]),&(bboxGlobal[3*maxtag]),3*maxtag,MPI_DOUBLE,MPI_MAX,scomm);
 //
 // test bounding boxes against OBB for 
 // possible intersection
 //
 wobb.vec[0][0]=1;wobb.vec[0][1]=0;wobb.vec[0][2]=0;
 wobb.vec[1][0]=0;wobb.vec[1][1]=1;wobb.vec[1][2]=0;
 wobb.vec[2][0]=0;wobb.vec[2][1]=0;wobb.vec[2][2]=1;
 //
 // figure out who may have a wall box that
 // the local mesh block has intersection with
 // and set comm patterns for those.
 //
 nrecv=nb=0;
 for(i=0;i<maxtag;i++)
   {
     if (holeMap[i].existWall )
       {
	 for(j=0;j<3;j++)
	   {
	     holeMap[i].extents[j]=bboxGlobal[3*i+j];
	     holeMap[i].extents[j+3]=bboxGlobal[3*i+j+3*maxtag];
	     ds[j]=holeMap[i].extents[j+3]-holeMap[i].extents[j];
	   }	 
	 dsmax=max(ds[0],ds[1]);
	 dsmax=max(dsmax,ds[2]);	 
	 dsbox=dsmax/NSAM;	 
	 for(j=0;j<3;j++)
	   {
	     holeMap[i].extents[j]-=(2*dsbox);
	     holeMap[i].extents[j+3]+=(2*dsbox);
	     holeMap[i].nx[j]=floor(max((holeMap[i].extents[j+3]-holeMap[i].extents[j])/dsbox,1));
	     wobb.xc[j]=(holeMap[i].extents[j+3]+holeMap[i].extents[j])*0.5;
	     wobb.dxc[j]=(holeMap[i].extents[j+3]-holeMap[i].extents[j])*0.5;
	   }
	 if (i!=meshtag-1) { 	 	 
	 if ( obbIntersectCheck(mb->obb->vec,mb->obb->xc,mb->obb->dxc,
				wobb.vec,wobb.xc,wobb.dxc) ||
	      obbIntersectCheck(wobb.vec,wobb.xc,wobb.dxc,
				mb->obb->vec,mb->obb->xc,mb->obb->dxc))
	   {
	     for(j=meshtagcft[i];j<meshtagcft[i+1];j++)
	       {
		 commcount[j]++;
		 procListRecv[nrecv]=meshtag2procmap[j];
		 nrecv++;
	       }
	     blockmap[nb]=i;
	     nb++;
	   }
         }
       }
   }
 //
 // all reduce the comm map to produce the
 // sender side info
 //
 mpirequest=(MPI_Request *) malloc(sizeof(MPI_Request)*nrecv);
 mpistatus=(MPI_Status *) malloc(sizeof(MPI_Status)*nrecv);
 MPI_Allreduce(commcount,commcountGlobal,pc->numprocs,MPI_INT,MPI_SUM,scomm);
 for(i=0;i<nrecv;i++)
   MPI_Isend(&i,1,MPI_INT,procListRecv[i],1,scomm,&mpirequest[i]);
 procListSend=(int *)malloc(sizeof(int)*commcountGlobal[pc->myid]);
 nsend=0;
 while(nsend < commcountGlobal[pc->myid])
   {
     MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,scomm,&iflag,&status2);
     if (iflag) {       
       MPI_Recv(&anyval,1,MPI_INT,status2.MPI_SOURCE,1,scomm,&status3);
       procListSend[nsend]=status2.MPI_SOURCE;
       nsend++;
     }
   }
 MPI_Waitall(nrecv,mpirequest,mpistatus);
 free(mpirequest);
 free(mpistatus);
 //
 pc->setMap(nsend,nrecv,procListSend,procListRecv);
 //
 i=meshtag-1;
 bufferSize=holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
 holeMap[i].samLocal=(int *)malloc(sizeof(int)*bufferSize);	  
 for(j=0;j<bufferSize;j++) holeMap[i].samLocal[j]=0;
 //
 for(i=0;i<nb;i++)
   {
     bufferSize=holeMap[blockmap[i]].nx[0]*holeMap[blockmap[i]].nx[1]*holeMap[blockmap[i]].nx[2];
     holeMap[blockmap[i]].sam=(int *)malloc(sizeof(int)*bufferSize);
     for(j=0;j<bufferSize;j++) holeMap[blockmap[i]].sam[j]=0;
   }
 //
 // mark the wall boundary cells in the holeMap
 //
 if (holeMap[meshtag-1].existWall) {
 mb->markWallBoundary(holeMap[meshtag-1].samLocal,holeMap[meshtag-1].nx,holeMap[meshtag-1].extents);
 }
 //
 // now send the information across of the
 // stuff we painted
 // this can be compacted to not send zeros
 // next level of optimization TODO
 //
 sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
 rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
 pc->initPackets(sndPack,rcvPack);
 //
 for(i=0;i<nsend;i++)
   {
     sndPack[i].nints=holeMap[meshtag-1].nx[0]*holeMap[meshtag-1].nx[1]*holeMap[meshtag-1].nx[2];
     sndPack[i].intData=(int *)malloc(sndPack[i].nints*sizeof(int));
     for(j=0;j<sndPack[i].nints;j++)
       sndPack[i].intData[j]=holeMap[meshtag-1].samLocal[j];
   }
 pc->sendRecvPackets(sndPack,rcvPack);
 //
 for(i=0;i<nrecv;i++)
   {
     recvmeshtag=proc2meshtagmap[procListRecv[i]];
     for(j=0;j<rcvPack[i].nints;j++)
       holeMap[recvmeshtag-1].sam[j]+=rcvPack[i].intData[j];
   }
 //
 pc->clearPackets(sndPack,rcvPack);
 free(sndPack);
 free(rcvPack);
 //
 // free sam local
 //
 i=meshtag-1;
 if (holeMap[i].existWall) free(holeMap[i].samLocal);
 //
 // set the global number of meshes to maxtag
 //
 nmesh=maxtag;
 //
 // now fill the holeMap
 //
 for(j=0;j<nb;j++)
   {
     i=blockmap[j];
     if (holeMap[i].existWall) fillHoleMap(holeMap[i].sam,holeMap[i].nx,isym);
   }
 //
 // output the hole map
 //
 //this->outputHoleMap();
 //
 // free local memory
 //
 free(existHoleLocal);
 free(existHole);
 free(bboxLocal);
 free(bboxGlobal);
 //
 free(proc2meshtagmap);
 free(meshtag2procmap);
 free(meshtagcft);
 free(commcount);
 free(commcountGlobal);
 free(procListSend);
 free(procListRecv);
 free(blockmap);
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
	 
