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

extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  double computeCellVolume(double xv[8][3],int nvert);
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);  
}

void MeshBlock::setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
			int *wbcnodei,int *obcnodei,
                        int ntypesi,int *nvi,int *nci,int **vconni,
                        uint64_t* cell_gid)
{
  int i;
  //
  // set internal pointers
  //
  meshtag=btag;
  nnodes=nnodesi;
  x=xyzi;
  iblank=ibli;
  nwbc=nwbci;
  nobc=nobci;
  wbcnode=wbcnodei;
  obcnode=obcnodei;
  //
  ntypes=ntypesi;
  //
  nv=nvi;
  nc=nci;
  vconn=vconni;
  cellGID = cell_gid;
  //
  //tracei(nnodes);
  //for(i=0;i<ntypes;i++) tracei(nc[i]);
  ncells=0;
  for(i=0;i<ntypes;i++) ncells+=nc[i];
}

void MeshBlock::preprocess(void)
{
  int i;
  //
  // set all iblanks = 1
  //
  for(i=0;i<nnodes;i++) iblank[i]=1;
  //
  // find oriented bounding boxes
  //
  if (obb) free(obb);
  obb=(OBB *) malloc(sizeof(OBB));
  findOBB(x,obb->xc,obb->dxc,obb->vec,nnodes);
  tagBoundary();
}

void MeshBlock::tagBoundary(void)
{
  int i,j,k,n,m,ii;
  int itag;
  int inode[8];
  double xv[8][3];
  double vol;
  int *iflag;
  int nvert,i3;
  FILE *fp;
  char intstring[7];
  char fname[80];  
  

  //
  // do this only once
  // i.e. when the meshblock is first 
  // initialized, cellRes would be NULL in this case
  //

      if(cellRes) free(cellRes);
      if(nodeRes) free(nodeRes);
    
      //
      cellRes=(double *) malloc(sizeof(double)*ncells);
      nodeRes=(double *) malloc(sizeof(double)*nnodes);
      //
      // this is a local array
      //
      iflag=(int *)malloc(sizeof(int)*nnodes);
      // 
      for(i=0;i<nnodes;i++) iflag[i]=0;
      //
      if (userSpecifiedNodeRes ==NULL && userSpecifiedCellRes ==NULL)
	{
	  for(i=0;i<nnodes;i++) iflag[i]=0;
	  for(i=0;i<nnodes;i++) nodeRes[i]=0.0;
	  //
	  k=0;
	  for(n=0;n<ntypes;n++)
	    {
	      nvert=nv[n];
	      for(i=0;i<nc[n];i++)
		{
		  for(m=0;m<nvert;m++)
		    {
		      inode[m]=vconn[n][nvert*i+m]-BASE;
		      i3=3*inode[m];
		      for(j=0;j<3;j++)
			xv[m][j]=x[i3+j];
		    }
		  vol=computeCellVolume(xv,nvert);
		  cellRes[k++]=(vol*resolutionScale);
		  for(m=0;m<nvert;m++)
		    {
		      iflag[inode[m]]++;
		      nodeRes[inode[m]]+=(vol*resolutionScale);
		    }	      
		}
	    }
	}
      else
	{
	  k=0;
	  for(n=0;n<ntypes;n++)
	    {
	      for(i=0;i<nc[n];i++)		
		{
		  cellRes[k]=userSpecifiedCellRes[k];
		  k++;
		}
	    }
	  for(k=0;k<nnodes;k++) nodeRes[k]=userSpecifiedNodeRes[k];
	}	  
      //
      // compute nodal resolution as the average of 
      // all the cells associated with it. This takes care
      // of partition boundaries as well.
      //
      for(i=0;i<nnodes;i++)
	{
          if (iflag[i]!=0)  nodeRes[i]/=iflag[i];
	  iflag[i]=0;
	}
      //
      // now tag the boundary nodes
      // reuse the iflag array
      //
      //tracei(nobc);
      for(i=0;i<nobc;i++)
        { 
          ii=(obcnode[i]-BASE);
	  iflag[(obcnode[i]-BASE)]=1;
        }
      //
      // now tag all the nodes of boundary cells
      // to be mandatory receptors
      //
      for(n=0;n<ntypes;n++)
	{
	  nvert=nv[n];
	  for(i=0;i<nc[n];i++)
	    {
	      itag=0;
	      for(m=0;m<nvert;m++)
		{
		  inode[m]=vconn[n][nvert*i+m]-BASE;
		  if (iflag[inode[m]]) itag=1;
		}
	      if (itag)
		{
		  for(m=0;m<nvert;m++)
		    {
		      //iflag[inode[m]]=1;
		      nodeRes[inode[m]]=BIGVALUE;
		    }
		}
	    }
	}
      /*
      sprintf(intstring,"%d",100000+myid);
      sprintf(fname,"nodeRes%s.dat",&(intstring[1]));
      fp=fopen(fname,"w");
      for(i=0;i<nnodes;i++)
       {
        if (nodeRes[i]==BIGVALUE) {
         fprintf(fp,"%e %e %e\n",x[3*i],x[3*i+1],x[3*i+2]);
         }
       }
      fclose(fp);
      */
      //
      // now tag all the cells which have 
      // mandatory receptors as nodes as not acceptable
      // donors
      k=0;
      for(n=0;n<ntypes;n++)
	{
	  nvert=nv[n];
	  for(i=0;i<nc[n];i++)
	    {
	      for(m=0;m<nvert;m++)
		{
		  inode[m]=vconn[n][nvert*i+m]-BASE;
		  if (nodeRes[inode[m]]==BIGVALUE) //(iflag[inode[m]]) 
		    {
		      cellRes[k]=BIGVALUE;
		      break;
		    }
		}
	      k++;
	    }
	}
      free(iflag);
}

void MeshBlock::writeGridFile(int bid)
{
  char fname[80];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"part%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
	  ncells);
  for(i=0;i<nnodes;i++)
    {
      fprintf(fp,"%.14e %.14e %.14e %d\n",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
    }

  ba=1-BASE;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
	  else if (nvert==8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}

void MeshBlock::writeCellFile(int bid)
{
  char fname[80];
  char qstr[2];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"cell%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\",\"IBLANK_CELL\" ");
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,
	  ncells);
  fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=NODAL,5=CELLCENTERED)\n");
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i+1]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i+2]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%d\n",iblank[i]);
  for(i=0;i<ncells;i++) fprintf(fp,"%d\n",iblank_cell[i]);
  ba=1-BASE;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
	  else if (nvert==8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}

void MeshBlock::writeFlowFile(int bid,double *q,int nvar,int type)
{
  char fname[80];
  char qstr[2];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"flow%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\" ");
  for(i=0;i<nvar;i++)
    {
      sprintf(qstr,"Q%d",i);
      fprintf(fp,"\"%s\",",qstr);
    }
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
	  ncells);

  if (type==0)
    {
      for(i=0;i<nnodes;i++)
	{
	  fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
	  for(j=0;j<nvar;j++)
	    fprintf(fp,"%lf ",q[i*nvar+j]);      
	  //for(j=0;j<nvar;j++)
	  //  fprintf(fp,"%lf ", x[3*i]+x[3*i+1]+x[3*i+2]);
          fprintf(fp,"\n");
	}
    }
  else
    {
      for(i=0;i<nnodes;i++)
        {
          fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
          for(j=0;j<nvar;j++)
            fprintf(fp,"%lf ",q[j*nnodes+i]);
          fprintf(fp,"\n");
        }
    }
  ba=1-BASE;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
	  else if (nvert==8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}
  
void MeshBlock::getWallBounds(int *mtag,int *existWall, double wbox[6])
{
  int i,j,i3;
  int inode;

  *mtag=meshtag+(1-BASE); 
  if (nwbc <=0) {
    *existWall=0;
    for(i=0;i<6;i++) wbox[i]=0;
    return;
  }

  *existWall=1;
  wbox[0]=wbox[1]=wbox[2]=BIGVALUE;
  wbox[3]=wbox[4]=wbox[5]=-BIGVALUE;

  for(i=0;i<nwbc;i++)
    {
      inode=wbcnode[i]-BASE;
      i3=3*inode;
      for(j=0;j<3;j++)
	{
	  wbox[j]=min(wbox[j],x[i3+j]);
	  wbox[j+3]=max(wbox[j+3],x[i3+j]);
	}
    }
  
}
  
void MeshBlock::markWallBoundary(int *sam,int nx[3],double extents[6])
{
  int i,j,k,m,n;
  int nvert;
  int ii,jj,kk,mm;
  int i3,iv;
  int *iflag;
  int *inode;
  char intstring[7];
  char fname[80];
  double ds[3];
  double xv;
  int imin[3];
  int imax[3];
  FILE *fp;
  //
  iflag=(int *)malloc(sizeof(int)*ncells);
  inode=(int *) malloc(sizeof(int)*nnodes);
  ///
  //sprintf(intstring,"%d",100000+myid);
  //sprintf(fname,"wbc%s.dat",&(intstring[1]));
  //fp=fopen(fname,"w");
  for(i=0;i<ncells;i++) iflag[i]=0;
  for(i=0;i<nnodes;i++) inode[i]=0;
  //
  for(i=0;i<nwbc;i++)
   {
    ii=wbcnode[i]-BASE;
    //fprintf(fp,"%e %e %e\n",x[3*ii],x[3*ii+1],x[3*ii+2]);
    inode[ii]=1;
   }
  //fclose(fp);
  //
  // mark wall boundary cells
  //
  m=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  for(j=0;j<nvert;j++)
	    {
	      ii=vconn[n][nvert*i+j]-BASE;
	      if (inode[ii]==1)
		{
		  iflag[m]=1;
		  break;
		}
	    }
	  m++;
	}
    }
  //
  // find delta's in each directions
  //
  for(k=0;k<3;k++) ds[k]=(extents[k+3]-extents[k])/nx[k];
  //
  // mark sam cells with wall boundary cells now
  //
   m=0;
   for(n=0;n<ntypes;n++)
     {
       nvert=nv[n];
       for(i=0;i<nc[n];i++)
 	{
	  if (iflag[m]==1) 
	    {
	      //
	      // find the index bounds of each wall boundary cell
	      // bounding box
	      //
	      imin[0]=imin[1]=imin[2]=BIGINT;
	      imax[0]=imax[1]=imax[2]=-BIGINT;
	      for(j=0;j<nvert;j++)
		{
		  i3=3*(vconn[n][nvert*i+j]-BASE);
		  for(k=0;k<3;k++)
		    {
		      xv=x[i3+k];
		      iv=floor((xv-extents[k])/ds[k]);
		      imin[k]=min(imin[k],iv);
		      imax[k]=max(imax[k],iv);
		    }
		}
	     for(j=0;j<3;j++)
              {
	       imin[j]=max(imin[j],0);
               imax[j]=min(imax[j],nx[j]-1);
              }
	      //
	      // mark sam to 1
	      //
	      for(kk=imin[2];kk<imax[2]+1;kk++)
	        for(jj=imin[1];jj<imax[1]+1;jj++)
		  for (ii=imin[0];ii<imax[0]+1;ii++)
		   {
		    mm=kk*nx[1]*nx[0]+jj*nx[0]+ii;
		     sam[mm]=2;
		   }
	    }
	  m++;
	}
     }
   free(iflag);
   free(inode);
}

void MeshBlock::getReducedOBB(OBB *obc,double *realData) 
{
  int i,j,k,m,n,i3;
  int nvert;
  bool iflag;
  double bbox[6],xd[3];

  for(j=0;j<3;j++)
    {
      realData[j]=BIGVALUE;
      realData[j+3]=-BIGVALUE;
    }
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  bbox[0]=bbox[1]=bbox[2]=BIGVALUE;
	  bbox[3]=bbox[4]=bbox[5]=-BIGVALUE;

	  for(m=0;m<nvert;m++)
	    {
	      i3=3*(vconn[n][nvert*i+m]-BASE);
	      for(j=0;j<3;j++) xd[j]=0;
	      for(j=0;j<3;j++)
		for(k=0;k<3;k++)
		  xd[j]+=(x[i3+k]-obc->xc[k])*obc->vec[j][k];
	      for(j=0;j<3;j++) bbox[j]=min(bbox[j],xd[j]);
	      for(j=0;j<3;j++) bbox[j+3]=max(bbox[j+3],xd[j]);
	    }
	  iflag=0;
	  for(j=0;j<3;j++) iflag=(iflag || (bbox[j] > obc->dxc[j]));
	  if (iflag) continue;
	  iflag=0;
	  for(j=0;j<3;j++) iflag=(iflag || (bbox[j+3] < -obc->dxc[j]));
	  if (iflag) continue;
	  for (m=0;m<nvert;m++)
	    {
	      i3=3*(vconn[n][nvert*i+m]-BASE);
	      for(j=0;j<3;j++) xd[j]=0;
	      for(j=0;j<3;j++)
		for(k=0;k<3;k++)
		  xd[j]+=(x[i3+k]-obb->xc[k])*obb->vec[j][k];
	      for(j=0;j<3;j++) realData[j]=min(realData[j],xd[j]);
	      for(j=0;j<3;j++) realData[j+3]=max(realData[j+3],xd[j]);
	    }
	}
    }
  for(j=0;j<6;j++) bbox[j]=realData[j];
  for(j=0;j<3;j++)
    {
      realData[j]=obb->xc[j];
      for(k=0;k<3;k++)
       realData[j]+=((bbox[k]+bbox[k+3])*0.5)*obb->vec[k][j];
      realData[j+3]=(bbox[j+3]-bbox[j])*0.51;
    }
  return;
}
	      
void MeshBlock::getQueryPoints(OBB *obc,
			       int *nints,int **intData,
			       int *nreals, double **realData)
{
  int i,j,k;
  int i3;
  double xd[3];
  int *inode;
  int iptr;
  int m;

  inode=(int *)malloc(sizeof(int)*nnodes);
  *nints=*nreals=0; 
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(x[i3+k]-obc->xc[k])*obc->vec[j][k];

      if (fabs(xd[0]) <= obc->dxc[0] &&
	  fabs(xd[1]) <= obc->dxc[1] &&
	  fabs(xd[2]) <= obc->dxc[2])
	{
	  inode[*nints]=i;
	  (*nints)++;
	  (*nreals)+=3;

	}
    }
  //
  (*intData)=(int *)malloc(sizeof(int)*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  //
  m=0;
  for(i=0;i<*nints;i++)
    {
      i3=3*inode[i];
      (*intData)[i]=inode[i];
      (*realData)[m++]=x[i3];
      (*realData)[m++]=x[i3+1];
      (*realData)[m++]=x[i3+2];
    }
  //
  free(inode);
}  
  
void MeshBlock::writeOBB(int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"box%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,
	  1);

  for(l=0;l<2;l++)
    {
      il=2*(l%2)-1;
      for(k=0;k<2;k++)
	{
	  ik=2*(k%2)-1;
	  for(j=0;j<2;j++)
	    {
	      ij=2*(j%2)-1;
	      xx[0]=xx[1]=xx[2]=0;
	      for(m=0;m<3;m++)
		xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
		  +ik*obb->vec[1][m]*obb->dxc[1]
		  +il*obb->vec[2][m]*obb->dxc[2];	      
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	    }
	}
    }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fprintf(fp,"%e %e %e\n",obb->xc[0],obb->xc[1],obb->xc[2]);
  for(k=0;k<3;k++)
   fprintf(fp,"%e %e %e\n",obb->vec[0][k],obb->vec[1][k],obb->vec[2][k]);
  fprintf(fp,"%e %e %e\n",obb->dxc[0],obb->dxc[1],obb->dxc[2]);
  fclose(fp);
}
//
// destructor that deallocates all the
// the dynamic objects inside
//
MeshBlock::~MeshBlock()
{
  int i;
  //
  // free all data that is owned by this MeshBlock
  // i.e not the pointers of the external code.
  //
  if (cellRes) free(cellRes);
  if (nodeRes) free(nodeRes);
  if (elementBbox) free(elementBbox);
  if (elementList) free(elementList);
  if (adt) delete[] adt;
  if (donorList) {
    for(i=0;i<nnodes;i++) deallocateLinkList(donorList[i]);
    free(donorList);
  }
  if (interpList) {
    for(i=0;i<interpListSize;i++)
      {
	if (interpList[i].inode) free(interpList[i].inode);
	if (interpList[i].weights) free(interpList[i].weights);
      }
    free(interpList);
  }
  return;
  if (interpList2) {
    for(i=0;i<interp2ListSize;i++)
      {
	if (interpList2[i].inode) free(interpList2[i].inode);
	if (interpList2[i].weights) free(interpList2[i].weights);
      }
    free(interpList2);
  }
  if (interpListCart) {
    for(i=0;i<interpListCartSize;i++)
      {
        if (interpListCart[i].inode) free(interpListCart[i].inode);
 	if (interpListCart[i].weights) free(interpListCart[i].weights);
      }
    free(interpListCart);
  }
  if (!ihigh) {
   if (iblank_cell) free(iblank_cell);
  }
  if (obb) free(obb);
  if (isearch) free(isearch);
  if (xsearch) free(xsearch);
  if (rst) free(rst);
  if (interp2donor) free(interp2donor);
  if (cancelList) deallocateLinkList2(cancelList);
  if (ctag) free(ctag);
  if (pointsPerCell) free(pointsPerCell);
  if (rxyz) free(rxyz);
  if (picked) free(picked);
  if (rxyzCart) free(rxyzCart);
  if (donorIdCart) free(donorIdCart);
  if (pickedCart) free(pickedCart);
  if (ctag_cart) free(ctag_cart);
  // need to add code here for other objects as and
  // when they become part of MeshBlock object  
};
//
// set user specified node and cell resolutions
//
void MeshBlock::setResolutions(double *nres,double *cres)
{
  userSpecifiedNodeRes=nres;
  userSpecifiedCellRes=cres;
}
