//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is TIOGA_FREE software; you can redistribute it and/or
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
#include "codetypes.h"
#include "MeshBlock.h"
extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  double computeCellVolume(double xv[8][3],int nvert);
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);  
  double tdot_product(double a[3],double b[3],double c[3]);
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
  //TRACEI(nnodes);
  //for(i=0;i<ntypes;i++) TRACEI(nc[i]);
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
  check_for_uniform_hex();
  if (uniform_hex) create_hex_cell_map();
  if (obb) TIOGA_FREE(obb);
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
  int *iextmp,*iextmp1; 
  int iex;
  //
  // do this only once
  // i.e. when the meshblock is first 
  // initialized, cellRes would be NULL in this case
  //

  if(cellRes) TIOGA_FREE(cellRes);
  if(nodeRes) TIOGA_FREE(nodeRes);
    
  //
  cellRes=(double *) malloc(sizeof(double)*ncells);
  nodeRes=(double *) malloc(sizeof(double)*nnodes);
  //
  // this is a local array
  //
  iflag=(int *)malloc(sizeof(int)*nnodes);
  iextmp=(int *) malloc(sizeof(double)*nnodes);
  iextmp1=(int *) malloc(sizeof(double)*nnodes);
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
      iextmp[i]=iextmp1[i]=0;
    }
  //
  // now tag the boundary nodes
  // reuse the iflag array
  //
  //TRACEI(nobc);
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
		  iextmp[inode[m]]=iextmp1[inode[m]]=1;
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
  //
  for(iex=0;iex<mexclude;iex++)
    {
      k=0;
      for(n=0;n<ntypes;n++)
	{
	  nvert=nv[n];
	  for(i=0;i<nc[n];i++)
	    {
	      for(m=0;m<nvert;m++)
		{
		  inode[m]=vconn[n][nvert*i+m]-BASE;
		  if (iextmp[inode[m]]==1) //(iflag[inode[m]]) 
		    {
		      cellRes[k]=BIGVALUE;
		      break;
		    }		    
		}
	      if (cellRes[k]==BIGVALUE) 
		{
		  for(m=0;m<nvert;m++) {
		    inode[m]=vconn[n][nvert*i+m]-BASE;
		    if (iextmp[inode[m]]!=1) iextmp1[inode[m]]=1;
		  }
		}
	      k++;
	    }
	}
      for(i=0;i<nnodes;i++) iextmp[i]=iextmp1[i];	
    }
  TIOGA_FREE(iflag);
  TIOGA_FREE(iextmp);
  TIOGA_FREE(iextmp1);
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
  int *ibl;
  //
  // if fringes were reduced use
  // iblank_reduced
  //
  if (iblank_reduced) 
    {
      ibl=iblank_reduced;
    }
  else
    {
      ibl=iblank;
    }
  //
  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"flow%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\",\"BTAG\" ");
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
	  fprintf(fp,"%lf %lf %lf %d %d ",x[3*i],x[3*i+1],x[3*i+2],ibl[i],meshtag);
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
          fprintf(fp,"%lf %lf %lf %d %d ",x[3*i],x[3*i+1],x[3*i+2],ibl[i],meshtag);
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
  fprintf(fp,"%d\n",nwbc);
  for(i=0;i<nwbc;i++)
    fprintf(fp,"%d\n",wbcnode[i]);
  fprintf(fp,"%d\n",nobc);
  for(i=0;i<nobc;i++)
    fprintf(fp,"%d\n",obcnode[i]);
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
	  wbox[j]=TIOGA_MIN(wbox[j],x[i3+j]);
	  wbox[j+3]=TIOGA_MAX(wbox[j+3],x[i3+j]);
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
		      imin[k]=TIOGA_MIN(imin[k],iv);
		      imax[k]=TIOGA_MAX(imax[k],iv);
		    }
		}
	     for(j=0;j<3;j++)
              {
	       imin[j]=TIOGA_MAX(imin[j],0);
               imax[j]=TIOGA_MIN(imax[j],nx[j]-1);
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
   TIOGA_FREE(iflag);
   TIOGA_FREE(inode);
}

void MeshBlock::getReducedOBB(OBB *obc,double *realData) 
{
  int i,j,k,m,n,i3;
  int nvert;
  bool iflag;
  double bbox[6],xd[3];

  for(j=0;j<3;j++)
    {
      realData[j]=obb->xc[j];
      realData[j+3]=obb->dxc[j];
    }
  return;

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
	      for(j=0;j<3;j++) bbox[j]=TIOGA_MIN(bbox[j],xd[j]);
	      for(j=0;j<3;j++) bbox[j+3]=TIOGA_MAX(bbox[j+3],xd[j]);
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
	      for(j=0;j<3;j++) realData[j]=TIOGA_MIN(realData[j],xd[j]);
	      for(j=0;j<3;j++) realData[j+3]=TIOGA_MAX(realData[j+3],xd[j]);
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
	  (*nreals)+=4;

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
      (*realData)[m++]=nodeRes[inode[i]];
    }
  //
  TIOGA_FREE(inode);
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
  // TIOGA_FREE all data that is owned by this MeshBlock
  // i.e not the pointers of the external code.
  //
  if (cellRes) TIOGA_FREE(cellRes);
  if (nodeRes) TIOGA_FREE(nodeRes);
  if (elementBbox) TIOGA_FREE(elementBbox);
  if (elementList) TIOGA_FREE(elementList);
  if (adt) delete[] adt;
  if (donorList) {
    for(i=0;i<nnodes;i++) deallocateLinkList(donorList[i]);
    TIOGA_FREE(donorList);
  }
  if (interpList) {
    for(i=0;i<interpListSize;i++)
      {
	if (interpList[i].inode) TIOGA_FREE(interpList[i].inode);
	if (interpList[i].weights) TIOGA_FREE(interpList[i].weights);
      }
    TIOGA_FREE(interpList);
  }
  if (interpList2) {
    for(i=0;i<interp2ListSize;i++)
      {
	if (interpList2[i].inode) TIOGA_FREE(interpList2[i].inode);
	if (interpList2[i].weights) TIOGA_FREE(interpList2[i].weights);
      }
    TIOGA_FREE(interpList2);
  }
  if (interpListCart) {
    for(i=0;i<interpListCartSize;i++)
      {
        if (interpListCart[i].inode) TIOGA_FREE(interpListCart[i].inode);
 	if (interpListCart[i].weights) TIOGA_FREE(interpListCart[i].weights);
      }
    TIOGA_FREE(interpListCart);
  }
  // For nalu-wind API the iblank_cell array is managed on the nalu side
  // if (!ihigh) {
  //  if (iblank_cell) TIOGA_FREE(iblank_cell);
  // }
  if (obb) TIOGA_FREE(obb);
  if (isearch) TIOGA_FREE(isearch);
  if (xsearch) TIOGA_FREE(xsearch);
  if (res_search) TIOGA_FREE(res_search);
  if (xtag) TIOGA_FREE(xtag);
  if (rst) TIOGA_FREE(rst);
  if (interp2donor) TIOGA_FREE(interp2donor);
  if (cancelList) deallocateLinkList2(cancelList);
  if (ctag) TIOGA_FREE(ctag);
  if (pointsPerCell) TIOGA_FREE(pointsPerCell);
  if (rxyz) TIOGA_FREE(rxyz);
  if (picked) TIOGA_FREE(picked);
  if (rxyzCart) TIOGA_FREE(rxyzCart);
  if (donorIdCart) TIOGA_FREE(donorIdCart);
  if (pickedCart) TIOGA_FREE(pickedCart);
  if (ctag_cart) TIOGA_FREE(ctag_cart);

  if (tagsearch) TIOGA_FREE(tagsearch);
  if (donorId) TIOGA_FREE(donorId);
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
//
// detect if a given meshblock is a uniform hex
// and create a data structure that will enable faster
// searching
//
void MeshBlock::check_for_uniform_hex(void)
{
  double xv[8][3];
  int hex_present=0;

  for(int n=0;n<ntypes;n++)
    {
      int nvert=nv[n];
      if (nvert==8) {
	hex_present=1;
	for(int i=0;i<nc[n];i++)
	  {
	    int vold=-1;
	    for(int m=0;m<nvert;m++)
	      {
		if (vconn[n][nvert*i+m]==vold) return; // degenerated hex are not uniform
		vold=vconn[n][nvert*i+m];
		int i3=3*(vconn[n][nvert*i+m]-BASE);
		for(int k=0;k<3;k++)
		  xv[m][k]=x[i3+k];
	      }
	    //
	    // check angles to see if sides are 
	    // rectangles
	    //
	    // z=0 side
	    //
	    if (fabs(tdot_product(xv[1],xv[3],xv[0])) > TOL) return;
	    if (fabs(tdot_product(xv[1],xv[3],xv[2])) > TOL) return;
	    //
	    // x=0 side
	    //
	    if (fabs(tdot_product(xv[3],xv[4],xv[0])) > TOL) return;
	    if (fabs(tdot_product(xv[3],xv[4],xv[7])) > TOL) return;
	    //
	    // y=0 side
	    //
	    if (fabs(tdot_product(xv[4],xv[1],xv[0])) > TOL) return;
	    if (fabs(tdot_product(xv[4],xv[1],xv[5])) > TOL) return;
	    //
	    // need to check just one more angle on 
	    // on the corner against 0 (6)
	    //
	    if (fabs(tdot_product(xv[5],xv[7],xv[6])) > TOL) return;
	    //
	    // so this is a hex
	    // check if it has the same size as the previous hex
	    // if not return
	    //
	    if (i==0){
	      dx[0]=tdot_product(xv[1],xv[1],xv[0]);
	      dx[1]=tdot_product(xv[3],xv[3],xv[0]);
	      dx[2]=tdot_product(xv[4],xv[4],xv[0]);
	    }
	    else {
	      if (fabs(dx[0]-tdot_product(xv[1],xv[1],xv[0])) > TOL) return;
	      if (fabs(dx[1]-tdot_product(xv[3],xv[3],xv[0])) > TOL) return;
	      if (fabs(dx[2]-tdot_product(xv[4],xv[4],xv[0])) > TOL) return;
	    }
	  }
      }
    }
  if (hex_present) {
    for(int j=0;j<3;j++) dx[j]=sqrt(dx[j]);
    uniform_hex=1;
    if (obh) TIOGA_FREE(obh);
    obh=(OBB *) malloc(sizeof(OBB));
    for(int k=0;k<3;k++)
      obh->vec[0][k]=(xv[1][k]-xv[0][k])/dx[0];
    for(int k=0;k<3;k++)
      obh->vec[1][k]=(xv[3][k]-xv[0][k])/dx[1];
    for(int k=0;k<3;k++)
      obh->vec[2][k]=(xv[4][k]-xv[0][k])/dx[2];

    for(int j=0;j<3;j++)
     for(int k=0;k<3;k++)
       obh->vec[j][k]=0;      
    obh->vec[0][0]=obh->vec[1][1]=obh->vec[2][2]=1;
    //
    double xd[3];
    double xmax[3];
    double xmin[3];
    xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
    xmin[0]=xmin[1]=xmin[2]=BIGVALUE;	
    //
    for(int i=0;i<nnodes;i++)
    {
      int i3=3*i;
      for(int j=0;j<3;j++) xd[j]=0;
      //
      for(int j=0;j<3;j++)
	for(int k=0;k<3;k++)
	  xd[j]+=x[i3+k]*obh->vec[j][k];
      //
      for(int j=0;j<3;j++)
	{
	  xmax[j]=TIOGA_MAX(xmax[j],xd[j]);
	  xmin[j]=TIOGA_MIN(xmin[j],xd[j]);
	}
    }
        //
    // find the extents of the box
    // and coordinates of the center w.r.t. xc
    // increase extents by 1% for tolerance
    //
    for(int j=0;j<3;j++)
      {
	xmax[j]+=TOL;
	xmin[j]-=TOL;
	obh->dxc[j]=(xmax[j]-xmin[j])*0.5;	
	xd[j]=(xmax[j]+xmin[j])*0.5;
      }
    //
    // find the center of the box in 
    // actual cartesian coordinates
    //
    for(int j=0;j<3;j++)
      {
	obh->xc[j]=0.0;
	for(int k=0;k<3;k++)
	  obh->xc[j]+=(xd[k]*obh->vec[k][j]);
      }    
   }
  return;
}
		
void MeshBlock::create_hex_cell_map(void)
{
  for(int j=0;j<3;j++)
    {
      xlow[j]=obh->xc[j]-obh->dxc[j];
      idims[j]=round(2*obh->dxc[j]/dx[j]);
      dx[j]=(2*obh->dxc[j])/idims[j];
    }
  //
  if (uindx) TIOGA_FREE(uindx);
  uindx=(int *)malloc(sizeof(int)*idims[0]*idims[1]*idims[2]);
  for(int i=0;i<idims[0]*idims[1]*idims[2];uindx[i++]=-1);
  //
  for(int i=0;i<nc[0];i++)
    {
      double xc[3];
      double xd[3];
      int idx[3];
      for(int j=0;j<3;j++)
	{
	  int lnode=vconn[0][8*i]-BASE;
	  int tnode=vconn[0][8*i+6]-BASE;
	  xc[j]=0.5*(x[3*lnode+j]+x[3*tnode+j]);
	}
      for(int j=0;j<3;j++)
	{
	  xd[j]=0;
	  for(int k=0;k<3;k++)
	    xd[j]+=(xc[k]-xlow[k])*obh->vec[j][k];
         idx[j]=xd[j]/dx[j];
	}
       uindx[idx[2]*idims[1]*idims[0]+idx[1]*idims[0]+idx[0]]=i;
    }
}
