// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"
#include "MeshBlock.h"
#include <cstring>
#include <stdexcept>

extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  double computeCellVolume(double xv[8][3],int nvert);
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);  
  double tdot_product(double a[3],double b[3],double c[3]);
  void getobbcoords(double xc[3],double dxc[3],double vec[3][3],double xv[8][3]);
  void transform2OBB(double xv[3],double xc[3],double vec[3][3],double xd[3]);
  void writebbox(OBB *obb,int bid);
  void writebboxdiv(OBB *obb,int bid);
}

void MeshBlock::setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
			int *wbcnodei,int *obcnodei,
                        int ntypesi,int *nvi,int *nci,int **vconni,
                        uint64_t* cell_gid, uint64_t* node_gid)
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
  nodeGID = node_gid;
  //
  //TRACEI(nnodes);
  //for(i=0;i<ntypes;i++) TRACEI(nc[i]);
  ncells=0;
  for(i=0;i<ntypes;i++) ncells+=nc[i];

#ifdef TIOGA_HAS_NODEGID
  if (nodeGID == NULL)
      throw std::runtime_error("#tioga: global IDs for nodes not provided");
#endif
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
  if (check_uniform_hex_flag) {
      check_for_uniform_hex();
      if (uniform_hex) create_hex_cell_map();
  }
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
  
  for(int j=0;j<3;j++)
    {
     mapdims[j]=10;
     mapdx[j]=2*obb->dxc[j]/mapdims[j];
    }
  //
  // compute nodal resolution as the average of 
  // all the cells associated with it. This takes care
  // of partition boundaries as well.
  //
  // Also create the inverse map of nodes
  // 
  if (icft) TIOGA_FREE(icft);
  icft=(int *)malloc(sizeof(int)*(mapdims[2]*mapdims[1]*mapdims[0]+1));
  if (invmap) TIOGA_FREE(invmap);
  invmap=(int *)malloc(sizeof(int)*nnodes);
  for(int i=0;i<mapdims[2]*mapdims[1]*mapdims[0]+1;i++) icft[i]=-1;
  icft[0]=0;	  
  int *iptr;
  iptr=(int *)malloc(sizeof(int)*nnodes);
  //
  for(i=0;i<nnodes;i++)
    {
      double xd[3];
      int idx[3];
      if (iflag[i]!=0)  nodeRes[i]/=iflag[i];
      iflag[i]=0;
      iextmp[i]=iextmp1[i]=0;
      
      for(int j=0;j<3;j++)
	{
	  xd[j]=obb->dxc[j];
	  for(int k=0;k<3;k++)
	    xd[j]+=(x[3*i+k]-obb->xc[k])*obb->vec[j][k];
	  idx[j]=xd[j]/mapdx[j];
	}
      int indx=idx[2]*mapdims[1]*mapdims[0]+idx[1]*mapdims[0]+idx[0];
      iptr[i]=icft[indx+1];
      icft[indx+1]=i;
    }
 
 int kc=0;
 for(int i=0;i<mapdims[2]*mapdims[1]*mapdims[0];i++)
  {
   int ip=icft[i+1];
   int m=0;
   while(ip != -1) 
   {
     invmap[kc++]=ip;
     ip=iptr[ip];
     m++;
   }
   icft[i+1]=icft[i]+m;
  }

  TIOGA_FREE(iptr);
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
  // also make the inverse map mask
  if (mapmask) TIOGA_FREE(mapmask);
  mapmask=(int *)malloc(sizeof(int)*mapdims[2]*mapdims[1]*mapdims[0]);
  for(int i=0;i<mapdims[2]*mapdims[1]*mapdims[0];i++) mapmask[i]=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  double xd[3],xc[3],xmin[3],xmax[3];	    
	  int idx[3];
	  itag=0;
	  for(int j=0;j<3;j++) { xmin[j]=BIGVALUE;xmax[j]=-BIGVALUE;}
	  for(m=0;m<nvert;m++)
	    {
	      inode[m]=vconn[n][nvert*i+m]-BASE;
	      if (iflag[inode[m]]) itag=1;
	      for(int j=0;j<3;j++)
		{
		  xd[j]=obb->dxc[j];
		  for(int k=0;k<3;k++) 
		    xd[j]+=(x[3*inode[m]+k]-obb->xc[k])*obb->vec[j][k];
		  xmin[j]=TIOGA_MIN(xd[j],xmin[j]);
		  xmax[j]=TIOGA_MAX(xd[j],xmax[j]);
		} 	    
	    }
          for(int j=0;j<3;j++) { xmin[j]-=TOL; xmax[j]+=TOL;}
	  for(int j=xmin[0]/mapdx[0];j<=xmax[0]/mapdx[0];j++)
            for(int k=xmin[1]/mapdx[1];k<=xmax[1]/mapdx[1];k++)
              for(int l=xmin[2]/mapdx[2];l<=xmax[2]/mapdx[2];l++)
                {
                  idx[0]=TIOGA_MAX(TIOGA_MIN(j,mapdims[0]-1),0);
                  idx[1]=TIOGA_MAX(TIOGA_MIN(k,mapdims[1]-1),0);
                  idx[2]=TIOGA_MAX(TIOGA_MIN(l,mapdims[2]-1),0);
                  mapmask[idx[2]*mapdims[1]*mapdims[0]+idx[1]*mapdims[0]+idx[0]]=1;
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
  char qstr[3];
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
  fprintf(fp, "TITLE =\"Tioga output\"\n");
  fprintf(fp, "VARIABLES = \"X\"\n");
  fprintf(fp, "\"Y\"\n");
  fprintf(fp, "\"Z\"\n");
  fprintf(fp, "\"IBLANK\"\n");
  fprintf(fp, "\"IBLANK_CELL\"\n");
  fprintf(fp, "ZONE T=\"VOL_MIXED\"\n");
  fprintf(fp, " Nodes=%d, Elements=%d, ZONETYPE=FEBrick\n", nnodes, ncells);
  fprintf(fp, " DATAPACKING=BLOCK\n");
  fprintf(fp, " VARLOCATION=([5]=CELLCENTERED)\n");
  fprintf(fp, " DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i+1]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%lf\n",x[3*i+2]);
  for(i=0;i<nnodes;i++) fprintf(fp,"%d.0\n",iblank[i]);
  for(i=0;i<ncells;i++) fprintf(fp,"%d.0\n",iblank_cell[i]);
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
  char qstr[3];
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
  /*
  for(j=0;j<3;j++)
    {
      realData[j]=obb->xc[j];
      realData[j+3]=obb->dxc[j];
    }
  return;
  */
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


void MeshBlock::getReducedOBB2(OBB *obc,double *realData) 
{
  int i,j,k,l,m,n,i3,jmin,kmin,lmin,jmax,kmax,lmax,indx;
  double bbox[6],xd[3];
  double xmin[3],xmax[3],xv[8][3];
  double delta;
  int imin[3],imax[3];

  getobbcoords(obc->xc,obc->dxc,obc->vec,xv);
  for(j=0;j<3;j++) {xmin[j]=BIGVALUE;xmax[j]=-BIGVALUE;};
  for(n=0;n<8;n++)
    {
      transform2OBB(xv[n],obb->xc,obb->vec,xd);
      for(j=0;j<3;j++)
 	{
          xmin[j]=TIOGA_MIN(xmin[j],xd[j]+obb->dxc[j]);
          xmax[j]=TIOGA_MAX(xmax[j],xd[j]+obb->dxc[j]);
        }
    }
  for(j=0;j<3;j++) 
    {
      delta=0.01*(xmax[j]-xmin[j]);
      xmin[j]-=delta;
      xmax[j]+=delta;
      imin[j]=TIOGA_MAX(xmin[j]/mapdx[j],0);
      imax[j]=TIOGA_MIN(xmax[j]/mapdx[j],mapdims[j]-1);
    }
  lmin=mapdims[2]-1;
  kmin=mapdims[1]-1;
  jmin=mapdims[0]-1;
  lmax=kmax=jmax=0;
  for(l=imin[2];l<=imax[2];l++)
    for(k=imin[1];k<=imax[1];k++)
      for(j=imin[0];j<=imax[0];j++)
	{
	  indx=l*mapdims[1]*mapdims[0]+k*mapdims[0]+j;
	  if (mapmask[indx]) {
	    lmin=TIOGA_MIN(lmin,l);
	    kmin=TIOGA_MIN(kmin,k);
	    jmin=TIOGA_MIN(jmin,j);
	    lmax=TIOGA_MAX(lmax,l);
	    kmax=TIOGA_MAX(kmax,k);
	    jmax=TIOGA_MAX(jmax,j);
	  }
	}
  bbox[0]=-obb->dxc[0]+jmin*mapdx[0];
  bbox[1]=-obb->dxc[1]+kmin*mapdx[1];
  bbox[2]=-obb->dxc[2]+lmin*mapdx[2];
  bbox[3]=-obb->dxc[0]+(jmax+1)*mapdx[0];
  bbox[4]=-obb->dxc[1]+(kmax+1)*mapdx[1];
  bbox[5]=-obb->dxc[2]+(lmax+1)*mapdx[2];
  for(j=0;j<3;j++)
    {
      realData[j]=obb->xc[j];
      for(k=0;k<3;k++)
	realData[j]+=((bbox[k]+bbox[k+3])*0.5)*obb->vec[k][j];
      realData[j+3]=(bbox[j+3]-bbox[j])*0.5;
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
  if (myid==0 && meshtag==1) {TRACEI(*nints);} 
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
	      
void MeshBlock::getQueryPoints2(OBB *obc,
			       int *nints,int **intData,
			       int *nreals, double **realData)
{
  int i,j,k,l,il,ik,ij,n,m,i3,iflag;
  int indx,iptr;
  int *inode;
  double delta;
  double xv[8][3],mdx[3],xd[3],xc[3];
  double xmax[3],xmin[3];
  int imin[3],imax[3];
  //
  inode=(int *)malloc(sizeof(int)*nnodes);
  *nints=*nreals=0;
  getobbcoords(obc->xc,obc->dxc,obc->vec,xv);
  for(j=0;j<3;j++) {xmin[j]=BIGVALUE;xmax[j]=-BIGVALUE;};
  //writebbox(obc,1);
  //writebbox(obb,2);
  //writebboxdiv(obb,1);
  //
  for(n=0;n<8;n++)
    {
      transform2OBB(xv[n],obb->xc,obb->vec,xd);
      for(j=0;j<3;j++)
 	{
          xmin[j]=TIOGA_MIN(xmin[j],xd[j]+obb->dxc[j]);
          xmax[j]=TIOGA_MAX(xmax[j],xd[j]+obb->dxc[j]);
        }
    }
  
  for(j=0;j<3;j++) 
    {
      delta=0.01*(xmax[j]-xmin[j]);
      xmin[j]-=delta;
      xmax[j]+=delta;
      imin[j]=TIOGA_MAX(xmin[j]/mapdx[j],0);
      imax[j]=TIOGA_MIN(xmax[j]/mapdx[j],mapdims[j]-1);
      mdx[j]=0.5*mapdx[j];
    }
  //
  // find min/max extends of a single sub-block
  // in OBC axes
  //
  getobbcoords(obc->xc,mdx,obb->vec,xv);
  for(j=0;j<3;j++) {xmin[j]=BIGVALUE;xmax[j]=-BIGVALUE;}
  for(m=0;m<8;m++)
    {
      transform2OBB(xv[m],obc->xc,obc->vec,xd);
      for(j=0;j<3;j++)
	{
	  xmin[j]=TIOGA_MIN(xmin[j],xd[j]);
	  xmax[j]=TIOGA_MAX(xmax[j],xd[j]);
	}
    }
  //printf("xmin :%f %f %f\n",xmin[0],xmin[1],xmin[2]);
  //printf("xmax :%f %f %f\n",xmax[0],xmax[1],xmax[2]);
  //
  // now find the actual number of points
  // that are within OBC using only the 
  // sub-blocks with potential bbox overlap
  //
  //FILE *fp,*fp2,*fp3,*fp4;
  //fp=fopen("v.dat","w");
  //fp2=fopen("v2.dat","w");
  //fp3=fopen("v3.dat","w");
  //fp4=fopen("v4.dat","w");
  for(l=imin[2];l<=imax[2];l++)
    for(k=imin[1];k<=imax[1];k++)
      for(j=imin[0];j<=imax[0];j++)
	{
	  //
	  // centroid of each sub-block 
	  // in OBC axes
	  //
	  xd[0]=-obb->dxc[0]+j*mapdx[0]+mapdx[0]*0.5;
	  xd[1]=-obb->dxc[1]+k*mapdx[1]+mapdx[1]*0.5;
	  xd[2]=-obb->dxc[2]+l*mapdx[2]+mapdx[2]*0.5;
	  for(n=0;n<3;n++)
	    {
	      xc[n]=obb->xc[n];
	      for(ij=0;ij<3;ij++)
		xc[n]+=(xd[ij]*obb->vec[ij][n]);
	    }
          //if (j==0 && k==0 & l==0) {
          //printf("%f %f %f\n",obc->vec[0][0],obc->vec[1][0],obc->vec[2][0]);
          //printf("%f %f %f\n",obc->vec[0][1],obc->vec[1][1],obc->vec[2][1]);
          //printf("%f %f %f\n",obc->vec[0][2],obc->vec[1][2],obc->vec[2][2]);
          //printf("%f %f %f\n",obc->xc[0],obc->xc[1],obc->xc[2]);
          //}
          //fprintf(fp,"%f %f %f\n",xc[0],xc[1],xc[2]);
	  transform2OBB(xc,obc->xc,obc->vec,xd);
          //fprintf(fp2,"%f %f %f\n",xd[0],xd[1],xd[2]);
	  //if (fabs(xd[0]) <= obc->dxc[0] &&
	  //    fabs(xd[1]) <= obc->dxc[1] &&
	  //   fabs(xd[2]) <= obc->dxc[2])
//		{
//                  fprintf(fp3,"%f %f %f\n",xc[0],xc[1],xc[2]);
//		}
          //
	  // check if this sub-block overlaps OBC
	  //
	  iflag=0;
	  for(ij=0;ij<3 && !iflag;ij++) iflag=(iflag || (xmin[ij]+xd[ij] > obc->dxc[ij]));
	  if (iflag) continue;
	  iflag=0;
	  for(ij=0;ij<3 && !iflag;ij++) iflag=(iflag || (xmax[ij]+xd[ij] < -obc->dxc[ij]));
	  if (iflag) continue;
          //fprintf(fp4,"%f %f %f\n",xc[0],xc[1],xc[2]);
	  //
	  // if there overlap
	  // go through points within the sub-block
	  // to figure out what needs to be send
	  //
	  indx=l*mapdims[1]*mapdims[0]+k*mapdims[0]+j;
	  for(m=icft[indx];m<icft[indx+1];m++)
	    {
	      i3=3*invmap[m];
	      for(ik=0;ik<3;ik++) xc[ik]=x[i3+ik];
	      transform2OBB(xc,obc->xc,obc->vec,xd);
	      if (fabs(xd[0]) <= obc->dxc[0] &&
		  fabs(xd[1]) <= obc->dxc[1] &&
		  fabs(xd[2]) <= obc->dxc[2])
		{
		  inode[*nints]=invmap[m];
		  (*nints)++;
		  (*nreals)+=4;
		}
	    }
	}
//  TRACEI(*nints);
//  fclose(fp);
//  fclose(fp2);
//  fclose(fp3);
//  fclose(fp4);
//  int ierr;
//  MPI_Abort(MPI_COMM_WORLD,ierr);
  //
#ifdef TIOGA_HAS_NODEGID
  int nintsPerNode = 3;
#else
  int nintsPerNode = 1;
#endif
  (*intData)=(int *)malloc(sizeof(int)*(*nints) * nintsPerNode);
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  //
  m=0;
  int iidx = 0;
  for(i=0;i<*nints;i++) {
      i3=3*inode[i];
      (*intData)[iidx++]=inode[i];
#ifdef TIOGA_HAS_NODEGID
      std::memcpy(&(*intData)[iidx],&nodeGID[inode[i]], sizeof(uint64_t));
      iidx += 2;
#endif
      (*realData)[m++]=x[i3];
      (*realData)[m++]=x[i3+1];
      (*realData)[m++]=x[i3+2];
      (*realData)[m++]=nodeRes[inode[i]];
  }

  // Adjust nints to the proper array size
  *nints *= nintsPerNode;
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
  if (obh) TIOGA_FREE(obh);
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
  if (receptorIdCart) TIOGA_FREE(receptorIdCart);
  if (icft) TIOGA_FREE(icft);
  if (mapmask) TIOGA_FREE(mapmask);
  if (uindx) TIOGA_FREE(uindx);
  if (invmap) TIOGA_FREE(invmap);
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
  if (ntypes > 1) return;
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
    for(int j=0;j<3;j++)
     for(int k=0;k<3;k++)
       obh->vec[j][k]=0;      
    for(int k=0;k<3;k++)
      obh->vec[0][k]=(xv[1][k]-xv[0][k])/dx[0];
    for(int k=0;k<3;k++)
      obh->vec[1][k]=(xv[3][k]-xv[0][k])/dx[1];
    for(int k=0;k<3;k++)
      obh->vec[2][k]=(xv[4][k]-xv[0][k])/dx[2];
    //obh->vec[0][0]=obh->vec[1][1]=obh->vec[2][2]=1;
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
      xlow[j]=obh->xc[j];
      for (int k=0;k<3;k++)
         xlow[j]-=(obh->dxc[k]*obh->vec[k][j]);
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
