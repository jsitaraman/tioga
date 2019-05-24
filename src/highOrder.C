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
#include "codetypes.h"
#include "MeshBlock.h"

#define ROW 0
#define COLUMN 1
#define NFRAC 1331

extern "C" 
{
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
  int checkHoleMap(double *x,int *nx,int *sam,double *extents);
}

void MeshBlock::getCellIblanks2(void)
{
  int i;
  int n,nvert,m;
  int icell;
  int inode[8];
  int ncount,flag;
  int verbose;

  icell=0;
  if (iblank_cell==NULL) iblank_cell=(int *)malloc(sizeof(int)*ncells);
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
        {
          flag=1;
          iblank_cell[icell]=1;
          verbose=0;
          //if (myid==763 && icell==7975) {
          // verbose=1;
          //}
          ncount=0;
          for(m=0;m<nvert && flag;m++)
            {
              inode[m]=vconn[n][nvert*i+m]-BASE;
              if (verbose) {
                TRACEI(m);
                TRACEI(inode[m]);
                TRACEI(iblank[inode[m]]);
              }
              if (iblank[inode[m]]==0)
                {
                  iblank_cell[icell]=0;
                  flag=0;
                }
              ncount=ncount+(iblank[inode[m]]==-1);
            }
          if (verbose) {
            TRACEI(icell);
            TRACEI(ncount);
            TRACEI(nvert);
	  }
          if (flag)
            {
              if (ncount ==nvert) iblank_cell[icell]=-1;
	      //            if (ncount > 0)  iblank_cell[icell]=-1;
	      //            if (ncount >= nvert/2) iblank_cell[icell]=-1;
            }
          icell++;
        }
    }
}


void MeshBlock::getCellIblanks(void)
{
  int i;
  int n,nvert,m;
  int icell;
  int inode[8];
  int ncount,flag;
  int verbose;
  int *ibl;
  
  if (iblank_reduced) 
    {
      ibl=iblank_reduced;
    }
  else 
    {
      ibl=iblank;
    }

  icell=0;
  if (iblank_cell==NULL) iblank_cell=(int *)malloc(sizeof(int)*ncells);
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  flag=1;
	  iblank_cell[icell]=1;
          verbose=0;
	  //if (myid==763 && icell==7975) {
	  // verbose=1;
          //}
	  ncount=0;
	  for(m=0;m<nvert && flag;m++)
	    {
	      inode[m]=vconn[n][nvert*i+m]-BASE;
	      if (verbose) {
                TRACEI(m);
                TRACEI(inode[m]);
                TRACEI(ibl[inode[m]]);
              }
	      if (ibl[inode[m]]==0)
		{
		  iblank_cell[icell]=0;
		  flag=0;
		}
	      ncount=ncount+(ibl[inode[m]]== -1);
	    }
	  if (verbose) {
	    TRACEI(icell);
            TRACEI(ncount);
            TRACEI(nvert);
           } 
	  if (flag) 
	    {
	      if (ncount ==nvert) iblank_cell[icell]=-1;
//	      if (ncount > 0)  iblank_cell[icell]=-1;
//	      if (ncount >= nvert/2) iblank_cell[icell]=-1;
	    }
	  icell++;
	}
    }
}

void MeshBlock::clearOrphans(HOLEMAP *holemap, int nmesh,int *itmp)
{
  int i,j,k,m;
  int reject;

  if (ihigh) 
    {
      m=0;
      for(i=0;i<nreceptorCells;i++)
	{
	  reject=0;
	  for(j=0;j<pointsPerCell[i];j++)
	    {
	      if (itmp[m]==0) 
		{
		  reject=2;
		  for(k=0;k<nmesh;k++)
		    if (k!=(meshtag-BASE) && holemap[k].existWall) 
		      {
			if (checkHoleMap(&rxyz[3*m],holemap[k].nx,holemap[k].sam,holemap[k].extents)) 
			  {
			    reject=1;
			    break;
			  }
		      }
		}
	      m++;
	    }
	  if (reject==1) 
	    {
	      iblank_cell[ctag[i]-1]=0; // changed to hole if inside hole map
	    }
	  else if (reject==2) 
	    {
	      iblank_cell[ctag[i]-1]=1; // changed to field if not inside hole map
	    }
	}
    }
  else
    {
      m=0;
      for(i=0;i<nnodes;i++)
	{
	  if (picked[i]) 
	    {
	      reject=0;
	      if (itmp[m]==0) reject=1;
	      if (reject) iblank[i]=1; // changed to field for near-body
                                       // perhaps not the right thing to do
              m++;
	    }
	}
    }
}


void MeshBlock::getInternalNodes(void)
{
  int i,m,j,icell,indx,isum,i3,n,nvert;
  //
  nreceptorCells=0;
  //
  if (ctag!=NULL) TIOGA_FREE(ctag);
  ctag=(int *)malloc(sizeof(int)*ncells);
  //
  for(i=0;i<ncells;i++)
      if (iblank_cell[i]==-1) ctag[nreceptorCells++]=i+BASE;
  //
  if (ihigh) 
    {
      if (pointsPerCell!=NULL) TIOGA_FREE(pointsPerCell);
      pointsPerCell=(int *)malloc(sizeof(int)*nreceptorCells);
      //
      maxPointsPerCell=0;
      ntotalPoints=0;
      //
      for(i=0;i<nreceptorCells;i++)
	{
	  get_nodes_per_cell(&(ctag[i]),&(pointsPerCell[i]));
	  ntotalPoints+=pointsPerCell[i];
	  maxPointsPerCell=TIOGA_MAX(maxPointsPerCell,pointsPerCell[i]);
      }
      //
      if (rxyz !=NULL) TIOGA_FREE(rxyz);
      //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
      rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);
      //
      m=0;
      for(i=0;i<nreceptorCells;i++)
	{
	  get_receptor_nodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
	  m+=(3*pointsPerCell[i]);
	}
    }
  else
    {
      ntotalPoints=0;      
      if (picked !=NULL) TIOGA_FREE(picked);
      picked=(int *) malloc(sizeof(int)*nnodes);
      for(i=0;i<nnodes;i++) {
        picked[i]=0;
        if (iblank[i]==-1) {
          picked[i]=1;
          ntotalPoints++;
         }
      }
  
      if (rxyz !=NULL) TIOGA_FREE(rxyz);
      rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);
      m=0;
      for(i=0;i<nnodes;i++)
	if (picked[i]) 
	  {
	    i3=3*i;
	    for(j=0;j<3;j++)
	      rxyz[m++]=x[i3+j];
	  }
    }
}

void MeshBlock::getExtraQueryPoints(OBB *obc,
			       int *nints,int **intData,
			       int *nreals, double **realData)
{
  int i,j,k;
  int i3;
  double xd[3];
  int *inode;
  int iptr;
  int m;

  inode=(int *)malloc(sizeof(int)*ntotalPoints);
  *nints=*nreals=0; 
  for(i=0;i<ntotalPoints;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(rxyz[i3+k]-obc->xc[k])*obc->vec[j][k];

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
      (*realData)[m++]=rxyz[i3];
      (*realData)[m++]=rxyz[i3+1];
      (*realData)[m++]=rxyz[i3+2];
      (*realData)[m++]=BIGVALUE;
    }
  //
  TIOGA_FREE(inode);
}  

void MeshBlock::processPointDonors(void)
{
  int i,j,m,n;
  int isum,nvert,i3,ivert;
  double *frac;
  int icell;
  int ndim;
  int ierr;
  double xv[8][3];
  double xp[3];
  double frac2[8];
  //
  ndim=NFRAC;
  frac=(double *) malloc(sizeof(double)*ndim);
  interp2ListSize = ninterp2;
  ninterp2=0;
  //
  for(i=0;i<nsearch;i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]]==1) ninterp2++;
  //
  if (interpList2) {
    for(i=0;i<interp2ListSize;i++)
      {
	if (interpList2[i].inode) TIOGA_FREE(interpList2[i].inode);
	if (interpList2[i].weights) TIOGA_FREE(interpList2[i].weights);
      }
    TIOGA_FREE(interpList2);
  } 
  interp2ListSize = ninterp2;
  interpList2=(INTERPLIST *)malloc(sizeof(INTERPLIST)*interp2ListSize);
  for(i=0;i<interp2ListSize;i++)
    {
      interpList2[i].inode=NULL;
      interpList2[i].weights=NULL;
   }
  //  
  //printf("nsearch=%d %d\n",nsearch,myid);
  m=0;
  for(i=0;i<nsearch;i++)
    {
      if (donorId[i] > -1 && iblank_cell[donorId[i]]==1) 
	{
	  if (ihigh) 
	    {
	      icell=donorId[i]+BASE;
	      interpList2[m].inode=(int *) malloc(sizeof(int));		  
	      interpList2[m].nweights=0;
	      donor_frac(&(icell),
			 &(xsearch[3*i]),
			 &(interpList2[m].nweights),
			 &(interpList2[m].inode[0]),
			 frac,
			 &(rst[3*i]),&ndim);
	      interpList2[m].weights=(double *)malloc(sizeof(double)*interpList2[m].nweights);
	      for(j=0;j<interpList2[m].nweights;j++)
		interpList2[m].weights[j]=frac[j];
	      interpList2[m].receptorInfo[0]=isearch[3*i];
	      interpList2[m].receptorInfo[1]=isearch[3*i+1];
	      interpList2[m].receptorInfo[2]=isearch[3*i+2];
	      m++;
	    }
	  else
	    {
	      icell=donorId[i];
	      isum=0;
	      for(n=0;n<ntypes;n++)
		{
		  isum+=nc[n];
		  if (icell < isum) 
		    {
		      icell=icell-(isum-nc[n]);
		      break;
		    }
		}
	      nvert=nv[n];
	      interpList2[m].inode=(int *) malloc(sizeof(int)*nvert);
	      interpList2[m].nweights=nvert;
	      interpList2[m].weights=(double *)malloc(sizeof(double)*interpList2[m].nweights);
	      for(ivert=0;ivert<nvert;ivert++)
		{
		  interpList2[m].inode[ivert]=vconn[n][nvert*icell+ivert]-BASE;
		  i3=3*interpList2[m].inode[ivert];
		  for(j=0;j<3;j++) xv[ivert][j]=x[i3+j];
		}
	      xp[0]=xsearch[3*i];
	      xp[1]=xsearch[3*i+1];
	      xp[2]=xsearch[3*i+2];
	      computeNodalWeights(xv,xp,frac2,nvert);
	      for(j=0;j<nvert;j++)
		interpList2[m].weights[j]=frac2[j];
	      interpList2[m].receptorInfo[0]=isearch[3*i];
	      interpList2[m].receptorInfo[1]=isearch[3*i+1];
	      interpList2[m].receptorInfo[2]=isearch[3*i+2];
	      m++;
	    }
	}
    }
  TIOGA_FREE(frac);
}

void MeshBlock::getInterpolatedSolutionAtPoints(int *nints,int *nreals,int **intData,
						double **realData,
						double *q,
						int nvar, int interptype)
{
  int i;
  int k,m,inode;
  double weight;
  double *qq;
  int icount,dcount;
  //
  qq=(double *)malloc(sizeof(double)*nvar);
  //
  (*nints)=ninterp2;
  (*nreals)=ninterp2*nvar;
  if ((*nints)==0) {
	TIOGA_FREE(qq);
	return;
  }
  //
  (*intData)=(int *)malloc(sizeof(int)*3*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  icount=dcount=0;
  //
  if (ihigh) 
    {
      if (interptype==ROW)
	{    
	  for(i=0;i<ninterp2;i++)
	    {
	      for(k=0;k<nvar;k++) qq[k]=0;
	      inode=interpList2[i].inode[0]-BASE;
	      for(m=0;m<interpList2[i].nweights;m++)
		{
		  weight=interpList2[i].weights[m];
		  //if (weight < 0 || weight > 1.0) {
		  //	TRACED(weight);
		  //	printf("warning: weights are not convex\n");
		  //    }
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[inode+m*nvar+k]*weight;
		}
	      (*intData)[icount++]=interpList2[i].receptorInfo[0];
	      (*intData)[icount++]=interpList2[i].receptorInfo[1];
	      (*intData)[icount++]=interpList2[i].receptorInfo[2];
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
    }
  else
    {
      if (interptype==ROW)
	{    
	  for(i=0;i<ninterp2;i++)
	    {
	      for(k=0;k<nvar;k++) qq[k]=0;
	      for(m=0;m<interpList2[i].nweights;m++)
		{
		  inode=interpList2[i].inode[m];
		  weight=interpList2[i].weights[m];
		  if (weight < -TOL || weight > 1.0+TOL) {
                    TRACED(weight);
                    printf("warning: weights are not convex 2\n");
                   }
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[inode*nvar+k]*weight;
		}
	      (*intData)[icount++]=interpList2[i].receptorInfo[0];
	      (*intData)[icount++]=interpList2[i].receptorInfo[1];
	      (*intData)[icount++]=interpList2[i].receptorInfo[2];
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
      else if (interptype==COLUMN)
	{
	  for(i=0;i<ninterp2;i++)
	    {
	      for(k=0;k<nvar;k++) qq[k]=0;
	      for(m=0;m<interpList2[i].nweights;m++)
		{
		  inode=interpList2[i].inode[m];
		  weight=interpList2[i].weights[m];
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[k*nnodes+inode]*weight;
		}
	      (*intData)[icount++]=interpList2[i].receptorInfo[0];
	      (*intData)[icount++]=interpList2[i].receptorInfo[1];
	      (*intData)[icount++]=interpList2[i].receptorInfo[2];
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
    }
  //
  // no column-wise storage for high-order data
  //
  TIOGA_FREE(qq);
}
	
void MeshBlock::updatePointData(double *q,double *qtmp,int nvar,int interptype)  
{
  int i,j,k,n,m;
  double *qout;
  int index_out;
  int npts;
  //
  if (ihigh) 
    {
      npts=NFRAC;
      qout=(double *)malloc(sizeof(double)*nvar*npts);	
      //
      m=0;
      for(i=0;i<nreceptorCells;i++)
	{
	  if (iblank_cell[ctag[i]-1]==-1) 
	    {
	      convert_to_modal(&(ctag[i]),&(pointsPerCell[i]),&(qtmp[m]),&npts,
			       &index_out,qout);
	      index_out-=BASE;
	      k=0;
	      for(j=0;j<npts;j++)
		for(n=0;n<nvar;n++)
		  {
		    q[index_out+j*nvar+n]=qout[k];
		    k++;
		  }
	    }
	  m+=(pointsPerCell[i]*nvar);
	}
      TIOGA_FREE(qout);
    }
  else
    {
      m=0;
      for(i=0;i<nnodes;i++)
	{
	  if (picked[i]) 
	    {
	      if (iblank[i]==-1) updateSolnData(i,&(qtmp[m]),q,nvar,interptype);
	      m+=nvar;
	    }
	}
    }
}
