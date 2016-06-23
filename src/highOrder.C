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

#define ROW 0
#define COLUMN 1
#define NFRAC 1331

#define FRINGE -1
#define HOLE 0
#define NORMAL 1

extern "C" 
{
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}

void MeshBlock::getCellIblanks(void)
{
  std::vector<int> inode;

  int icell = 0;
  if (iblank_cell==NULL) iblank_cell=(int *)malloc(sizeof(int)*ncells);

  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    inode.resize(nvert);
    for (int i = 0; i < nc[n]; i++)
    {
      int flag=1;
      iblank_cell[icell]=1;
      int ncount=0;
      for (int m = 0; m < nvert && flag; m++)
      {
        inode[m] = vconn[n][nvert*i+m]-BASE;
        if (iblank[inode[m]] == 0)
        {
          iblank_cell[icell] = 0;
          flag = 0;
        }
        ncount = ncount + (iblank[inode[m]]==-1);
      }

      if (flag)
      {
        if (ncount ==nvert) iblank_cell[icell] = -1;
        // if (ncount >= nvert/2) iblank_cell[icell]=-1;
      }
      icell++;
    }
  }
}

void MeshBlock::setArtificialBoundaries(void)
{
  std::vector<int> iblankCell(ncells);

  if (iblank_cell == NULL) iblank_cell = new int(ncells);

  /* Initialized iblank_cell to all normal cells and
   * blank all cells containing a hole node */
  unsigned int ic = 0;
  for (int n = 0; n < ntypes; n++)
  {
    for (int i = 0; i < nc[n]; i++)
    {
      unsigned int nFringe = 0;
      iblank_cell[ic] = NORMAL;
      iblankCell[ic] = NORMAL;
      for (int j = 0; j < nv[n] && flag; j++)
      {
        int iv = vconn[n][nvert*i+j]-BASE;
        if (iblank[iv] == HOLE)
        {
          iblank_cell[ic] = HOLE;
          iblankCell[ic] = HOLE;
          break;
        }
        else if (iblank[iv] == FRINGE)
        {
          nFringe++;
        }
      }

      if (nFringe == nv[n])
      {
        //iblank_cell[ic] = FRINGE;
        iblank_cell[ic] = HOLE; /// TODO
        iblankCell[ic] = FRINGE;
      }

      ic++;
    }
  }

  /// TODO: implement more sophisticated algorithm (i.e. Galbraith's method)
//  /* ---- If nDims == 3 ---- */
//  ic = 0;
//  for (int n = 0; n < ntypes; n++)
//  {
//    for (int i = 0; i < nc[n]; i++)
//    {
//      if (iblank_cell[ic] != NORMAL)
//      {
//        for (int j = 0; j < nf[n]; j++) {
//          int ff = c2f(ic,j);
//          int ic2 = (f2c(ff,0) != ic) ? f2c(ff,0) : f2c(ff,1);
//          if (ic2 > -1) {
//            if (iblankCell[ic2] == NORMAL) {
//              iblankCell[ic2] = FRINGE;
//            }
//          } else {
//            // MPI Boundary
//            int F = findFirst(mpiFaces,c2f(ic,j));
//            if (F > -1) {
//              mpiFringeFaces.push_back(mpiFaces[F]);
//            }
//          }
//        }
//      }
//    }
//  }
//  /* ---- End if nDims == 3 ---- */
}

void MeshBlock::clearOrphans(int *itmp)
{
  int i,j,m;
  int reject;

  if (ihigh) 
    {
      m=0;
      for(i=0;i<nreceptorCells;i++)
	{
	  reject=0;
	  for(j=0;j<pointsPerCell[i];j++)
	    {
	      if (itmp[m]==0) reject=1;
	      m++;
	    }
	  if (reject) 
	    {
	      iblank_cell[ctag[i]-1]=0; // changed to hole for off-body
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
  nreceptorCells=0;

  if (ctag!=NULL) free(ctag);
  ctag = (int *)malloc(sizeof(int)*ncells);

  /* Gather a list of cell IDs for all receptor (fringe) cells */
  for (int i = 0; i < ncells; i++)
    if (iblank_cell[i] == -1)
      ctag[nreceptorCells++] = i+1;

  if (ihigh)
  {
    /* Get a list of positions of all internal DOFs (solution points)
     * for all fringe cells */

    if (pointsPerCell != NULL) free(pointsPerCell);
    pointsPerCell = (int *)malloc(sizeof(int)*nreceptorCells);

    /* Get total number of internal nodes (solution points) for our fringe cells */
    maxPointsPerCell=0;
    ntotalPoints=0;

    for (int i = 0; i < nreceptorCells; i++)
    {
      get_nodes_per_cell(&(ctag[i]),&(pointsPerCell[i]));
      ntotalPoints += pointsPerCell[i];
      maxPointsPerCell = max(maxPointsPerCell,pointsPerCell[i]);
    }

    if (rxyz != NULL) free(rxyz);
    //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
    rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);

    int m = 0;
    for (int i = 0; i < nreceptorCells; i++)
    {
      get_receptor_nodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
      m += (3*pointsPerCell[i]);
    }
  }
  else
  {
    /* Gather a list of positions of all mesh points tagged as 'fringe'
     * Will this ever be executed? */

    ntotalPoints=0;
    if (picked !=NULL) free(picked);
    picked = (int *) malloc(sizeof(int)*nnodes);

    for (int i = 0; i < nnodes; i++) {
      picked[i] = 0;
      if (iblank[i] == -1) {
        picked[i] = 1;
        ntotalPoints++;
      }
    }

    if (rxyz !=NULL) free(rxyz);
    rxyz = (double *)malloc(sizeof(double)*ntotalPoints*3);

    int m = 0;
    for (int i = 0; i < nnodes; i++)
    {
      if (picked[i])
      {
        int i3 = 3*i;
        for (int j = 0; j < 3; j++)
          rxyz[m++] = x[i3+j];
      }
    }
  }
}

void MeshBlock::getBoundaryNodes(void)
{
  /// TODO
  /** Use cell/face connectivity + callback funcs
   *  to get locations of AB flux points */

  nreceptorFaces=0;

  if (ftag!=NULL) free(ctag);
  ftag = (int *)malloc(sizeof(int)*nfaces);

  /* Gather a list of cell IDs for all receptor (fringe) cells */
  for (int i = 0; i < nfaces; i++)
    if (iblank_face[i] == -1)
      ftag[nreceptorFaces++] = i+1;

  /* Get a list of positions of all internal DOFs (solution points)
   * for all fringe cells */

  if (pointsPerFace != NULL) free(pointsPerFace);
  pointsPerFace = (int *)malloc(sizeof(int)*nreceptorCells);

  /* Get total number of internal nodes (solution points) for our fringe cells */
  maxPointsPerFace=0;
  ntotalPoints=0;

  for (int i = 0; i < nreceptorCells; i++)
  {
    get_nodes_per_face(&(ftag[i]),&(pointsPerFace[i])); /// Access by face? or by cell+local face ID?
    ntotalPoints += pointsPerFace[i];
    maxPointsPerFace = max(maxPointsPerFace,pointsPerFace[i]);
  }

  if (rxyz != NULL) free(rxyz);
  //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
  rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);

  int m = 0;
  for (int i = 0; i < nreceptorCells; i++)
  {
    /// Access by face ID, or cell ID + local-face ID?
    get_face_nodes(&(ftag[i]),&(pointsPerFace[i]),&(rxyz[m]));
    m += (3*pointsPerFace[i]);
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
      (*realData)[m++]=rxyz[i3];
      (*realData)[m++]=rxyz[i3+1];
      (*realData)[m++]=rxyz[i3+2];
    }
  //
  free(inode);
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
	if (interpList2[i].inode) free(interpList2[i].inode);
	if (interpList2[i].weights) free(interpList2[i].weights);
      }
    free(interpList2);
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

	      /* Use High-Order callback function to get interpolation weights */
	      donor_frac(&(icell),
			 &(xsearch[3*i]),
			 &(interpList2[m].nweights),
			 &(interpList2[m].inode[0]),
			 frac,
			 &(rst[3*i]),&ndim);

	      interpList2[m].weights=(double *)malloc(sizeof(double)*interpList2[m].nweights);
	      for(j=0;j<interpList2[m].nweights;j++)
		interpList2[m].weights[j]=frac[j];
	      interpList2[m].receptorInfo[0]=isearch[2*i];
	      interpList2[m].receptorInfo[1]=isearch[2*i+1];
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
	      interpList2[m].receptorInfo[0]=isearch[2*i];
	      interpList2[m].receptorInfo[1]=isearch[2*i+1];
	      m++;
	    }
	}
    }
  free(frac);
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
	free(qq);
	return;
  }
  //
  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
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
		  //	traced(weight);
		  //	printf("warning: weights are not convex\n");
		  //    }
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[inode+m*nvar+k]*weight;
		}
	      (*intData)[icount++]=interpList2[i].receptorInfo[0];
	      (*intData)[icount++]=interpList2[i].receptorInfo[1];
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
		  if (weight < 0 || weight > 1.0) {
                    traced(weight);
                    printf("warning: weights are not convex\n");
                   }
		  for(k=0;k<nvar;k++)
		    qq[k]+=q[inode*nvar+k]*weight;
		}
	      (*intData)[icount++]=interpList2[i].receptorInfo[0];
	      (*intData)[icount++]=interpList2[i].receptorInfo[1];
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
	      for(k=0;k<nvar;k++)
		(*realData)[dcount++]=qq[k];
	    }
	}
    }
  //
  // no column-wise storage for high-order data
  //
  free(qq);
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
      free(qout);
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
