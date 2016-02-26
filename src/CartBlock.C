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
#include "CartBlock.h"
#include "CartGrid.h"
extern "C" {
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);
  void deallocateLinkList3(INTEGERLIST2 *temp);
  void deallocateLinkList4(INTERPLIST2 *temp);
  void insertInList(DONORLIST **donorList,DONORLIST *temp1);
  void get_amr_index_xyz( int nq,int i,int j,int k,
			  int pBasis,
			  int nX,int nY,int nZ,
			  int nf,
			  double xlo[3],double dx[3],
			  double qnodes[],
			  int* index, double* xyz);
    void amr_index_to_ijklmn(int pBasis,int nX,int nY,int nZ, int nf, int nq,
			     int index, int* ijklmn);
    int checkHoleMap(double *x,int *nx,int *sam,double *extents);
}
void CartBlock::getInterpolatedData(int *nints,int *nreals,int **intData,
				    double **realData,
				    double *q,
				    int nvar, int interptype)
{
}


void CartBlock::update(double *qval, int index,int nq)
{
  int i;
  for(i=0;i<nq;i++)
    q[index+d3*i]=qval[i];
}

  
void CartBlock::preprocess(CartGrid *cg)
  {
    int nf;    
    for(int n=0;n<3;n++) xlo[n]=cg->xlo[3*global_id+n];
    for(int n=0;n<3;n++) dx[n]=cg->dx[3*global_id+n];
    dims[0]=cg->ihi[3*global_id]  -cg->ilo[3*global_id  ]+1;
    dims[1]=cg->ihi[3*global_id+1]-cg->ilo[3*global_id+1]+1;
    dims[2]=cg->ihi[3*global_id+2]-cg->ilo[3*global_id+2]+1;
    pdegree=cg->porder[global_id];
    p3=(pdegree+1)*(pdegree+1)*(pdegree+1);
    nf=cg->nf;
    qstride=cg->qstride;
    donor_frac=cg->donor_frac;
    qnode=cg->qnode;
    d1=dims[0];
    d2=dims[0]*dims[1];
    d3=d2*dims[2];
  };

void CartBlock::initializeLists(void)
{
  int i;
  for(i=0;i<ndof;i++) deallocateLinkList(donorList[i]);
  deallocateLinkList4(interpList);
  interpList=NULL;
}


void CartBlock::insertInInterpList(int procid,int remoteid,double *xtmp)
{
  int i,n;
  int ix[3];
  double *rst;
  rst=(double *)malloc(sizeof(double)*3);
  if (interpList==NULL) 
    {
      interpList=(INTERPLIST2 *)malloc(sizeof(INTERPLIST));
      listptr=interpList;
    }
  else
    {
      listptr->next=(INTERPLIST2 *)malloc(sizeof(INTERPLIST));
      listptr=listptr->next;
    }

  listptr->next=NULL;
  listptr->inode=NULL;
  listptr->weights=NULL;
  listptr->receptorInfo[0]=procid;
  listptr->receptorInfo[1]=remoteid;

  for(n=0;n<3;n++)
    {
      ix[n]=(xtmp[n]-xlo[n])/dims[n];
      rst[n]=(xtmp[n]-xlo[n]+ix[n]*dx[n])/dx[n];
    }

  listptr->nweights=(pdegree+1)*(pdegree+1)*(pdegree+1);
  listptr->inode=(int *)malloc(sizeof(int)*3);
  listptr->inode[0]=ix[1];
  listptr->inode[1]=ix[2];
  listptr->inode[2]=ix[3];
  donor_frac(&pdegree,rst,&(listptr->nweights),(listptr->weights));
  
  free(rst);
}
  
void CartBlock::insertInDonorList(int senderid,int index,int meshtagdonor,int remoteid,double cellRes)
{
  DONORLIST *temp1;
  int ijklmn[6];
  int pointid;

  amr_index_to_ijklmn(pdegree,dims[0],dims[1],dims[2],nf,qstride,index,ijklmn);
  pointid=ijklmn[5]*(pdegree+1)*(pdegree+1)*d3+
          ijklmn[4]*(pdegree+1)*d3+
          ijklmn[3]*d3+
          ijklmn[2]*d2+
          ijklmn[1]*d1+
          ijklmn[0];
  temp1->donorData[0]=senderid;
  temp1->donorData[1]=meshtagdonor;
  temp1->donorData[2]=remoteid;
  temp1->donorRes=cellRes;
  insertInList(&donorList[pointid],temp1);
}

void CartBlock::processDonors(HOLEMAP *holemap, int nmesh)
{
  int i,j,k,l,m,n,h,p;
  int ibcount,idof,meshtagdonor,icount;
  DONORLIST *temp;
  int *iflag;
  double *xtmp;
  int *index;
  //
  // first mark hole points
  //
  iflag=(int *)malloc(sizeof(int)*nmesh);
  index=(int *)malloc(sizeof(int)*p3);
  xtmp=(double *)malloc(sizeof(int)*p3);
  //
  ibcount=-1;
  idof=-1;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	{
	  ibcount++;
	  ibl[ibcount]=1;
	  get_amr_index_xyz(qstride,i,j,k,p,dims[0],dims[1],dims[2],nf,
			    xlo,dx,qnode,index,xtmp);
	  for(p=0;p<p3;p++)
	    {
	      idof++;
	      if (donorList[idof]==NULL)
		{
		  for(h=0;h<nmesh;h++)
		    if (holemap[h].existWall)
		      {
			if (checkHoleMap(&xtmp[3*p],holemap[h].nx,holemap[h].sam,holemap[h].extents))
			  {
			    ibl[ibcount]=0;
			    break;
			  }
		      }
		}
	      else
		{
		  temp=donorList[i];
		  for(h=0;h<nmesh;h++) iflag[h]=0;
		  while(temp!=NULL) 
		    {
		      meshtagdonor=temp->donorData[1]-BASE;
		      iflag[meshtagdonor]=1;
		      temp=temp->next;
		    }
		  for(h=0;h<nmesh;h++)
		    {
		      if (holemap[h].existWall)
			{
			  if (!iflag[h])
			    if (checkHoleMap(&xtmp[3*p],holemap[h].nx,holemap[h].sam,holemap[h].extents))
			      {
				ibl[ibcount]=0;
				break;
			      }
			}
		    }
		}
	    }
	}

  ibcount=-1;
  idof=-1;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	{
	  ibcount++;
	  if (ibl[ibcount]==0) 
	    {
	      for(p=0;p<p3;p++)
		{
		  idof++;
		  if (donorList[idof]!=NULL)
		    {
		      temp=donorList[idof];
		      while(temp!=NULL)
			{
			  temp->cancel=1;
			  temp=temp->next;
			}
		    }
		}

	    }
	  else
	    {
	      icount=0;
	      for(p=0;p<p3;p++)
		{
		  idof++;
		  if (donorList[idof]!=NULL)
		    {
		      temp=donorList[idof];
		      while(temp!=NULL)
			{
			  if (temp->donorRes < BIGVALUE)
			    {
			      icount++;
			      break;
			    }
			  temp=temp->next;
			}
		    }
		}
	      if (icount==p3) 
		{
		  ibl[ibcount]=-1;
		}
	      else
		{
		  for(p=0;p<p3;p++)
		    {
		      idof++;
		      if (donorList[idof]!=NULL)
			{
			  temp=donorList[idof];
			  while(temp!=NULL)
			    {
			      temp->cancel=1;
			      temp=temp->next;
			    }
			}
		    } 
		}
	    }
	}
}
			      

void CartBlock::getCancellationData(int *cancelledData, int *ncancel)
{
  int i,j,k,p,m;
  int idof;
  DONORLIST *temp;
  idof=-1;
  m=0;
  *ncancel=0;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	for(p=0;p<p3;p++)
	  {
	    idof++;
	    if (donorList[idof]!=NULL)
	      {
		temp=donorList[idof];
		while(temp!=NULL)
		  {
		    if (temp->cancel==1) {
		      *ncancel++;
		      cancelledData[m++]=temp->donorData[0];
		      cancelledData[m++]=1;
		      cancelledData[m++]=temp->donorData[2];
		    }
		  }
	      }
	  }
}

