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
#include "linCartInterp.h"
#include <assert.h>
#include <stdexcept>
extern "C" {
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);
  void deallocateLinkList3(INTEGERLIST2 *temp);
  void deallocateLinkList4(INTERPLIST2 *temp);
  void insertInList(DONORLIST **donorList,DONORLIST *temp1);
  int checkHoleMap(double *x,int *nx,int *sam,double *extents);
}
void CartBlock::getInterpolatedData(int *nints,int *nreals,int **intData,
				    double **realData,
				    int nvar)
{
  int i,n;
  int index;
  double *qq;
  int *tmpint;
  double *tmpreal;
  int icount,dcount;
  int nintold,nrealold;
  int interpCount=0;
  double weight;
  listptr=interpList;
  while(listptr!=NULL)
    {
      interpCount++;
      listptr=listptr->next;
    }
  if (interpCount > 0) 
    {
      nintold=(*nints);
      nrealold=(*nreals);
      if (nintold > 0) {
       tmpint=(int *)malloc(sizeof(int)*3*(*nints));
       tmpreal=(double *)malloc(sizeof(double)*(*nreals));
       for(i=0;i<(*nints)*3;i++) tmpint[i]=(*intData)[i];
       for(i=0;i<(*nreals);i++) tmpreal[i]=(*realData)[i];
       //
       TIOGA_FREE((*intData));
       TIOGA_FREE((*realData)); // didnt free this before ??
       //  
      }
      (*nints)+=interpCount;
      (*nreals)+=(interpCount*nvar);
      (*intData)=(int *)malloc(sizeof(int)*3*(*nints));
      (*realData)=(double *)malloc(sizeof(double)*(*nreals));
      if (nintold > 0) {
       for(i=0;i<nintold*3;i++) (*intData)[i]=tmpint[i];
       for(i=0;i<nrealold;i++) (*realData)[i]=tmpreal[i];      
       TIOGA_FREE(tmpint);
       TIOGA_FREE(tmpreal);
      }
      listptr=interpList;
      icount=3*nintold;
      dcount=nrealold;
      qq=(double *)malloc(sizeof(double)*nvar);
      while(listptr!=NULL)
      {
        (*intData)[icount++]=listptr->receptorInfo[0];
        (*intData)[icount++]=-1-listptr->receptorInfo[2];
        (*intData)[icount++]=listptr->receptorInfo[1];

        for(n=0;n<nvar;n++) qq[n]=0; // zero out solution

        for(i=0;i<listptr->nweights;i++)
        {
          //Q[nq,nZ+2*nf,nY+2*nf,nX+2*nf]--> C++ Cell storage
          index = (dims[1]+2*nf)*(dims[0]+2*nf)*(listptr->inode[3*i+2]+nf)
              + (dims[0]+2*nf)*(listptr->inode[3*i+1]+nf)
              + (listptr->inode[3*i]+nf);

          for(n=0;n<nvar;n++)
          {
            weight=listptr->weights[i];
            qq[n]+=qcell[index+ncell_nf*n]*weight;
          }
        }

        for(n=0;n<nvar;n++) (*realData)[dcount++]=qq[n]; // update solution

        listptr=listptr->next;
      }
      TIOGA_FREE(qq);
    }
}


void CartBlock::update(double *qval, int index,int nq)
{
  if(index >= ncell_nf){
    int stop = 0;
    stop = stop -1;
    return;
  }

  int i;
  for(i=0;i<nq;i++)
    qcell[index+ncell_nf*i]=qval[i];
}

  
void CartBlock::preprocess(CartGrid *cg)
  {
    int nfrac;
    for(int n=0;n<3;n++) xlo[n]=cg->xlo[3*global_id+n];
    for(int n=0;n<3;n++) dx[n]=cg->dx[3*global_id+n];
    dims[0]=cg->ihi[3*global_id]  -cg->ilo[3*global_id  ]+1;
    dims[1]=cg->ihi[3*global_id+1]-cg->ilo[3*global_id+1]+1;
    dims[2]=cg->ihi[3*global_id+2]-cg->ilo[3*global_id+2]+1;
    nf=cg->nf;
    myid=cg->myid;
    donor_frac=cg->donor_frac;
    d1=dims[0];
    d2=dims[0]*dims[1];
    d3=d2*dims[2];
    ncell=d3;
    ncell_nf=(dims[0]+2*nf)*(dims[1]+2*nf)*(dims[2]+2*nf);
    nnode=(d1+1)*(d2+1)*(d3+1);
    nnode_nf=(dims[0]+1+2*nf)*(dims[1]+1+2*nf)*(dims[2]+1+2*nf);
  };

void CartBlock::initializeLists(void)
{
 donorList=(DONORLIST **)malloc(sizeof(DONORLIST *)*ncell);
 for(int i=0;i<ncell;i++) donorList[i]=NULL;
}

void CartBlock::clearLists(void)
{
  int i;
  if (donorList) {
  for(i=0;i<ncell;i++) { deallocateLinkList(donorList[i]); donorList[i]=NULL;}
  TIOGA_FREE(donorList);
  }
  deallocateLinkList4(interpList);
  interpList=NULL;
}


void CartBlock::insertInInterpList(int procid,int remoteid,int remoteblockid,double *xtmp)
{
  int i,n;
  int ix[3];
  double *rst;
  rst=(double *)malloc(sizeof(double)*3);
  if (interpList==NULL) 
    {
      interpList=(INTERPLIST2 *)malloc(sizeof(INTERPLIST2));
      listptr=interpList;
    }
  else
    {
      listptr->next=(INTERPLIST2 *)malloc(sizeof(INTERPLIST2));
      listptr=listptr->next;
    }
  listptr->next=NULL;
  listptr->inode=NULL;
  listptr->weights=NULL;
  listptr->receptorInfo[0]=procid;
  listptr->receptorInfo[1]=remoteid;
  listptr->receptorInfo[2]=remoteblockid;
  for(n=0;n<3;n++)
    {
      ix[n]=(xtmp[n]-xlo[n])/dx[n];
      rst[n]=(xtmp[n]-xlo[n]-ix[n]*dx[n])/dx[n];
      if (ix[n]==dims[n])
       {
        if (fabs(rst[n]) < TOL)
         {
         ix[n]--;
          rst[n]=(xtmp[n]-xlo[n]-ix[n]*dx[n])/dx[n];
         }
       }
      // if (!(ix[n] >=0 && ix[n] < dims[n]) && myid==77) {
      //  TRACEI(procid);
      //  TRACEI(global_id);
      //  TRACEI(local_id);
      //  TRACEI(remoteid);
      //  TRACEI(myid);
      //  TRACED(xtmp[0]);
      //  TRACED(xtmp[1]);
      //  TRACED(xtmp[2]);
      //  TRACED(xlo[0]);
      //  TRACED(xlo[1]);
      //  TRACED(xlo[2]);
      //  TRACED(dx[0]);
      //  TRACED(dx[1]);
      //  TRACED(dx[2]);
      //  TRACEI(ix[n]);
      //  TRACEI(n);
      //  TRACEI(dims[n]);
      //  printf("--------------------------\n");
      // }
     assert((ix[n] >=0 && ix[n] < dims[n]));
    }
  if (donor_frac == nullptr) {
    listptr->nweights=8;
    listptr->weights=(double *)malloc(sizeof(double)*listptr->nweights);
    listptr->inode=(int *)malloc(sizeof(int)*(listptr->nweights*3));
    cart_interp::linear_interpolation(nf,ix,dims,rst,&(listptr->nweights),
      listptr->inode,listptr->weights);
  }
  TIOGA_FREE(rst);
}
  
void CartBlock::insertInDonorList(int senderid,int index,int meshtagdonor,int remoteid,int remoteblockid, double cellRes)
{
  if(index >= ncell_nf)
    return;

  DONORLIST *temp1;
  int i,j,k,x_stride,xy_stride;
  int pointid;
  temp1=(DONORLIST *)malloc(sizeof(DONORLIST));

  // Get point-id accounting for nf
  x_stride = (dims[0]+2*nf);
  xy_stride = x_stride*(dims[1]+2*nf);
  k = index / xy_stride;
  index %= xy_stride;
  j = index / x_stride;
  index %= x_stride;
  i = index;
  pointid=(k-nf)*d2+(j-nf)*d1+(i-nf);

  if (!(pointid >= 0 && pointid < ncell)) {
    TRACEI(index);
    TRACEI(nf);
    TRACEI(dims[0]);
    TRACEI(dims[1]);
    TRACEI(dims[2]);
  }
  assert((pointid >= 0 && pointid < ncell));
    
  temp1->donorData[0]=senderid;
  temp1->donorData[1]=meshtagdonor;
  temp1->donorData[2]=remoteid;
  temp1->donorData[3]=remoteblockid;
  temp1->donorRes=cellRes;
  temp1->cancel=0;
  insertInList(&(donorList[pointid]),temp1);
}

void CartBlock::processDonors(HOLEMAP *holemap, int nmesh)
{
  int i,j,k,l,m,n,h;
  int ibcount,idof,meshtagdonor,icount;
  DONORLIST *temp;
  int *iflag;
  int holeFlag;
  double *xtmp;
  int ibindex;
  //FILE*fp;
  char fname[80];
  char qstr[2];
  char intstring[7];
  int ni,nj,nk,ibcheck;
  //sprintf(intstring,"%d",100000+myid);
  //sprintf(fname,"fringes_%s.dat",&(intstring[1]));
  //if (local_id==0) 
  //  {
  //    fp=fopen(fname,"w");
  //  }
  //else
  //  {
  //    fp=fopen(fname,"a");
  //  }  

  //
  // first mark hole points
  //
  iflag=(int *)malloc(sizeof(int)*nmesh);
  xtmp=(double *)malloc(sizeof(double)*3);
  //
  ibcount=-1;
  idof=-1;
  //for(i=0;i<(dims[0]+2*nf)*(dims[1]+2*nf)*(dims[2]+2*nf);i++) ibl[i]=1;

  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	{
	  ibcount++;
    xtmp[0] = xlo[0] + (i+0.5)*dx[0];
    xtmp[1] = xlo[1] + (j+0.5)*dx[1];
    xtmp[2] = xlo[2] + (k+0.5)*dx[2];

          holeFlag=1;
          idof=ibcount-1;
	    {
	      idof++;
	      if (donorList[idof]==NULL)
		{
		  for(h=0;h<nmesh;h++)
		    if (holemap[h].existWall)
		      {
			if (checkHoleMap(xtmp,holemap[h].nx,holemap[h].sam,holemap[h].extents))
			  {
	                    ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
			    ibl[ibindex]=0;
                            holeFlag=0;
			    break;
			  }
		      }
		}
	      else
		{
		  temp=donorList[idof];
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
			    if (checkHoleMap(xtmp,holemap[h].nx,holemap[h].sam,holemap[h].extents))
			      {
	                        ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
				ibl[ibindex]=0;
                                holeFlag=0;
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
          ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;

	  if (ibl[ibindex]==0) 
	    {
              idof=ibcount-1;
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
              idof=ibcount-1;
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
	      if (icount==1)
		{
		  ibl[ibindex]=-1;
		}
	      else
		{
                  idof=ibcount-1;
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

  
  /*
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
        {
          ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
          if (ibl[ibindex]==1)
            {
              ibcheck=1;
              for(nk=-1;nk<2;nk++)
                for(nj=-1;nj<2;nj++)
                  for(ni=-1;ni<2;ni++)
                    {
                      ibindex=(k+nf+nk)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf+nj)*(dims[0]+2*nf)+i+nf+ni;
                      if ((ibindex < 0) || (ibindex >= dims[0]*dims[1]*dims[2])) continue;
                      ibcheck=ibcheck && (ibl[ibindex]!=0);
                    }
              if (!ibcheck)
                {
                  printf("fixing orphan: myid/globalid/localid/(i,j,k)=%d %d %d %d %d %d \n",
                         myid,global_id,local_id,i,j,k);
		  ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
		  ibl[ibindex]=0;
                }
            }
        }
  */

  if (iflag) TIOGA_FREE(iflag);
  if (xtmp)  TIOGA_FREE(xtmp);
  // fclose(fp);
}
			      

void CartBlock::getCancellationData(int *cancelledData, int *ncancel)
{
  int i,j,k,m;
  int idof;
  DONORLIST *temp;
  idof=-1;
  m=0;
  *ncancel=0;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	  {
	    idof++;
	    if (donorList[idof]!=NULL)
	      {
		temp=donorList[idof];
		while(temp!=NULL)
		  {
		    if (temp->cancel==1) {
		      (*ncancel)++;
		      cancelledData[m++]=temp->donorData[0];
		      cancelledData[m++]=1;
		      cancelledData[m++]=temp->donorData[2];
		      cancelledData[m++]=temp->donorData[3];
		    }
	            temp=temp->next;
		  }
	      }
	  }
}


void CartBlock::writeCellFile(int bid)
{
  int ibmin,ibmax;
  char fname[80];
  char qstr[2];
  char intstring[7];
  char hash,c;
  int i,n,j,k,ibindex;
  int bodytag;
  FILE *fp;
  int ba,id;
  int nvert;
  int nnodes,ncells;
  int dd1,dd2;
  
  ibmin=30000000;
  ibmax=-30000000;
  nnodes=(dims[1]+1)*(dims[0]+1)*(dims[2]+1);
  ncells=dims[0]*dims[1]*dims[2];
  sprintf(intstring,"%d",100000+myid);
  sprintf(fname,"cart_cell%s.dat",&(intstring[1]));
  if (bid==0) 
    {
      fp=fopen(fname,"w");
    }
  else
    {
      fp=fopen(fname,"a");
    }
  if (bid==0) {
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK_CELL\" ");
  fprintf(fp,"\n");
  }
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,
	  ncells);
  fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");

  for(k=0;k<dims[2]+1;k++)
    for(j=0;j<dims[1]+1;j++)
      for(i=0;i<dims[0]+1;i++)
	fprintf(fp,"%lf\n",xlo[0]+dx[0]*i);
  for(k=0;k<dims[2]+1;k++)
    for(j=0;j<dims[1]+1;j++)
      for(i=0;i<dims[0]+1;i++)
	fprintf(fp,"%lf\n",xlo[1]+dx[1]*j);
  for(k=0;k<dims[2]+1;k++)
    for(j=0;j<dims[1]+1;j++)
      for(i=0;i<dims[0]+1;i++)
	fprintf(fp,"%lf\n",xlo[2]+dx[2]*k);
		
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	{
	  ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+
	    (j+nf)*(dims[0]+2*nf)+(i+nf);
          ibmin=TIOGA_MIN(ibmin,ibl[ibindex]);
          ibmax=TIOGA_MAX(ibmax,ibl[ibindex]);
	  fprintf(fp,"%d\n", ibl[ibindex]);
	}

  //printf("proc %d , block %d, ibmin/ibmax=%d %d\n",myid,bid,ibmin,ibmax);
  id=0;
  dd1=(dims[0]+1);
  dd2=dd1*(dims[1]+1);
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
        {
	  id=k*dd2+j*dd1+i+1;
          fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		 id,
		 id+1,
		 id+1+dd1,
		 id+dd1,
		 id+dd2,
		 id+1+dd2,
		 id+1+dd1+dd2,
		 id+dd1+dd2);
	}
                                               
  fclose(fp);
  return;
}
