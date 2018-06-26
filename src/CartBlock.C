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
#include "funcs.hpp"
#include "linklist.h"
#include "math_funcs.h"
#include "utils.h"

#include <assert.h>

extern "C" {
  void get_amr_index_xyz( int nq,int i,int j,int k,
			  int pBasis,
			  int nX,int nY,int nZ,
			  int nf,
			  double xlo[3],double dx[3],
			  double qnodes[],
			  int* index, double* xyz);
    void amr_index_to_ijklmn(int pBasis,int nX,int nY,int nZ, int nf, int nq,
			     int index, int* ijklmn);
}

void CartBlock::getInterpolatedData(int *nints,int *nreals,int **intData,
				    double **realData,
				    int nvar,int itype)
{
  int i,n,npnode;
  int ploc;
  double *xtmp;
  int *index;
  double qq[7];
  int *tmpint;
  double *tmpreal;
  int icount,dcount;
  int nintold,nrealold;
  int interpCount=0;
  int ix[3];

  listptr=interpList;
  while(listptr!=NULL)
    {
      interpCount++;
      listptr=listptr->next;
    }

  if (interpCount > 0) 
  {
    // Copy over current data before resizing allocation
    nintold=(*nints);
    nrealold=(*nreals);
    if (nintold > 0) {
      tmpint=(int *)malloc(sizeof(int)*3*(*nints));
      tmpreal=(double *)malloc(sizeof(double)*(*nreals));
      for(i=0;i<(*nints)*3;i++) tmpint[i]=(*intData)[i];
      for(i=0;i<(*nreals);i++) tmpreal[i]=(*realData)[i];

      free(*intData);
      free(*realData);
    }

    // Resize data arrays
    (*nints)+=interpCount;
    (*nreals)+=(interpCount*nvar);
    (*intData)=(int *)malloc(sizeof(int)*3*(*nints));
    (*realData)=(double *)malloc(sizeof(double)*(*nreals));

    // Copy over existing data
    if (nintold > 0) {
      for(i=0;i<nintold*3;i++) (*intData)[i]=tmpint[i];
      for(i=0;i<nrealold;i++) (*realData)[i]=tmpreal[i];
      free(tmpint);
      free(tmpreal);
    }

    listptr=interpList;
    icount=3*nintold;
    dcount=nrealold;

    index = (int *)malloc(sizeof(int)*listptr->nweights);
    ploc=(pdegree)*(pdegree+1)/2;

    while(listptr!=NULL)
    {
      if (itype==0)
      {
	      xtmp = (double *)malloc(sizeof(double)*p3*3);
        get_amr_index_xyz(qstride,listptr->inode[0],listptr->inode[1],listptr->inode[2],
            pdegree,dims[0],dims[1],dims[2],nf,
            xlo,dx,&qnode[ploc],index,xtmp);
	      free(xtmp);
      }
      else
      {
        // inode contains index of lower-left corner of Lagrange reconstruction
        const int I = listptr->inode[0];
        const int J = listptr->inode[1];
        const int K = listptr->inode[2];

        // Shorthands for reconstructing node index
        const int N = pdegree + 1;
        const int dj = (dims[0]+2*nf);
        const int dk = (dims[0]+2*nf)*(dims[1]+2*nf);
        const int ind = I + nf + dj*(J+nf) + dk*(K+nf);

        // Grab the cube of p3 node indices starting from I,J,K
        int m = 0;
        for (int k = 0; k < N; k++)
          for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++)
            {
              index[m] = ind + i + dj*j + dk*k;
              m++;
            }
      }

      (*intData)[icount++]=listptr->receptorInfo[0];
      (*intData)[icount++]=-1;
      (*intData)[icount++]=listptr->receptorInfo[1];

      for (int n = 0; n < nvar; n++)
      {
        qq[n] = 0;
        for (int i = 0; i < listptr->nweights; i++)
        {
          const double weight = listptr->weights[i];
          qq[n] += q[index[i]+d3nf*n] * weight;
        }
      }

      //writeqnode_(&myid,qq,&nvar);
      for(int n=0;n<nvar;n++)
        (*realData)[dcount++]=qq[n];
      listptr=listptr->next;
    }
    free(index);
  }
}


void CartBlock::update(double *qval, int index,int nq,int itype)
{
  int i;
  if (itype==0) 
   {
    for(i=0;i<nq;i++)
      q[index+d3nf*i]=qval[i];
   }
  else
   {
     int itmp=index;
     int k=itmp/(dims[1]*dims[0]);
     itmp=itmp%(dims[1]*dims[0]);
     int j=itmp/dims[0];
     int i=itmp%dims[0];
     itmp=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
     for(i=0;i<nq;i++)
       q[itmp+d3nf*i]=qval[i];
   }
}

  
void CartBlock::preprocess(CartGrid *cg,int itype)
{
  int nfrac;
  FILE *fp;
  char intstring[7];
  char fname[20];
  myid = cg->myid;

  for(int n=0;n<3;n++) xlo[n]=cg->xlo[3*global_id+n];
  for(int n=0;n<3;n++) dx[n]=cg->dx[3*global_id+n];
  dims[0] = cg->ihi[3*global_id]  -cg->ilo[3*global_id  ]+1;
  dims[1] = cg->ihi[3*global_id+1]-cg->ilo[3*global_id+1]+1;
  dims[2] = cg->ihi[3*global_id+2]-cg->ilo[3*global_id+2]+1;

  pdegree = cg->porder[global_id];
  p3 = (pdegree+1)*(pdegree+1)*(pdegree+1);
  nf = cg->nf;
  qstride = cg->qstride;
  donor_frac = cg->donor_frac;
  qnode = cg->qnode;
  d1 = dims[0];
  d2 = dims[0]*dims[1];
  d3 = d2*dims[2];
  d3nf =(dims[0]+2*nf)*(dims[1]+2*nf)*(dims[2]+2*nf);

  free(ibstore);
  ibstore = (int *)malloc(sizeof(int)*d3nf);
  for (int i=0;i<d3nf;i++) ibstore[i] = ibl[i];

  ndof = (itype==0) ? (d3*p3):d3;
};

void CartBlock::initializeLists(void)
{
  donorList = (DONORLIST **)malloc(sizeof(DONORLIST *)*ndof);
  for(int i=0;i<ndof;i++) donorList[i]=NULL;
}

void CartBlock::clearLists(void)
{
  if (donorList)
  {
    for (int i = 0; i < ndof; i++)
    {
      deallocateLinkList(donorList[i]);
      donorList[i] = NULL;
    }
    free(donorList);
    donorList = NULL;
  }

  deallocateLinkList4(interpList);
  interpList = NULL;
}


void CartBlock::insertInInterpList(int procid,int remoteid,double *xtmp,int itype)
{
  int i,j,k,n;
  int i1,j1,k1;
  int ix[3];
  int index[8];
  double *rst = (double *)malloc(sizeof(double)*3);
  if (interpList==NULL) 
  {
    interpList = (INTERPLIST2 *)malloc(sizeof(INTERPLIST2));
    listptr = interpList;
  }
  else
  {
    listptr->next = (INTERPLIST2 *)malloc(sizeof(INTERPLIST2));
    listptr = listptr->next;
  }

  listptr->next=NULL;
  listptr->inode=NULL;
  listptr->weights=NULL;
  listptr->receptorInfo[0]=procid;
  listptr->receptorInfo[1]=remoteid;

  for (int n = 0; n < 3; n++)
  {
    ix[n] = (xtmp[n] - xlo[n]) / dx[n];
    rst[n] = (xtmp[n] - xlo[n] - ix[n]*dx[n]) / dx[n];
      // if (!(ix[n] >=0 && ix[n] < dims[n]) && myid==77) {
      //  tracei(procid);
      //  tracei(global_id);
      //  tracei(local_id);
      //  tracei(remoteid);
      //  tracei(myid);
      //  traced(xtmp[0]);
      //  traced(xtmp[1]);
      //  traced(xtmp[2]);
      //  traced(xlo[0]);
      //  traced(xlo[1]);
      //  traced(xlo[2]);
      //  traced(dx[0]);
      //  traced(dx[1]);
      //  traced(dx[2]);
      //  tracei(ix[n]);
      //  tracei(n);
      //  tracei(dims[n]);
      //  printf("--------------------------\n");
      // }
     assert((ix[n] >=0 && ix[n] < dims[n]));
  }

  listptr->inode=(int *)malloc(sizeof(int)*3);
  listptr->inode[0]=ix[0];
  listptr->inode[1]=ix[1];
  listptr->inode[2]=ix[2];

  if (itype==0) {
    listptr->nweights = p3;
    listptr->weights=(double *)malloc(sizeof(double)*listptr->nweights);
    donor_frac(&pdegree,rst,&(listptr->nweights),(listptr->weights));  
  }
  else
  {
    /// TODO: maybe write separate function? (need another to get indices as well)
    listptr->nweights = p3;
    listptr->weights = (double *)malloc(sizeof(double)*listptr->nweights);

    // -------------------- from cart_interp.cpp --------------------
    const int N = pdegree + 1;
    double dn = 2. / pdegree;

    // Create a local 'reference element' and perform Lagrange interpolation
    // Element is centered around node i0,j0,k0
    unsigned int I = std::min(std::max(ix[0] - N/2 + 1, 0), dims[0]-N);
    unsigned int J = std::min(std::max(ix[1] - N/2 + 1, 0), dims[1]-N);
    unsigned int K = std::min(std::max(ix[2] - N/2 + 1, 0), dims[2]-N);

    // Change inode to point to the lower-left corner of our constructed element
    listptr->inode[0] = I;
    listptr->inode[1] = J;
    listptr->inode[2] = K;

    double xi  = -1. + 2. * (xtmp[0] - (xlo[0] + I*dx[0])) / (pdegree * dx[0]);
    double eta = -1. + 2. * (xtmp[1] - (xlo[1] + J*dx[1])) / (pdegree * dx[1]);
    double nu  = -1. + 2. * (xtmp[2] - (xlo[2] + K*dx[2])) / (pdegree * dx[2]);

    std::vector<double> xiGrid(N); // Equidistant grid from -1 to 1
    for (int i = 0; i < N; i++)
      xiGrid[i] = -1. + i * dn;

    std::vector<double> lag_i(N), lag_j(N), lag_k(N);
    for (int i = 0; i < N; i++)
    {
      lag_i[i] = tg_funcs::Lagrange(xiGrid.data(), N, xi, i);
      lag_j[i] = tg_funcs::Lagrange(xiGrid.data(), N, eta, i);
      lag_k[i] = tg_funcs::Lagrange(xiGrid.data(), N, nu, i);
    }

    // Shorthands for reconstructing node index
    unsigned int dj = dims[0]+2*nf;
    unsigned int dk = dj*(dims[1]+2*nf);
    unsigned int ind = (I+nf) + (J+nf)*dj + (K+nf)*dk;

    int m = 0;
    for (int k = 0; k < N; k++)
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
        {
          listptr->weights[m] =  lag_i[i] * lag_j[j] * lag_k[k];
          // Set iblank of this node to be mandatory donor
          ibl[ind + i + dj*j + dk*k] = -2;
          m++;
        }
    // --------------------------------------------------------------------
  }

  free(rst);
}
  
void CartBlock::insertInDonorList(int senderid,int index,int meshtagdonor,int remoteid,double cellRes,
                                 double receptorRes,int itype)
{
  int ijklmn[6];
  int pointid;

  DONORLIST* temp1 = (DONORLIST *)malloc(sizeof(DONORLIST));

  if (itype==0) 
  {
    amr_index_to_ijklmn(pdegree,dims[0],dims[1],dims[2],nf,qstride,index,ijklmn);
    pointid = (ijklmn[2]*d2 + ijklmn[1]*d1 + ijklmn[0])*p3 +
              ijklmn[5]*(pdegree+1)*(pdegree+1) + ijklmn[4]*(pdegree+1)+ijklmn[3];

    if (!(pointid >= 0 && pointid < ndof)) {
      tracei(index);
      tracei(nf);
      tracei(pdegree);
      tracei(dims[0]);
      tracei(dims[1]);
      tracei(dims[2]);
      tracei(qstride);
      printf("%d %d %d %d %d %d\n",ijklmn[0],ijklmn[1],ijklmn[2],
          ijklmn[3],ijklmn[4],ijklmn[5]);
    }
    assert((pointid >= 0 && pointid < ndof));
  }
  else
  {
    pointid = index;
  }
    
  temp1->donorData[0] = senderid;
  temp1->donorData[1] = meshtagdonor;
  temp1->donorData[2] = remoteid;
  temp1->donorRes     = cellRes;
  temp1->receptorRes  = receptorRes;
  temp1->cancel=0;
  assert((pointid >= 0 && pointid < ndof));
  insertInList(&(donorList[pointid]),temp1);
}

void CartBlock::processDonors(HOLEMAP *holemap, int nmesh,int itype)
{
  DONORLIST *temp;
  char fname[80];
  char qstr[2];
  char intstring[7];
  int ni,nj,nk,ibcheck;

  // first mark hole points
  int p3t = (itype==0)?p3:1;
  int* iflag = (int *)malloc(sizeof(int)*nmesh);
  int* index = (int *)malloc(sizeof(int)*p3t);
  double* xtmp = (double *)malloc(sizeof(double)*p3t*3);
  
  int ibcount = -1;
  int idof = -1;
  
  //for(i=0;i<(dims[0]+2*nf)*(dims[1]+2*nf)*(dims[2]+2*nf);i++) ibl[i]=1;

  int ploc = pdegree*(pdegree+1)/2;
  //iadd=(itype==0)?0:1;
  int iadd = 0;
   
  for (int k=0; k < dims[2]+iadd; k++)
  {
    for (int j=0; j < dims[1]+iadd; j++)
    {
      for (int i=0; i < dims[0]+iadd; i++)
      {
        ibcount++;
        if (itype==0) 
        {
          get_amr_index_xyz(qstride,i,j,k,pdegree,dims[0],dims[1],dims[2],nf,
              xlo,dx,&qnode[ploc],index,xtmp);
        }
        else 
        {
          xtmp[0]=xlo[0]+dx[0]*i;
          xtmp[1]=xlo[1]+dx[1]*j;
          xtmp[2]=xlo[2]+dx[2]*k;
        }

        int holeFlag = 1;
        idof = ibcount*p3t-1;
        for (int p = 0; p < p3t && holeFlag; p++)
        {
          idof++;
          if (donorList[idof]==NULL)
          {
            for (int h = 0; h < nmesh; h++)
            {
              if (holemap[h].existWall)
              {
                if (checkHoleMap(&xtmp[3*p],holemap[h].nx,holemap[h].sam,holemap[h].extents))
                {
                  int ibindex = (dims[0]+2*nf)*((dims[1]+2*nf)*(k+nf)+ j+nf) + i+nf;
                  // check hole-map
                  ibl[ibindex]=0;
                  holeFlag=0;
                  break;
                }
              }
            }
          }
          else
          {
            temp = donorList[idof];
            for (int h = 0; h < nmesh; h++) iflag[h]=0;
            while(temp!=NULL) 
            {
              int meshtagdonor=temp->donorData[1]-BASE;
              iflag[meshtagdonor]=1;
              temp=temp->next;
            }
            for (int h = 0; h < nmesh; h++)
            {
              if (holemap[h].existWall)
              {
                if (!iflag[h])
                  if (checkHoleMap(&xtmp[3*p],holemap[h].nx,holemap[h].sam,holemap[h].extents))
                  {
                    int dj = dims[0] + 2*nf;
                    int dk = dj * (dims[1] + 2*nf);
                    int ibindex = (k+nf)*dk + (j+nf)*dj + i+nf;
                    // check hole-map
                    ibl[ibindex]=0;
                    holeFlag=0;
                    break;
                  }
              }
            }
          }
        }
      }
    }
  }

  /// DEBUGGING / HACK FOR INCONSISTENT IBLANK VALUES BETWEEN PROCESSES
  /*for(k=0;k<dims[2]+iadd;k++)
    for(j=0;j<dims[1]+iadd;j++)
      for(i=0;i<dims[0]+iadd;i++)
	{
    int ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+(j+nf)*(dims[0]+2*nf)+i+nf;
    xtmp[0]=xlo[0]+dx[0]*i;
    xtmp[1]=xlo[1]+dx[1]*j;
    xtmp[2]=xlo[2]+dx[2]*k;
    double R = std::sqrt(xtmp[0]*xtmp[0] + xtmp[1]*xtmp[1] + xtmp[2]*xtmp[2]);
    if (R > 2.)
      ibl[ibindex] = -2; /// Force 'normal' status below
  }*/

  ibcount = -1;
  idof = -1;    
  for (int k = 0; k < dims[2]+iadd; k++)
  {
    for (int j = 0; j < dims[1]+iadd; j++)
    {
      for (int i = 0; i < dims[0]+iadd; i++)
      {
        int dj = dims[0] + 2*nf;
        int dk = dj * (dims[1] + 2*nf);
        int ibindex = (k+nf)*dk + (j+nf)*dj + i+nf;
        ibcount++;

        if (itype==0) 
        {
          get_amr_index_xyz(qstride,i,j,k,pdegree,dims[0],dims[1],dims[2],nf,
              xlo,dx,&qnode[ploc],index,xtmp);
        }

        if (ibl[ibindex] == 0)  // Hole node
        {
          idof = ibcount*p3t-1;
          for (int p = 0; p < p3t; p++)
          {
            idof++;
            if (donorList[idof]!=NULL)
            {
              temp = donorList[idof];
              while (temp != NULL)
              {
                temp->cancel = 1; // Cancel receptor status
                temp = temp->next;
              }
            }
          }
        }
        else if (ibl[ibindex] == -2)  // Part of a Cartesian donor cell [Mandatory Donor]
        {
          idof = ibcount*p3t-1;
          for (int p = 0; p < p3t; p++)
          {
            idof++;
            if (donorList[idof]!=NULL)
            {
              temp = donorList[idof];
              while (temp != NULL)
              {
                temp->cancel = 1; // Cancel receptor status
                temp = temp->next;
              }
            }
          }
          ibl[ibindex] = 1; //min(ibstore[ibindex],1);
        }
        else
        {
          int icount=0;
          idof=ibcount*p3t-1;
          for (int p = 0; p < p3t; p++)
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
          if (icount==p3t) 
          {
            if (ibstore[ibindex]==1) ibl[ibindex]=-1;
            //for(p=0;p<p3t;p++)
            // fprintf(fp,"%f %f %f\n",xtmp[3*p],xtmp[3*p+1],xtmp[3*p+2]);
          }
          else
          {
            idof=ibcount*p3t-1;
            for (int p = 0; p < p3t; p++)
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
  }

  if (itype==0) 
  {
    for (int k = 0; k < dims[2]+iadd; k++)
    {
      for (int j = 0; j < dims[1]+iadd; j++)
      {
        for (int i = 0; i < dims[0]+iadd; i++)
        {
          int dj = dims[0] + 2*nf;
          int dk = dj * (dims[1] + 2*nf);
          int ibindex = (k+nf)*dk + (j+nf)*dj + i+nf;
          if (ibl[ibindex]==1)
          {
            int ibcheck = 1;
            for (int nk = -1; nk < 2; nk++)
              for (int nj = -1; nj < 2; nj++)
                for (int ni = -1; ni < 2; ni++)
                {
                  ibindex = (k+nf+nk)*dk + (j+nf+nj)*dj + i+nf+ni;
                  ibcheck = ibcheck && (ibl[ibindex]!=0);
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
      }
    }
  }
  //for(i=0;i<d3nf;i++) ibl[i]=ibstore[i];
  free(xtmp);
  free(index);
  free(iflag);
  // fclose(fp);
}
			      

void CartBlock::getCancellationData(int *cancelledData, int *ncancel, int itype)
{
  int i,j,k,p,m,p3t;
  int idof;
  DONORLIST *temp;

  p3t=(itype==0)?p3:1;
  idof=-1;
  m=0;
  *ncancel=0;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	for(p=0;p<p3t;p++)
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
		    }
	            temp=temp->next;
		  }
	      }
	  }
}


void CartBlock::writeCellFile(int bid,int itype)
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
  int dd1,dd2,iadd;
  
  ibmin=30000000;
  ibmax=-30000000;
  if (itype==0) 
    {
      nnodes=(dims[1]+1)*(dims[0]+1)*(dims[2]+1);
      ncells=dims[0]*dims[1]*dims[2];
    }
  else
    {
      nnodes=dims[0]*dims[1]*dims[2];
      ncells=(dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    }
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
  if (itype==0) 
   {
    fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK_CELL\" ");
   }
  else
   {
   fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK \" ");
   }
  
  fprintf(fp,"\n");
  }
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,
	  ncells);
  if (itype==0) 
    {
      fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");
    }
  else
    {
      fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=NODAL)\n");
    }
  iadd=(itype==0)?1:0;

  for(k=0;k<dims[2]+iadd;k++)
    for(j=0;j<dims[1]+iadd;j++)
      for(i=0;i<dims[0]+iadd;i++)
	fprintf(fp,"%lf\n",xlo[0]+dx[0]*i);
  for(k=0;k<dims[2]+iadd;k++)
    for(j=0;j<dims[1]+iadd;j++)
      for(i=0;i<dims[0]+iadd;i++)
	fprintf(fp,"%lf\n",xlo[1]+dx[1]*j);
  for(k=0;k<dims[2]+iadd;k++)
    for(j=0;j<dims[1]+iadd;j++)
      for(i=0;i<dims[0]+iadd;i++)
	fprintf(fp,"%lf\n",xlo[2]+dx[2]*k);


  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
	{
	  ibindex=(k+nf)*(dims[1]+2*nf)*(dims[0]+2*nf)+
	    (j+nf)*(dims[0]+2*nf)+(i+nf);
          ibmin=min(ibmin,ibl[ibindex]);
          ibmax=max(ibmax,ibl[ibindex]);
          //if (ibindex > (dims[0]+2*nf)*(dims[1]+2*nf)*(dims[2]+2*nf)) {
          // printf("problem: %d %d %d\n",myid,bid,ibindex);
          //}
	  fprintf(fp,"%d\n", ibl[ibindex]);
	}

  //printf("proc %d , block %d, ibmin/ibmax=%d %d\n",myid,bid,ibmin,ibmax);
  id=0;
  dd1=(dims[0]+iadd);
  dd2=dd1*(dims[1]+iadd);
  for(k=0;k<dims[2]-1+iadd;k++)
    for(j=0;j<dims[1]-1+iadd;j++)
      for(i=0;i<dims[0]-1+iadd;i++)
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
