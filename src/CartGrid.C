// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

# include "codetypes.h"
# include "CartGrid.h"

void CartGrid::registerData(int nfin,int qstridein,double *qnodein,int *idata,
			    double *rdata,int ngridsin,int qnodesize)
{
  int i,i3,i6,iloc,n;
  FILE *fp;
  ngrids = ngridsin;
  global_id=(int *) malloc(sizeof(int)*ngrids);
  level_num=(int *) malloc(sizeof(int)*ngrids);
  proc_id=(int *) malloc(sizeof(int)*ngrids);
  ilo=(int *) malloc(sizeof(int)*3*ngrids);
  ihi=(int *) malloc(sizeof(int)*3*ngrids);
  xlo=(double *) malloc(sizeof(double)*3*ngrids);
  dx=(double *) malloc(sizeof(double)*3*ngrids);
  porder=(int *) malloc(sizeof(int)*ngrids);
  local_id=(int *)malloc(sizeof(int)*ngrids);
  qnode=(double *)malloc(sizeof(double)*qnodesize);
  dims=(int *)malloc(sizeof(dims)*3*ngrids);
  for(i=0;i<qnodesize;i++)  { qnode[i]=qnodein[i];}
  //                            if (myid==0) printf("qnode[%d]= %f\n",i,qnode[i]);}
  nf=nfin;
  qstride=qstridein;
  if (myid==0) fp=fopen("cartGrid.dat","w");
  for(i=0;i<ngrids;i++)
    {
      i3=3*i;
      i6=2*i3;
      iloc=11*i;

      global_id[i]=idata[iloc];
      level_num[i]=idata[iloc+1];
      proc_id[i]=idata[iloc+2];
      porder[i]=idata[iloc+3];
      local_id[i]=idata[iloc+4];
      for(n=0;n<3;n++)
	{
	  ilo[i3+n]=idata[iloc+5+n];
	  ihi[i3+n]=idata[iloc+8+n];
	  dims[i3+n]=ihi[i3+n]-ilo[i3+n]+1;
	}
      xlo[i3]=rdata[i6];
      xlo[i3+1]=rdata[i6+1];
      xlo[i3+2]=rdata[i6+2];
      dx[i3]=rdata[i6+3];
      dx[i3+1]=rdata[i6+4];
      dx[i3+2]=rdata[i6+5];
      if (myid==0) 
        fprintf(fp,"%d %d %d %d %d %f %f %f\n",global_id[i],level_num[i],proc_id[i],
                                   porder[i],local_id[i],dx[i3],dx[i3+1],
                                   dx[i3+2]);
    }
   if (myid==0) fclose(fp);
};

//
// Bare bone preprocessor now
// willl add more support data structures
// to promote efficient search once the concept
// works
//
void CartGrid::preprocess(void)
{
  int i,n;
  //
  // find the global minimum coord location
  //
  xlosup[0]=xlosup[1]=xlosup[2]=BIGVALUE;
  maxlevel=-1;
  for (i=0;i<ngrids;i++)
    {
      for(n=0;n<3;n++)
	xlosup[n]=TIOGA_Min(xlosup[n],xlo[3*i+n]);
      maxlevel=TIOGA_Max(maxlevel,level_num[i]);
    }
    maxlevel++;
  lcount=(int *)malloc(sizeof(int)*maxlevel);
  dxlvl=(double *)malloc(sizeof(double)*3*maxlevel);
  for(i=0;i<maxlevel;i++) lcount[i]=0;
  for(i=0;i<ngrids;i++)
    {
      lcount[level_num[i]]++;
      for(n=0;n<3;n++)
	dxlvl[3*level_num[i]+n]=dx[3*i+n];
    }
}
//
// Basic search routine now
// will improve efficiency once it works
//
void CartGrid::search(double *x,int *donorid,int npts)
{
  int i,j,k,l,n,il[3];
  bool flag;
  int dcount;
  dcount=0;
  for(i=0;i<npts;i++)
    {
      flag=0;
      donorid[i]=-1;
      for(l=maxlevel-1;l>=0 && flag==0;l--)
	{
	  for(n=0;n<3;n++)
	    il[n]=floor((x[3*i+n]-xlosup[n])/dxlvl[3*l+n]);
	  for(j=0;j<ngrids && flag==0;j++)
	    {
	      if (level_num[j]==l) 
		{
		  flag=1;
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] >=xlo[3*j+n]);
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] <=xlo[3*j+n]+
		  //   			 dx[3*j+n]*(dims[3*j+n]));
		  //for(n=0;n<3;n++) flag = flag && (il[n] >=ilo[3*j+n]);
		  //for(n=0;n<3;n++) flag = flag && (il[n] <=ihi[3*j+n]);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]-xlo[3*j+n]) > -TOL);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]- (xlo[3*j+n]+
                                                     dx[3*j+n]*(dims[3*j+n]))) < TOL);
		  if (flag) { 
		    dcount++; 
		    donorid[i]=j; 
		}
	    }
	}
    }
  }
 //printf("CartGrid::search Processor %d located %d of %d points\n",myid,dcount,npts);
}
