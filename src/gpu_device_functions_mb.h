#ifndef GPU_DEVICE_FUNCTIONS_MB_H
#define GPU_DEVICE_FUNCTIONS_MB_H

#include "tioga_gpu.h"

#include <cmath>

#define NEQNS 3
#define maxStackSize 2048
#define d_fabs(a) (a > 0 ? a:-a)

// leave this as pifus native for now for 
// posterity sake, fix later
// 
//# define TIOGA_DEVICE_MIN(x,y)  (x) < (y) ? (x) : (y)
//# define TIOGA_DEVICE_MAX(x,y)  (x) > (y) ? (x) : (y)
# define TIOGA_DEVICE_BIGVALUE 1e15
# define TIOGA_DEVICE_TOL 1e-10
# define TIOGA_DEVICE_BASE 1

TIOGA_GPU_DEVICE
void d_invertMat3(double A[3][3],
  double f[3])
{
  double b11,b21,b22,b31,b32,b33;
  double u12,u13,u23;
  double d1,d2,d3;

  b11=1./A[0][0];
  u12=A[0][1]*b11;
  u13=A[0][2]*b11;
  b21=A[1][0];
  b22=1./(A[1][1]-b21*u12);
  u23=(A[1][2]-b21*u13)*b22;
  b31=A[2][0];
  b32=A[2][1]-b31*u12;
  b33=1./(A[2][2]-b31*u13-b32*u23);
  d1=f[0]*b11;
  d2=(f[1]-b21*d1)*b22;
  d3=(f[2]-b31*d1-b32*d2)*b33;
  f[2]=d3;
  f[1]=d2-u23*f[2];
  f[0]=d1-u12*f[1]-u13*f[2];
}

TIOGA_GPU_DEVICE
void d_solvec(double a[NEQNS][NEQNS],double b[NEQNS],int *iflag,int n)
{
  int i,j,k,l,flag;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;

  for(i=0;i<n;i++)
  {
    if (d_fabs(a[i][i]) < eps)
    {
      flag=1;
      for(k=i+1;k<n && flag;k++)
      {
        if (a[k][i]!=0)
        {
          flag=0;
          for(l=0;l<n;l++)
          {
            temp=a[k][l];
            a[k][l]=a[i][l];
            a[i][l]=temp;
          }
          temp=b[k];
          b[k]=b[i];
          b[i]=temp;
        }
      }
      if (flag) {*iflag=0;return;}
    }
    for(k=i+1;k<n;k++)
    {
      if (i!=k)
      {
        fact=-a[k][i]/a[i][i];
        for(j=0;j<n;j++)
        {
          a[k][j]+=fact*a[i][j];
        }
        b[k]+=fact*b[i];
      }
    }
  }

  for(i=n-1;i>=0;i--)
  {
    sum=0;
    for(j=i+1;j<n;j++) sum+=a[i][j]*b[j];
    b[i]=(b[i]-sum)/a[i][i];
  }

  *iflag=1;
  return;
}

TIOGA_GPU_DEVICE
void d_newtonSolve(double f[7][3],double *u1,double *v1,double *w1)
{
  int j;
  int iter,itmax,isolflag;
  double u,v,w;
  double uv,wu,vw,uvw,norm,convergenceLimit;
  double alph;
  double lhs[3][3];
  double rhs[3];
  //
  itmax=500;
  convergenceLimit=1e-14;
  alph=1.0;
  isolflag=1.0;
  //
  u=v=w=0.5;
  //
  for(iter=0;iter<itmax;iter++)
  {
    uv=u*v;
    vw=v*w;
    wu=w*u;
    uvw=u*v*w;

    for(j=0;j<3;j++)
      rhs[j]=f[0][j]+f[1][j]*u+f[2][j]*v+f[3][j]*w+
      f[4][j]*uv + f[5][j]*vw + f[6][j]*wu +
      f[7][j]*uvw;

    norm=rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2];
    if (sqrt(norm) <= convergenceLimit) break;

    for(j=0;j<3;j++)
    {
      lhs[j][0]=f[1][j]+f[4][j]*v+f[6][j]*w+f[7][j]*vw;
      lhs[j][1]=f[2][j]+f[5][j]*w+f[4][j]*u+f[7][j]*wu;
      lhs[j][2]=f[3][j]+f[6][j]*u+f[5][j]*v+f[7][j]*uv;
    }

    d_invertMat3(lhs,rhs);
    if (isolflag==0) break;

    u-=(rhs[0]*alph);
    v-=(rhs[1]*alph);
    w-=(rhs[2]*alph);
  }
  if (iter==itmax) {u=2.0;v=w=0.;}
  if (isolflag==0) {
    u=2.0;
    v=w=0.;
  }
  *u1=u;
  *v1=v;
  *w1=w;
  return;
}

TIOGA_GPU_DEVICE
void d_computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert)
{
  int j,k,isolflag;
  double f[8][3];
  double u,v,w;
  double oneminusU,oneminusV,oneminusW,oneminusUV;
  double lhs[3][3];
  double rhs[3];

  switch(nvert)
  {
    case 4:
      //
      // tetrahedron
      //
      for(k=0;k<3;k++)
      {
        for(j=0;j<3;j++)
          lhs[j][k]=xv[k][j]-xv[3][j];
        rhs[k]=xp[k]-xv[3][k];
      }
      //
      // invert the 3x3 matrix
      //
      d_solvec(lhs,rhs,&isolflag,3);
      //
      // check if the solution is not degenerate
      //
      if (isolflag)
      {
        for(k=0;k<3;k++) frac[k]=rhs[k];
        frac[3]=1.-frac[0]-frac[1]-frac[2];
      }
      else
      {
        frac[0]=1.0;
        frac[1]=frac[2]=frac[3]=0;
      }
      break;
    case 5:
      //
      // pyramid
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[3][j]-xv[0][j];
        f[3][j]=xv[4][j]-xv[0][j];
        //
        f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
        f[5][j]=xv[0][j]-xv[3][j];
        f[6][j]=xv[0][j]-xv[1][j];
        f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j];
      }
      //
      d_newtonSolve(f,&u,&v,&w);
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=w;
      //
      break;
    case 6:
      //
      // prizm
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[2][j]-xv[0][j];
        f[3][j]=xv[3][j]-xv[0][j];
        //
        f[4][j]=0;
        f[5][j]=xv[0][j]-xv[2][j]-xv[3][j]+xv[5][j];
        f[6][j]=xv[0][j]-xv[1][j]-xv[3][j]+xv[4][j];
        f[7][j]=0.;
      }
      //
      d_newtonSolve(f,&u,&v,&w);
      //
      oneminusUV=1.0-u-v;
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusUV*oneminusW;
      frac[1]=u*oneminusW;
      frac[2]=v*oneminusW;
      frac[3]=oneminusUV*w;
      frac[4]=u*w;
      frac[5]=v*w;
      //
      break;
    case 8:
      //
      // hexahedra
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[3][j]-xv[0][j];
        f[3][j]=xv[4][j]-xv[0][j];
        //
        f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
        f[5][j]=xv[0][j]-xv[3][j]+xv[7][j]-xv[4][j];
        f[6][j]=xv[0][j]-xv[1][j]+xv[5][j]-xv[4][j];
        f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j]+
            xv[4][j]-xv[5][j]+xv[6][j]-xv[7][j];
      }
      //
      d_newtonSolve(f,&u,&v,&w);
      //
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;
      //
      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=oneminusU*oneminusV*w;
      frac[5]=u*oneminusV*w;
      frac[6]=u*v*w;
      frac[7]=oneminusU*v*w;
      //
      break;
    default:
      // printf("Interpolation not implemented for polyhedra with %d vertices\n",nvert);
      break;
  }
}

TIOGA_GPU_DEVICE
void d_checkContainment(double *x, int **vconn,int *nc, int *nv, int ntypes, int *elementList,
  int cellIndex[2], int adtElement, double *xsearch)
{
  int i,j,m,n,i3;
  int nvert;
  int icell;
  int isum;
  double xv[8][3];
  double frac[8];
  int ihigh=0;
  //
  icell=elementList[adtElement];
  if (ihigh==0)
  {
    //
    // locate the type of the cell
    //
    isum=0;
    for(n=0;n<ntypes;n++)
    {
      isum+=nc[n];
      if (icell < isum)
      {
        i=icell-(isum-nc[n]);
        break;
      }
    }
    //
    // now collect all the vertices in the
    // array xv
    //
    nvert=nv[n];
    for(m=0;m<nvert;m++)
    {
      i3=3*(vconn[n][nvert*i+m]-TIOGA_DEVICE_BASE);
      for(j=0;j<3;j++)
        xv[m][j]=x[i3+j];
    }
    //
    d_computeNodalWeights(xv,xsearch,frac,nvert);
    //
    cellIndex[0]=icell;
    cellIndex[1]=0;
    //
    // if any of the nodal weights are
    // not in between [-TOL 1+TOL] discard
    // cell
    //
    for(m=0;m<nvert;m++)
    {
      if ((frac[m]+TIOGA_DEVICE_TOL)*(frac[m]-1.0-TIOGA_DEVICE_TOL) > 0)
      {
        cellIndex[0]=-1;
        return;
      }
      // need to find another way to implement this
      //if (fabs(frac[m]) < TOL && cellRes[icell]==TIOGA_DEVICE_BIGVALUE) cellIndex[1]=1;
    }

    return;
  }
  else
  {
    // not implemented;
  }
}

TIOGA_GPU_DEVICE
void d_searchIntersections_containment(int idx, int cellIndex[2],
  int *adtIntegers,
  double *adtReals,
  double *coord,
  double *xsearch,
  double *x, int *nc, int *nv,int ntypes,
  int *elementList,
  int **vconn,
  int nelem,
  int ndim,
  int *nchecks)
{
  int i,is;
  int d,nodeChild,flag;
  double element[6];
  int nodeStack[maxStackSize];
  int mm=0;
  int nstack=1;
  int bruteSearch=0;
  int node=0;
  int level=0;

  nodeStack[0]=node;
  cellIndex[0]=cellIndex[1]=-1;
  while(nstack > 0 && !bruteSearch)
  {
    mm=0;
    for(is=0;is<nstack;is++)
    {
      node=nodeStack[is];
      for(i=0;i<ndim;i++)
        element[i]=coord[ndim*(adtIntegers[4*node])+i];
      //
      flag=1;
      for(i=0;i<ndim/2;i++)
        flag = (flag && (xsearch[i] >=element[i]-TIOGA_DEVICE_TOL));
      for(i=ndim/2;i<ndim;i++)
        flag = (flag && (xsearch[i-ndim/2] <=element[i]+TIOGA_DEVICE_TOL));
      if (flag)
      {
        d_checkContainment(x,vconn,nc,nv,ntypes,elementList,
          cellIndex,adtIntegers[4*node],xsearch);
        if (cellIndex[0] > -1 && cellIndex[1]==0) return;
      }
      //
      // check the left and right children
      // now and sort based on distance if it is
      // within the given octant
      //
      for(d=1;d<3;d++)
      {
        nodeChild=adtIntegers[4*node+d];
        if (nodeChild > -1) {
          nodeChild=adtIntegers[4*nodeChild+3];
          for(i=0;i<ndim;i++)
          {
            element[i]=adtReals[ndim*nodeChild+i];
          }
          flag=1;
          for(i=0;i<ndim/2;i++)
            flag = (flag && (xsearch[i] >=element[i]-TIOGA_DEVICE_TOL));
          for(i=ndim/2;i<ndim;i++)
            flag = (flag && (xsearch[i-ndim/2] <=element[i]+TIOGA_DEVICE_TOL));
          if (flag)
          {
            if (mm+nstack > maxStackSize-1)  {
              cellIndex[0]=cellIndex[1]=-1;
              bruteSearch=1;
            }
            else {
              nodeStack[mm+nstack]=nodeChild;
              (*nchecks)=(mm+nstack > (*nchecks)) ? (mm+nstack):(*nchecks);
              mm++;
            }
          }
        }
      }
    }
    if (!bruteSearch) {
      for(int j=0;j<mm;j++)
        nodeStack[j]=nodeStack[j+nstack];
      nstack=mm;
      level++;
    }
  }
  if (bruteSearch) {
    for(node=0;node<nelem;node++)
    {
      for(i=0;i<ndim;i++)
        element[i]=coord[ndim*(adtIntegers[4*node])+i];
      //
      flag=1;
      for(i=0;i<ndim/2;i++)
        flag = (flag && (xsearch[i] >=element[i]-TIOGA_DEVICE_TOL));
      for(i=ndim/2;i<ndim;i++)
        flag = (flag && (xsearch[i-ndim/2] <=element[i]+TIOGA_DEVICE_TOL));
      if (flag)
      {
        d_checkContainment(x,vconn,nc,nv,ntypes,elementList,
          cellIndex,adtIntegers[4*node],xsearch);
        if (cellIndex[0] > -1 && cellIndex[1]==0) return;
      }
    }
  }
  return;
}

#endif /* GPU_DEVICE_FUNCTIONS_MB_H */
