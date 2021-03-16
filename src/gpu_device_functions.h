
#include <cmath>

#define NEQNS 12
#define maxStackSize 256
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
int faceEdgeIntersectCheck(double *xface, // face coordinate array
			   double *xx, // search point
			   double *xc, // cell center
			   double *xi, // intersection pont
			   double EPS1, // tol. for ray 
			   double EPS2, // tol. in plane
			   int nvert)   // number of vertices
{
  double mat[3][3];
  double rhs[3];
  double sol[3];
  double ONE1=1+EPS1;
  double ONE2=1+EPS2;
  int i,j;
  double det,deti;
  int intersects;
//  FILE *fp;

//   fp=fopen("cell.dat","a");
//   fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=1 ET=QUADRILATERAL, F=FEPOINT\n",nvert);  
//   for(i=0;i<nvert;i++)
//     fprintf(fp,"%f %f %f\n",xface[i][0],xface[i][1],xface[i][2]);
//   fprintf(fp,"1 2 3 %d\n",(nvert==4)?nvert:3);
//   fclose(fp);
  intersects=0;
  for(i=1;i<nvert-1;i++)
    {
      
      for (j=0;j<3;j++)
	{
	  mat[j][0]=xc[j]-xx[j];
	  mat[j][1]=xface[i*3+j]-xface[0*3+j];
	  mat[j][2]=xface[(i+1)*3+j]-xface[0*3+j];
	  rhs[j]=xc[j]-xface[0*3+j];
	}
      //
      // invert
      // 
      det=( mat[0][0]*(mat[2][2]*mat[1][1]-mat[2][1]*mat[1][2])-
	    mat[1][0]*(mat[2][2]*mat[0][1]-mat[2][1]*mat[0][2])+
	    mat[2][0]*(mat[1][2]*mat[0][1]-mat[1][1]*mat[0][2]));

      if (d_fabs(det) > 1E-10) {
	deti=1./det;
	sol[0]=deti*
	  ( (mat[2][2]*mat[1][1]-mat[2][1]*mat[1][2])*rhs[0]
	    -(mat[2][2]*mat[0][1]-mat[2][1]*mat[0][2])*rhs[1]
	    +(mat[1][2]*mat[0][1]-mat[1][1]*mat[0][2])*rhs[2]);
	
	sol[1]=deti*
	  (-(mat[2][2]*mat[1][0]-mat[2][0]*mat[1][2])*rhs[0]
	   +(mat[2][2]*mat[0][0]-mat[2][0]*mat[0][2])*rhs[1]
	   -(mat[1][2]*mat[0][0]-mat[1][0]*mat[0][2])*rhs[2]);
	
	sol[2]=deti*
	  ((mat[2][1]*mat[1][0]-mat[2][0]*mat[1][1])*rhs[0]
	   -(mat[2][1]*mat[0][0]-mat[2][0]*mat[0][1])*rhs[1]
	   +(mat[1][1]*mat[0][0]-mat[1][0]*mat[0][1])*rhs[2]);
	
	if ((sol[0]+EPS1)*(sol[0]-ONE1) < 0) {
	  if ((sol[1]+EPS2) >= 0 && (sol[2]+EPS2) >=0 && 
	      (sol[1]+sol[2]) <= ONE2) 
	    {
	      xi[0]=xc[0]+sol[0]*(xx[0]-xc[0]);
	      xi[1]=xc[1]+sol[0]*(xx[1]-xc[1]);
	      xi[2]=xc[2]+sol[0]*(xx[2]-xc[2]);
	      intersects=1;
	    }			  
	}
      }
    }
  return intersects;
}

TIOGA_GPU_DEVICE
void invertMat3(double A[3][3],
                double f[3])
{
  int i;
  double b11,b21,b22,b31,b32,b33,b41,b42,b43,b44;
  double u12,u13,u23;
  double d1,d2,d3;
  double c1,c2,c3;
  double a1,a2,a3;
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
  return;
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
      for(j=i+1;j<n;j++)
	sum+=a[i][j]*b[j];
      b[i]=(b[i]-sum)/a[i][i];
    }
  *iflag=1;
  return;

}

TIOGA_GPU_DEVICE
void d_newtonSolve(double *f,double *u1,double *v1,double *w1)
{
  int i,j,k;
  int iter,itmax,isolflag;
  double u,v,w;
  double uv,wu,vw,uvw,norm,convergenceLimit;
  double alph;
  double lhs[9];
  double rhs[3],sol[3];
  double det,deti;
  //
  itmax=500;
  //itmax=0; // uncommenting this works, i.e. by pass this routine
  convergenceLimit=1e-8;
  alph=1.0;
  isolflag=1;
  norm=1.0;
  //
  u=v=w=0.5;
  //
  for(iter=0;iter<itmax && sqrt(norm) > convergenceLimit;iter++)
    {
      uv=u*v;
      vw=v*w;
      wu=w*u;
      uvw=u*v*w;
      
      for(j=0;j<3;j++)
	rhs[j]=f[0*3+j]+f[1*3+j]*u+f[2*3+j]*v+f[3*3+j]*w+
	  f[4*3+j]*uv + f[5*3+j]*vw + f[6*3+j]*wu +
	  f[7*3+j]*uvw;

      
      norm=rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2];
      //if (sqrt(norm) <= convergenceLimit) break;

      for(j=0;j<3;j++)
	{
	  lhs[j*3+0]=f[1*3+j]+f[4*3+j]*v+f[6*3+j]*w+f[7*3+j]*vw;
	  lhs[j*3+1]=f[2*3+j]+f[5*3+j]*w+f[4*3+j]*u+f[7*3+j]*wu;
	  lhs[j*3+2]=f[3*3+j]+f[6*3+j]*u+f[5*3+j]*v+f[7*3+j]*uv;
	}      

      det=( lhs[0*3+0]*(lhs[2*3+2]*lhs[1*3+1]-lhs[2*3+1]*lhs[1*3+2])-
	    lhs[1*3+0]*(lhs[2*3+2]*lhs[0*3+1]-lhs[2*3+1]*lhs[0*3+2])+
	    lhs[2*3+0]*(lhs[1*3+2]*lhs[0*3+1]-lhs[1*3+1]*lhs[0*3+2]));
      sol[0]=sol[1]=sol[2]=0;
      //det+=TIOGA_DEVICE_TOL;
      if (d_fabs(det) > 1E-10 ) { 
	deti=1.0/det;
 	sol[0]=deti*
	  ( (lhs[2*3+2]*lhs[1*3+1]-lhs[2*3+1]*lhs[1*3+2])*rhs[0]
	    -(lhs[2*3+2]*lhs[0*3+1]-lhs[2*3+1]*lhs[0*3+2])*rhs[1]
	    +(lhs[1*3+2]*lhs[0*3+1]-lhs[1*3+1]*lhs[0*3+2])*rhs[2]);

	sol[1]=deti*
	  (-(lhs[2*3+2]*lhs[1*3+0]-lhs[2*3+0]*lhs[1*3+2])*rhs[0]
	   +(lhs[2*3+2]*lhs[0*3+0]-lhs[2*3+0]*lhs[0*3+2])*rhs[1]
	   -(lhs[1*3+2]*lhs[0*3+0]-lhs[1*3+0]*lhs[0*3+2])*rhs[2]);
	
	sol[2]=deti*
	  ((lhs[2*3+1]*lhs[1*3+0]-lhs[2*3+0]*lhs[1*3+1])*rhs[0]
	   -(lhs[2*3+1]*lhs[0*3+0]-lhs[2*3+0]*lhs[0*3+1])*rhs[1]
	   +(lhs[1*3+1]*lhs[0*3+0]-lhs[1*3+0]*lhs[0*3+1])*rhs[2]);
      }
      else {
	isolflag=0;
      }
      //printf("%f\n",det);
      /*
      if (d_fabs(det) > 1E-10) {
	deti=1.0/det;
 	sol[0]=deti*
	  ( (lhs[2][2]*lhs[1][1]-lhs[2][1]*lhs[1][2])*rhs[0]
	    -(lhs[2][2]*lhs[0][1]-lhs[2][1]*lhs[0][2])*rhs[1]
	    +(lhs[1][2]*lhs[0][1]-lhs[1][1]*lhs[0][2])*rhs[2]);


	}
      //else {
      //isolflag=0;
      //}
      //d_solvec(lhs,rhs,&isolflag,3);
      //invertMat3(lhs,rhs);
      //if (isolflag==0) break;
      */
      u-=(sol[0]*alph);
      v-=(sol[1]*alph);
      w-=(sol[2]*alph);
    }
  if (iter==itmax) {u=2.0;v=w=0.;}
  if (u < -TIOGA_DEVICE_TOL || u > 1.0+TIOGA_DEVICE_TOL) u=2.0;
  if (v < -TIOGA_DEVICE_TOL || v > 1.0+TIOGA_DEVICE_TOL) v=2.0;
  if (w < -TIOGA_DEVICE_TOL || w > 1.0+TIOGA_DEVICE_TOL) w=2.0;
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
void d_interprbf(double *xcloud, double *P,double *weights,int np,int itype, int *iflag)
{
  int i,j,n;
  double M[NEQNS][NEQNS], B[NEQNS];
  double dd;
  int neqns=np+4;
  for(i=0;i<np;i++)
    {
      for(j=0;j<np;j++)
        {
          dd=0;
          for(n=0;n<3;n++) dd+=(xcloud[3*i+n]-xcloud[3*j+n])*(xcloud[3*i+n]-xcloud[3*j+n]);
          M[i][j]=std::exp(-dd);
        }
      M[i][np]=M[np][i]=1.0;
      for(j=np+1,n=0;j<np+4;j++,n++) M[i][j]=M[j][i]=xcloud[3*i+n];
      dd=0;
      for(n=0;n<3;n++) dd+=(xcloud[3*i+n]-P[n])*(xcloud[3*i+n]-P[n]);
      B[i]=std::exp(-dd);
    }
  for(i=np;i<np+4;i++)
    {
      for(j=np;j<np+4;j++)
        M[i][j]=0.0;
      B[i]=(i==np)?1.0:P[i-(np+1)];
    }
  d_solvec(M,B,iflag,neqns);
  for(i=0;i<np;i++) weights[i]=B[i];
}

TIOGA_GPU_DEVICE
void d_interpls1(double *xcloud, double *P, double *weights, int np, int *iflag)
{

 int i,j,ii,i3,m;
 double B[NEQNS][3],JB[NEQNS][3],S[3],X[3][3],Z[3],Y[NEQNS],dtest;
 double c,fac,M[NEQNS][3],w;
 m=0;
 for(i=0;i<np;i++)
   for(j=0;j<3;j++)
      M[i][j]=xcloud[m++];
 fac=1.0/np;
 for(j=0;j<3;j++)
   {
    S[j]=0;
    for(i=0;i<np;i++)
      {
	B[i][j]=M[i][j]-P[j];
	S[j]=S[j]+B[i][j];
      }
    S[j]=fac*S[j];
   }
 for(j=0;j<3;j++)
   for(i=0;i<np;i++)
     JB[i][j]=B[i][j]-S[j];

 for(i=0;i<3;i++)
   for(j=0;j<3;j++)
     for(ii=0;ii<np;ii++)
       X[i][j]+=B[ii][i]*JB[ii][j];
 Z[0]=S[0];
 Z[1]=S[1];
 Z[2]=S[2];

 invertMat3(X,Z);
 for(i=0;i<np;i++)
   {
     Y[i]=0;
     for(j=0;j<3;j++)
       Y[i]+=B[i][j]*Z[j];
   }
 c=0;
 for(j=0;j<3;j++)
   c+=S[j]*Z[j];
 for(i=0;i<np;i++)
   {
     Y[i]-=c;
     weights[i]=fac-Y[i];
   }
 return;

}

TIOGA_GPU_DEVICE
void d_computeNodalWeights(double *xv,double *xp,double *frac,int nvert)
{
  int i,j,k,isolflag;
  double f[24];
  double u,v,w;
  double oneminusU,oneminusV,oneminusW,oneminusUV;
  double lhs[3][3];
  double rhs[3];
  double frac2[8];
  //for(i=0;i<8;i++) frac[i]=0.0;
  //return; 
  switch(nvert)
    {
    case 4:
      //
      // tetrahedron
      //
      for(k=0;k<3;k++)
	{
	  for(j=0;j<3;j++)
	    lhs[j][k]=xv[k*3+j]-xv[3*3+j];
	  rhs[k]=xp[k]-xv[3*3+k];
	}
      //
      // invert the 3x3 matrix
      //
     invertMat3(lhs,rhs);
     isolflag=1;
     //d_solvec(lhs,rhs,&isolflag,3);
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
	  f[0*3+j]=xv[0*3+j]-xp[j];
	  f[1*3+j]=xv[1*3+j]-xv[0*3+j];
	  f[2*3+j]=xv[3*3+j]-xv[0*3+j];
	  f[3*3+j]=xv[4*3+j]-xv[0*3+j];
	  //
	  f[4*3+j]=xv[0*3+j]-xv[1*3+j]+xv[2*3+j]-xv[3*3+j];
	  f[5*3+j]=xv[0*3+j]-xv[3*3+j];
	  f[6*3+j]=xv[0*3+j]-xv[1*3+j];
	  f[7*3+j]=-xv[0*3+j]+xv[1*3+j]-xv[2*3+j]+xv[3*3+j];
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
	  f[0*3+j]=xv[0*3+j]-xp[j];
	  f[1*3+j]=xv[1*3+j]-xv[0*3+j];
	  f[2*3+j]=xv[2*3+j]-xv[0*3+j];
	  f[3*3+j]=xv[3*3+j]-xv[0*3+j];
	  //
	  f[4*3+j]=0;
	  f[5*3+j]=xv[0*3+j]-xv[2*3+j]-xv[3*3+j]+xv[5*3+j];
	  f[6*3+j]=xv[0*3+j]-xv[1*3+j]-xv[3*3+j]+xv[4*3+j];
	  f[7*3+j]=0.;
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
      double xcloud[24];
      int m;
      for(i=0;i<8;i++)
        for(j=0;j<3;j++)
         xcloud[m++]=xv[i*3+j];
      
      for(j=0;j<3;j++)
	{
	  f[0*3+j]=xv[0*3+j]-xp[j];
	  f[1*3+j]=xv[1*3+j]-xv[0*3+j];
	  f[2*3+j]=xv[3*3+j]-xv[0*3+j];
	  f[3*3+j]=xv[4*3+j]-xv[0*3+j];
	  //
	  f[4*3+j]=xv[0*3+j]-xv[1*3+j]+xv[2*3+j]-xv[3*3+j];
	  f[5*3+j]=xv[0*3+j]-xv[3*3+j]+xv[7*3+j]-xv[4*3+j];
	  f[6*3+j]=xv[0*3+j]-xv[1*3+j]+xv[5*3+j]-xv[4*3+j];
	  f[7*3+j]=-xv[0*3+j]+xv[1*3+j]-xv[2*3+j]+xv[3*3+j]+
	    xv[4*3+j]-xv[5*3+j]+xv[6*3+j]-xv[7*3+j];
	}
     
      //
      //d_interpls1(xcloud, xp, frac, 8 ,&isolflag);
      if (0) { 
      d_newtonSolve(f,&u,&v,&w);
      oneminusU=1-u;
      oneminusV=1-v;
      oneminusW=1-w;
      //
      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=oneminusU*oneminusV*w;
      frac[5]=u*oneminusV*w;
      frac[6]=u*v*w;
      frac[7]=oneminusU*v*w;     
      }
      //
      break;
    default:
      //printf("Interpolation not implemented for polyhedra with %d vertices\n",nvert);
      break;
    }
return;
}

TIOGA_GPU_DEVICE	   
void d_checkContainment(double *x, int *vconn,int *nc, int *nv, int ntypes, int *elementList,
			int *cellIndex, int adtElement, double *xsearch)
{
  int i,j,k,m,n,i3;
  int nvert;
  int icell,icell1;
  int passFlag;
  int isum;
  double xv[24];
  double frac[8];
  int ihigh=0;
  //
  icell=elementList[adtElement];
  if (ihigh==0) 
    {
      //
      // locate the type of the cell
      //
      /*
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
      */
      //
      // now collect all the vertices in the
      // array xv
      //
      i=icell;
      n=0;
      nvert=nv[n];
      nvert=8;
      for(m=0;m<nvert;m++)
	{
	  i3=3*(vconn[nvert*i+m]-TIOGA_DEVICE_BASE);
          //i3=0;
	  for(j=0;j<3;j++)
	    xv[3*m+j]=x[i3+j];
	}
      //
      d_computeNodalWeights(xv,xsearch,frac,nvert);

      //
      cellIndex[0]=icell;
      cellIndex[1]=0;
      //frac[0]=2;
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
              cellIndex[1]=0;
	    }
          // need to find another way to implement this
          //if (fabs(frac[m]) < TOL && cellRes[icell]==TIOGA_DEVICE_BIGVALUE) cellIndex[1]=1;
	}
    }
  else
    {
     // not implemented;
    }
 return;
}

TIOGA_GPU_DEVICE	   
void d_checkContainment2(double *x, int *vconn,int *nc, int *nv, int ntypes, int *elementList,
			 int *cellIndex, int adtElement, double *xsearch, 
			 int numverts[4][6], int faceInfo[4][24])
{
  int i,j,k,m,n,i3,ix;
  int nvert,nfaces;
  int icell,icell1;
  int passFlag;
  int isum;
  int e[8];
  double xface[12],xc[3],xi[3];
  int ihigh=0;
  //
  int ipoly;
  int intersect;
  //
  icell=elementList[adtElement];
  //
  // locate the type of the cell
  //
  /*
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
  */
  //
  // now collect all the vertices in the
  // array xv
  //
  i=icell;
  n=0;
  nvert=nv[n];
  xc[0]=xc[1]=xc[2]=0;
  for(m=0;m<nvert;m++)
    {
      e[m]=vconn[nvert*i+m];
      xc[0]+=x[3*(e[m]-TIOGA_DEVICE_BASE)];
      xc[1]+=x[3*(e[m]-TIOGA_DEVICE_BASE)+1];
      xc[2]+=x[3*(e[m]-TIOGA_DEVICE_BASE)+2];
    }
  //
  //  cell center
  //   
  xc[0]/=nvert;
  xc[1]/=nvert;
  xc[2]/=nvert;
  //
  cellIndex[0]=icell;
  cellIndex[1]=0;
  //
  if (nvert==4) { nfaces=4; ipoly=0;}
  if (nvert==5) { nfaces=5; ipoly=1;}
  if (nvert==6) { nfaces=5; ipoly=2;}
  if (nvert==8) { nfaces=6; ipoly=3;}
  //
  for(k=0;k<nfaces;k++)
    {
      for(j=0;j<numverts[ipoly][k];j++)
	{
	  ix=e[faceInfo[ipoly][4*k+j]-1]-TIOGA_DEVICE_BASE;
	  xface[3*j+0]=x[3*ix];
	  xface[3*j+1]=x[3*ix+1];
	  xface[3*j+2]=x[3*ix+2];
	}
      intersect=faceEdgeIntersectCheck(xface,xsearch,xc,xi,-TIOGA_DEVICE_TOL,TIOGA_DEVICE_TOL,
				       numverts[ipoly][k]);
      if (intersect) cellIndex[0]=-1;
    }
  
  return;
}

TIOGA_GPU_DEVICE
void d_searchIntersections_containment(int *donorId,
			               int *adtIntegers,
				       double *adtReals,
				       double *coord,
				       double *xsearch,
                                       double *x, int *nc, int *nv,int ntypes,
                                       int *elementList,
                                       int *vconn,
				       int nelem,
				       int ndim,
				       int *nchecks,
				       int numverts[4][6],
				       int faceInfo[4][24])
{
  int i,k,is;
  int d,nodeChild,flag;
  double element[6];
  double dv[3];
  double dtest,bdist;
  int nodeStack[maxStackSize];
  int mm=0;
  int nstack=1;
  double dmin[2];
  int cellIndex[2];
  int found;

  //typedef struct nodestack{
  //  int val;
  //  struct nodestack* next;
  //} nstack;  
  //
  int node=0;
  int level=0;

  nodeStack[0]=node;
  dmin[0]=dmin[1]=TIOGA_DEVICE_BIGVALUE;
  cellIndex[0]=cellIndex[1]=-1;
  *donorId=-1;
  found=0;
  while(nstack > 0 && !found) 
    {
      mm=0;
      for(is=0;is<nstack && !found ;is++)
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
            d_checkContainment2(x,vconn,nc,nv,ntypes,elementList,
				cellIndex,adtIntegers[4*node],xsearch,
				numverts,faceInfo);
            if (cellIndex[0] > -1 && cellIndex[1]==0) { *donorId=cellIndex[0]; found=1; }
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
		     nodeStack[mm+nstack]=nodeChild;
		     mm++;
		  }
	      }
	    }
	}
      //printf("mm=%d\n",mm);
      for(int j=0;j<mm;j++)
         nodeStack[j]=nodeStack[j+nstack];
      nstack=mm;
      level++;
    }
  return;
}


