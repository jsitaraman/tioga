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

extern "C"{
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}
			   
void MeshBlock::checkContainment(int *cellIndex, int adtElement, double *xsearch)
{
  int i,j,k,m,n,i3;
  int nvert;
  int icell,icell1;
  int passFlag;
  int isum;
  double xv[8][3];
  double frac[8];
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
	  i3=3*(vconn[n][nvert*i+m]-BASE);
	  for(j=0;j<3;j++)
	    xv[m][j]=x[i3+j];
	}
      //
      computeNodalWeights(xv,xsearch,frac,nvert);
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
	  if ((frac[m]+TOL)*(frac[m]-1.0-TOL) > 0) 
	    {
	      cellIndex[0]=-1;
	      return;
	    }
          if (fabs(frac[m]) < TOL && cellRes[icell]==BIGVALUE) cellIndex[1]=1;
	}
      
      return;
    }
  else
    {
      icell1=icell+BASE;
      cellIndex[0]=-1;
      cellIndex[1]=0;
      donor_inclusion_test(&icell1,xsearch,&passFlag,&(rst[ipoint]));
      if (passFlag) cellIndex[0]=icell;
      return;
    }

}
