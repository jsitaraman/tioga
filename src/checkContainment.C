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

extern "C"{
  void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert);
}

using namespace tg_funcs;
			   
void MeshBlock::checkContainment(int *cellIndex, int adtElement, double *xsearch, double *rst)
{
  int icell = elementList[adtElement];

  if (ihigh == 0)
  {
    //
    // locate the type of the cell
    //
    int isum = 0;
    int ic = -1;
    int n = 0;
    for (n = 0; n < ntypes; n++)
    {
      isum += nc[n];
      if (icell < isum)
      {
        ic = icell-(isum-nc[n]);
        break;
      }
    }

    if (ic < 0)
    {
      printf("TIOGA: invalid icell (adtElement) in checkContainment\n");
      exit(1);
    }

    int nvert = nv[n];
    if (nvert > 8)
      nvert = nNodesToFirstOrder(ncf[n], nv[n]);
//    if (nvert > 8)
//    {
//      /* --- Quadratic or higher-order shape functions - use general func --- */

//      std::vector<double> xv2(nvert*3);
//      for (int m = 0; m < nvert; m++)
//      {
//        int i3 = 3*(vconn[n][nvert*ic+m]-BASE);
//        for (int j = 0; j < 3; j++)
//          xv2[m*3+j] = x[i3+j];
//      }

//      double refloc[3];
//      bool isInEle = getRefLocNewton(xv2.data(), xsearch, &refloc[0], nvert, 3);

//      if (!isInEle)
//        *cellIndex = -1;
//      else
//        *cellIndex = icell;

//      return;
//    }
//    else
//    {
      /* --- Linear shape functions --- */

      double xv[8][3];
      double frac[8];

      for (int m = 0; m < nvert; m++)
      {
        int i3 = 3*(vconn[n][nv[n]*ic+m]-BASE);
        for (int j = 0; j < 3; j++)
          xv[m][j] = x[i3+j];
      }

      computeNodalWeights(xv,xsearch,frac,nvert);

      // if any of the nodal weights are not in between [-TOL 1+TOL] discard cell
      for (int m = 0; m < nvert; m++)
      {
        if ((frac[m]+TOL)*(frac[m]-1.0-TOL) > 0)
        {
          *cellIndex=-1;
          return;
        }
      }
      return;
//    }
  }
  else
  {
    int icell1 = icell+BASE;
    *cellIndex = -1;
    int passFlag;
    donor_inclusion_test(&icell1,xsearch,&passFlag,rst);
    if (passFlag) *cellIndex = icell;
  }

}
