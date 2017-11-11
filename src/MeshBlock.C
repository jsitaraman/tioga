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

extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  double computeCellVolume(double xv[8][3],int nvert);
  void deallocateLinkList(DONORLIST *temp);
  void deallocateLinkList2(INTEGERLIST *temp);  
}

using namespace tg_funcs;

void MeshBlock::setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
			int *wbcnodei,int *obcnodei,
			int ntypesi,int *nvi,int *nci,int **vconni)
{
  int i;
  //
  // set internal pointers
  //
  meshtag=btag;
  nnodes=nnodesi;
  x=xyzi;
  iblank=ibli;
  nwbc=nwbci;
  nobc=nobci;
  wbcnode=wbcnodei;
  obcnode=obcnodei;
  //
  ntypes=ntypesi;
  //
  nv=nvi;
  nc=nci;
  vconn=vconni;
  //
  //tracei(nnodes);
  //for (int i = 0; i < ntypes; i++) tracei(nc[i]);
  ncells=0;
  for (int i = 0; i < ntypes; i++) ncells+=nc[i];
}

void MeshBlock::setFaceData(int _gtype, int _nftype, int *_nf, int *_nfv,
    int **_f2v, int *_f2c, int *_c2f, int *_ib_face, int nOver, int nWall, int nMpi,
    int *oFaces, int *wFaces, int *mFaces, int *procR, int *idR)
{
  gridType = _gtype;
  nftype = _nftype;
  nf = _nf;
  nfv = _nfv;
  fconn = _f2v;
  f2c = _f2c;
  c2f = _c2f;
  iblank_face = _ib_face;
  nOverFaces = nOver;
  nWallFaces = nWall;
  nMpiFaces = nMpi;
  overFaces = oFaces;
  wallFaces = wFaces;
  mpiFaces = mFaces;
  mpiProcR = procR;
  mpiFidR = idR;

  nfaces = 0;
  for (int i = 0; i < nftype; i++)
    nfaces += nf[i];
}

#ifdef _GPU
void MeshBlock::setDeviceData(double* xyz, double* coords, int* ibc, int* ibf)
{
  mb_d.setDeviceData(xyz,coords,ibc,ibf);
}

#endif

void MeshBlock::preprocess(void)
{
  // set all iblanks = 1
  for (int i = 0; i < nnodes; i++) iblank[i] = NORMAL;

  // find oriented bounding boxes
  free(obb);

  obb = (OBB *) malloc(sizeof(OBB));

  findOBB(x,obb->xc,obb->dxc,obb->vec,nnodes);

  tagBoundary();
}

void MeshBlock::updateOBB(void)
{
  free(obb);
  obb = (OBB *) malloc(sizeof(OBB));

  findOBB(x,obb->xc,obb->dxc,obb->vec,nnodes);
}

/** Calculate 'cellRes' / 'nodeRes' (cell volume) for each cell / node
 *  and tag overset-boundary nodes by setting their nodeRes to a big value
 */
void MeshBlock::tagBoundary(void)
{
  std::vector<int> inode;
  double xv[8][3];
  std::vector<double> xv2;
  std::vector<int> iflag(nnodes, 0);

  // Do this only once
  // i.e. when the meshblock is first initialized, cellRes would be NULL in
  // this case
  cellRes.resize(ncells);
  nodeRes.resize(nnodes);

  if (userSpecifiedNodeRes == NULL && userSpecifiedCellRes == NULL)
  {
    for (int i = 0; i < nnodes; i++) nodeRes[i] = 0.0;

    int k = 0;
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      inode.resize(nvert);
      if (nvert > 8) xv2.resize(nvert*3);

      for (int i = 0; i < nc[n]; i++)
      {
        double vol = 0.;

        if (nvert > 8)
        {
          for (int m = 0; m < nvert; m++)
          {
            inode[m] = vconn[n][nvert*i+m]-BASE;
            int i3 = 3*inode[m];
            for (int j = 0; j < 3; j++)
              xv2[m*3+j] = x[i3+j];
          }
          vol = computeVolume(xv2.data(), nvert, 3);
        }
        else
        {
          for (int m = 0; m < nvert; m++)
          {
            inode[m] = vconn[n][nvert*i+m]-BASE;
            int i3 = 3*inode[m];
            for (int j = 0; j < 3; j++)
              xv[m][j] = x[i3+j];
          }
          vol = computeVolume(&xv[0][0], nvert, 3);
        }

        cellRes[k++] = vol*resolutionScale;
        for (int m = 0; m < nvert; m++)
        {
          iflag[inode[m]]++;
          nodeRes[inode[m]] += vol*resolutionScale;
        }
      }
    }
  }
  else
  {
    int k = 0;
    for (int n = 0; n < ntypes; n++)
    {
      for (int i = 0; i < nc[n]; i++)
      {
        cellRes[k] = userSpecifiedCellRes[k];
        k++;
      }
    }
    for (int k = 0; k < nnodes; k++) nodeRes[k] = userSpecifiedNodeRes[k];
  }

  // compute nodal resolution as the average of all the cells associated with it.
  // This takes care of partition boundaries as well.
  for (int i = 0; i < nnodes; i++)
  {
    if (iflag[i] != 0)
      nodeRes[i] /= iflag[i];
    iflag[i] = 0;
  }

  // now tag the boundary nodes
  // reuse the iflag array
  if (iartbnd)
  {
    for (int i = 0; i < nobc; i++)
      nodeRes[obcnode[i]-BASE] = BIGVALUE;
  }
  else
  {
    for (int i = 0; i < nobc; i++)
    {
      iflag[(obcnode[i]-BASE)] = 1;
    }

    // now tag all the nodes of boundary cells
    // to be mandatory receptors
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      inode.resize(nvert);
      for (int i = 0; i < nc[n]; i++)
      {
        int itag = 0;
        for (int m = 0; m < nvert; m++)
        {
          inode[m] = vconn[n][nvert*i+m]-BASE;
          if (iflag[inode[m]]) itag = 1;
        }
        if (itag)
        {
          for (int m = 0; m < nvert; m++)
          {
            nodeRes[inode[m]]=BIGVALUE;
          }
        }
      }
    }

    // now tag all the cells which have mandatory receptors as nodes as not
    // acceptable donors
    int k = 0;
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      inode.resize(nvert);
      for (int i = 0; i < nc[n]; i++)
      {
        for (int m = 0; m < nvert; m++)
        {
          inode[m] = vconn[n][nvert*i+m]-BASE;
          if (nodeRes[inode[m]] == BIGVALUE) //iflag[inode[m]])
          {
            cellRes[k]=BIGVALUE;
            break;
          }
        }
        k++;
      }
    }
  }
}

void MeshBlock::setupADT(void)
{
  PUSH_NVTX_RANGE("BuildADT", 2);
  for (int d = 0; d < 3; d++)
  {
    aabb[d]   =  BIG_DOUBLE;
    aabb[d+3] = -BIG_DOUBLE;
  }

  /// TODO: take in a bounding box of the region we're interested in searching (search.c)
  elementBbox.resize(ncells*6);
  elementList.resize(ncells);

  double xmin[3], xmax[3];
  for (int i = 0; i < ncells; i++)
  {
    int isum = 0;
    int ic = i;
    int n;
    for (n = 0; n < ntypes; n++)
    {
      isum += nc[n];
      if (ic < isum)
      {
        ic = i - (isum - nc[n]);
        break;
      }
    }

    int nvert = nv[n];
    xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
    xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
    for (int m = 0; m < nvert; m++)
    {
      int i3 = 3*(vconn[n][nvert*ic+m]-BASE);
      for (int j = 0; j < 3; j++)
      {
        xmin[j] = min(xmin[j], x[i3+j]);
        xmax[j] = max(xmax[j], x[i3+j]);
        // Overall partition bounding box
        aabb[j]   = min(aabb[j], x[i3+j]);
        aabb[j+3] = max(aabb[j+3], x[i3+j]);
      }
    }

    elementBbox[6*i] = xmin[0];
    elementBbox[6*i+1] = xmin[1];
    elementBbox[6*i+2] = xmin[2];
    elementBbox[6*i+3] = xmax[0];
    elementBbox[6*i+4] = xmax[1];
    elementBbox[6*i+5] = xmax[2];

    elementList[i] = i;
  }

  if (adt)
  {
    adt->clearData();
  }
  else
  {
    adt=new ADT[1];
  }

  adt->buildADT(2*nDims, ncells, elementBbox.data());
  POP_NVTX_RANGE;

#ifdef _GPU
  ncells_adt = ncells;
  mb_d.dataToDevice(nDims,nnodes,ncells,ncells_adt,nsearch,nv,nc,elementList.data(),
                    elementBbox.data(),isearch.data(),xsearch.data(),myid);

  adt_d.copyADT(adt);

  /* ---- Direct Cut Setup ---- */

  mb_d.extraDataToDevice(vconn[0]);
#endif
}


void MeshBlock::rebuildADT(void)
{
  /// TODO: make sure we're including all possible future donors at boundaries where other grids may move in
  /// [Include any elements assigned to 'DC_CUT' from hole cutting]?
  PUSH_NVTX_RANGE("ReBuildADT", 1);

  std::set<int> donorEles, adtEles;

  if (haveDonors)
  {
    // Collect all current donor elements & their nearest neighbors
    for (int i = 0; i < donorId.size(); i++)
    {
      int ic = donorId[i];
      donorEles.insert(ic);
    }

    for (auto ic : donorEles)
    {
      if (ic < 0) continue;

      int N = 2*nDims;
      adtEles.insert(ic);
      for (int j = 0; j < N; j++)
      {
        adtEles.insert(c2c[N*ic+j]);
      }
    }

    // Also add in elements from boundaries [incl. overset/cut boundaries]
    for (int i = 0; i < nreceptorFaces; i++)
    {
      int ff = ftag[i];
      int ic1 = f2c[2*ff+0];
      int ic2 = f2c[2*ff+1];
      adtEles.insert(ic1);
      adtEles.insert(ic2);
    }

    for (int i = 0; i < nMpiFaces; i++)
    {
      int ff = mpiFaces[i];
      int ic = f2c[2*ff+0];
      adtEles.insert(ic);
    }
  }
  else
  {
    // We're probably still in initialization - just add all elements [Need to search full grid at least once]
    for (int i = 0; i < ncells; i++)
      adtEles.insert(i);
  }

  adtEles.erase(-1);

  // Remove eles which don't overlap with the search point oriented bounding box

  OBB obq;
  findOBB(xsearch.data(),obq.xc,obq.dxc,obq.vec,nsearch);

  std::set<int> newAdtEles; // = adtEles;

  for (auto ic : adtEles)
  {
    int nvert = nv[0];

    double xmin[3];
    double xmax[3];
    double xd[3], dxc[3];
    xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
    xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
    for (int m = 0; m < nvert; m++)
    {
      int i3 = 3*(vconn[0][nvert*ic+m]-BASE);
      for (int j = 0; j < 3; j++)
      {
        xd[j] = 0;
        for (int k = 0; k < 3; k++)
          xd[j] += (x[i3+k]-obq.xc[k])*obq.vec[j][k];
        xmin[j] = min(xmin[j],xd[j]);
        xmax[j] = max(xmax[j],xd[j]);
      }
      for (int j = 0; j < 3; j++)
      {
        xd[j] = (xmax[j]+xmin[j])*0.5;
        dxc[j] = (xmax[j]-xmin[j])*0.5;
      }
    }

    if (fabs(xd[0]) <= (dxc[0]+obq.dxc[0]) &&
        fabs(xd[1]) <= (dxc[1]+obq.dxc[1]) &&
        fabs(xd[2]) <= (dxc[2]+obq.dxc[2]))
    {
      newAdtEles.insert(ic);
    }
  }

  adtEles = newAdtEles;
  ncells_adt = adtEles.size();
//printf("Rank %d - ncells_adt %d\n",myid,ncells_adt); /// DEBUGGING
  elementList.resize(ncells_adt);
  elementBbox.resize(ncells_adt*6);

  int I = 0;
  for (auto ic : adtEles)
  {
    elementList[I] = ic;
    I++;
  }

  double xmin[3], xmax[3];
  for (int i = 0; i < ncells_adt; i++)
  {
    int isum = 0;
    int ic = elementList[i];
    int n;
    for (n = 0; n < ntypes; n++)
    {
      isum += nc[n];
      if (ic < isum)
      {
        ic = ic - (isum - nc[n]);
        break;
      }
    }

    int nvert = nv[n];
    xmin[0] = xmin[1] = xmin[2] =  BIGVALUE;
    xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
    for (int m = 0; m < nvert; m++)
    {
      int i3 = 3*(vconn[n][nvert*ic+m]-BASE);
      for (int j = 0; j < 3; j++)
      {
        xmin[j] = min(xmin[j], x[i3+j]);
        xmax[j] = max(xmax[j], x[i3+j]);
      }
    }

    elementBbox[6*i] = xmin[0];
    elementBbox[6*i+1] = xmin[1];
    elementBbox[6*i+2] = xmin[2];
    elementBbox[6*i+3] = xmax[0];
    elementBbox[6*i+4] = xmax[1];
    elementBbox[6*i+5] = xmax[2];
  }

  if (adt)
  {
    adt->clearData();
  }
  else
  {
    adt=new ADT[1];
  }

  adt->buildADT(2*nDims, ncells_adt, elementBbox.data());

#ifdef _GPU
  mb_d.updateADTData(ncells_adt,elementList.data(),elementBbox.data());

  adt_d.copyADT(adt);
#endif
  POP_NVTX_RANGE;
}

void MeshBlock::writeGridFile(int bid)
{
  char fname[80];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"part%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
	  ncells);
  for (int i = 0; i < nnodes; i++)
    {
      fprintf(fp,"%.14e %.14e %.14e %d\n",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
    }

  ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
    {
      nvert=nv[n];
      for (int i = 0; i < nc[n]; i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
	  else if (nvert==8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}

void MeshBlock::writeCellFile(int bid, int* flag)
{
  char fname[80];
  char qstr[2];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"cell%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\",\"IBLANK_CELL\" ");
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,
	  ncells);
  fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=NODAL,5=CELLCENTERED)\n");
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i+1]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i+2]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%d\n",iblank[i]);
  if (flag != NULL)
    for (int i = 0; i < ncells; i++) fprintf(fp,"%d\n",flag[i]);
  else
    for (int i = 0; i < ncells; i++) fprintf(fp,"%d\n",iblank_cell[i]);
  ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
    {
      nvert=nv[n];
      for (int i = 0; i < nc[n]; i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
    else if (nvert>=8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}

void MeshBlock::writeFlowFile(int bid,double *q,int nvar,int type)
{
  char fname[80];
  char qstr[2];
  char intstring[7];
  char hash,c;
  int i,n,j;
  int bodytag;
  FILE *fp;
  int ba;
  int nvert;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"flow%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\" ");
  for (int i = 0; i < nvar; i++)
    {
      sprintf(qstr,"Q%d",i);
      fprintf(fp,"\"%s\",",qstr);
    }
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
	  ncells);

  if (type==0)
    {
      for (int i = 0; i < nnodes; i++)
	{
	  fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
	  for (int j = 0; j < nvar; j++)
	    fprintf(fp,"%lf ",q[i*nvar+j]);      
	  //for (int j = 0; j < nvar; j++)
	  //  fprintf(fp,"%lf ", x[3*i]+x[3*i+1]+x[3*i+2]);
          fprintf(fp,"\n");
	}
    }
  else
    {
      for (int i = 0; i < nnodes; i++)
        {
          fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
          for (int j = 0; j < nvar; j++)
            fprintf(fp,"%lf ",q[j*nnodes+i]);
          fprintf(fp,"\n");
        }
    }
  ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
    {
      nvert=nv[n];
      for (int i = 0; i < nc[n]; i++)
	{
	  if (nvert==4)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+3]+ba);
	    }
	  else if (nvert==5) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+4]+ba);
	    }
	  else if (nvert==6) 
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+5]+ba);
	    }
	  else if (nvert==8)
	    {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      vconn[n][nvert*i]+ba,
		      vconn[n][nvert*i+1]+ba,
		      vconn[n][nvert*i+2]+ba,
		      vconn[n][nvert*i+3]+ba,
		      vconn[n][nvert*i+4]+ba,
		      vconn[n][nvert*i+5]+ba,
		      vconn[n][nvert*i+6]+ba,
		      vconn[n][nvert*i+7]+ba);
	    }
	}
    }
  fclose(fp);
  return;
}
  
void MeshBlock::getWallBounds(int *mtag,int *existWall, double wbox[6])
{
  *mtag = meshtag - BASE; // + 1?
  if (nwbc <= 0)
  {
    *existWall = 0;
    for (int i = 0; i < 3; i++)
    {
      wbox[i]   = BIGVALUE;
      wbox[3+i] = -BIGVALUE;
    }
    return;
  }

  *existWall=1;
  wbox[0] = wbox[1] = wbox[2] =  BIGVALUE;
  wbox[3] = wbox[4] = wbox[5] = -BIGVALUE;

  for (int i = 0; i < nwbc; i++)
  {
    int inode = wbcnode[i]-BASE;
    int i3 = 3*inode;
    for (int j = 0; j < 3; j++)
    {
      wbox[j] = min(wbox[j],x[i3+j]);
      wbox[j+3] = max(wbox[j+3],x[i3+j]);
    }
  }
}

void MeshBlock::getOversetBounds(int *mtag, int *existOver, double obox[6])
{
  *mtag = meshtag - BASE;
  if (nobc <= 0)
  {
    *existOver=0;
    for (int i = 0; i < 3; i++)
    {
      obox[i]   =  BIGVALUE;
      obox[3+i] = -BIGVALUE;
    }
    return;
  }

  *existOver = 1;
  obox[0] = obox[1] = obox[2] =  BIGVALUE;
  obox[3] = obox[4] = obox[5] = -BIGVALUE;

  for (int i = 0; i < nobc; i++)
  {
    int ind = 3*obcnode[i]-BASE;
    for (int j = 0; j < 3; j++)
    {
      obox[j]   = min(obox[j],   x[ind+j]);
      obox[j+3] = max(obox[j+3], x[ind+j]);
    }
  }
}
  
void MeshBlock::markWallBoundary(int *sam,int nx[3],double extents[6])
{
  std::vector<int> iflag(ncells);
  std::vector<int> inode(nnodes);

  // Mark all wall boundary nodes
  for (int i = 0; i < nwbc; i++)
  {
    int ii = wbcnode[i]-BASE;
    inode[ii] = 1;
  }

  // Mark all wall boundary cells (Any cell with wall-boundary node)
  int m = 0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      for (int j = 0; j < nvert; j++)
      {
        int ii = vconn[n][nvert*i+j]-BASE;
        if (inode[ii] == 1)
        {
          iflag[m] = 1;
          break;
        }
      }
      m++;
    }
  }

  // find delta's in each directions
  double ds[3];
  for (int k = 0; k < 3; k++) ds[k] = (extents[k+3]-extents[k])/nx[k];

  // mark sam cells with wall boundary cells now
  int imin[3];
  int imax[3];
  m = 0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      if (iflag[m] == 1)
      {
        // find the index bounds of each wall boundary cell bounding box
        imin[0] = imin[1] = imin[2] =  BIGINT;
        imax[0] = imax[1] = imax[2] = -BIGINT;
        for (int j = 0; j < nvert; j++)
        {
          int i3 = 3*(vconn[n][nvert*i+j]-BASE);
          for (int k = 0; k < 3; k++)
          {
            double xv = x[i3+k];
            int iv = floor((xv-extents[k])/ds[k]);
            imin[k] = min(imin[k],iv);
            imax[k] = max(imax[k],iv);
          }
        }

        for (int j = 0; j < 3; j++)
        {
          imin[j] = max(imin[j],0);
          imax[j] = min(imax[j],nx[j]-1);
        }

        // mark sam to 2
        for (int kk = imin[2]; kk < imax[2]+1; kk++)
        {
          for (int jj = imin[1]; jj < imax[1]+1; jj++)
          {
            for (int ii = imin[0]; ii < imax[0]+1; ii++)
            {
              int mm = (kk*nx[1] + jj)*nx[0] + ii;
              sam[mm] = 2;
            }
          }
        }
      }
      m++;
    }
  }
}

void MeshBlock::markOversetBoundary(int *sam,int nx[3],double extents[6])
{
  std::vector<int> iflag(ncells);
  std::vector<int> inode(nnodes);

  // Mark all wall boundary nodes
  for (int i = 0; i < nobc; i++) /// TESTING
  {
    int ii = obcnode[i]-BASE; /// TESTING
    inode[ii] = 1;
  }

  // Mark all wall boundary cells (Any cell with wall-boundary node)
  int m = 0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      for (int j = 0; j < nvert; j++)
      {
        int ii = vconn[n][nvert*i+j]-BASE;
        if (inode[ii] == 1)
        {
          iflag[m] = 1;
          break;
        }
      }
      m++;
    }
  }

  // find delta's in each directions
  double ds[3];
  for (int k = 0; k < 3; k++) ds[k] = (extents[k+3]-extents[k])/nx[k];

  // mark sam cells with wall boundary cells now
  int imin[3];
  int imax[3];
  m = 0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      if (iflag[m] == 1)
      {
        // find the index bounds of each wall boundary cell bounding box
        imin[0] = imin[1] = imin[2] =  BIGINT;
        imax[0] = imax[1] = imax[2] = -BIGINT;
        for (int j = 0; j < nvert; j++)
        {
          int i3 = 3*(vconn[n][nvert*i+j]-BASE);
          for (int k = 0; k < 3; k++)
          {
            double xv = x[i3+k];
            int iv = floor((xv-extents[k])/ds[k]);
            imin[k] = min(imin[k],iv);
            imax[k] = max(imax[k],iv);
          }
        }

        for (int j = 0; j < 3; j++)
        {
          imin[j] = max(imin[j],0);
          imax[j] = min(imax[j],nx[j]-1);
        }

        // mark sam to 2
        for (int kk = imin[2]; kk < imax[2]+1; kk++)
        {
          for (int jj = imin[1]; jj < imax[1]+1; jj++)
          {
            for (int ii = imin[0]; ii < imax[0]+1; ii++)
            {
              int mm = (kk*nx[1] + jj)*nx[0] + ii;
              sam[mm] = 2;
            }
          }
        }
      }
      m++;
    }
  }
}

void MeshBlock::getReducedOBB(OBB *obc,double *realData) 
{
  double bbox[6],xd[3];

  for (int j = 0; j < 3; j++)
  {
    realData[j]   = BIGVALUE;
    realData[j+3] = -BIGVALUE;
  }

  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      bbox[0] = bbox[1] = bbox[2] =  BIGVALUE;
      bbox[3] = bbox[4] = bbox[5] = -BIGVALUE;

      for (int m = 0; m < min(nvert,8); m++) /// HACK-ish...
      {
        int i3 = 3*(vconn[n][nvert*i+m]-BASE);
        for (int j = 0; j < 3; j++)
          xd[j] = 0;

        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            xd[j] += (x[i3+k]-obc->xc[k])*obc->vec[j][k];

        for (int j = 0; j < 3; j++)
        {
          bbox[j]   = min(bbox[j],xd[j]);
          bbox[j+3] = max(bbox[j+3],xd[j]);
        }
      }

      int iflag = 0;
      for (int j = 0; j < 3; j++)
      {
        iflag = (iflag || (bbox[j] > obc->dxc[j]));
        iflag = (iflag || (bbox[j+3] < -obc->dxc[j]));
      }

      if (iflag) continue;

      for (int m = 0; m < min(nvert,8); m++) /// HACK-ish...
      {
        int i3 = 3*(vconn[n][nvert*i+m]-BASE);
        for (int j = 0; j < 3; j++)
          xd[j] = 0;

        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            xd[j] += (x[i3+k]-obb->xc[k])*obb->vec[j][k];

        for (int j = 0; j < 3; j++)
        {
          realData[j]   = min(realData[j],xd[j]);
          realData[j+3] = max(realData[j+3],xd[j]);
        }
      }
    }
  }

  for (int j = 0; j < 6; j++)
    bbox[j] = realData[j];

  for (int j = 0; j < 3; j++)
  {
    realData[j] = obb->xc[j];

    for (int k = 0; k < 3; k++)
      realData[j] += ((bbox[k]+bbox[k+3])*0.5)*obb->vec[k][j];

    realData[j+3] = (bbox[j+3]-bbox[j])*0.51;
  }
}
	      
void MeshBlock::getQueryPoints(OBB *obc,
			       int *nints,int **intData,
			       int *nreals, double **realData)
{
  int i,j,k;
  int i3;
  double xd[3];
  int *inode;
  int iptr;
  int m;

  inode=(int *)malloc(sizeof(int)*nnodes);
  *nints=*nreals=0; 
  for (int i = 0; i < nnodes; i++)
    {
      i3=3*i;
      for (int j = 0; j < 3; j++) xd[j]=0;
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
	  xd[j]+=(x[i3+k]-obc->xc[k])*obc->vec[j][k];

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
  for (int i = 0; i < *nints; i++)
    {
      i3=3*inode[i];
      (*intData)[i]=inode[i];
      (*realData)[m++]=x[i3];
      (*realData)[m++]=x[i3+1];
      (*realData)[m++]=x[i3+2];
    }
  //
  free(inode);
}  
  
void MeshBlock::writeOBB(int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"box%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,
	  1);

  for(l=0;l<2;l++)
    {
      il=2*(l%2)-1;
      for (int k = 0; k < 2; k++)
	{
	  ik=2*(k%2)-1;
	  for (int j = 0; j < 2; j++)
	    {
	      ij=2*(j%2)-1;
	      xx[0]=xx[1]=xx[2]=0;
	      for (int m = 0; m < 3; m++)
		xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
		  +ik*obb->vec[1][m]*obb->dxc[1]
		  +il*obb->vec[2][m]*obb->dxc[2];	      
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	    }
	}
    }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fprintf(fp,"%e %e %e\n",obb->xc[0],obb->xc[1],obb->xc[2]);
  for (int k = 0; k < 3; k++)
   fprintf(fp,"%e %e %e\n",obb->vec[0][k],obb->vec[1][k],obb->vec[2][k]);
  fprintf(fp,"%e %e %e\n",obb->dxc[0],obb->dxc[1],obb->dxc[2]);
  fclose(fp);
}

int MeshBlock::findPointDonor(double *x_pt)
{
  int foundCell;
  double rst[3] = {0.0};
  adt->searchADT(this,&foundCell,x_pt,rst);
  return foundCell;
}

std::unordered_set<int> MeshBlock::findCellDonors(double *bbox)
{
  std::unordered_set<int> foundCells;
  adt->searchADT_box(elementList.data(),foundCells,bbox);
  return foundCells;
}

//
// destructor that deallocates all the
// the dynamic objects inside
//
MeshBlock::~MeshBlock()
{
  //
  // free all data that is owned by this MeshBlock
  // i.e not the pointers of the external code.
  //
  if (adt) delete[] adt;
  if (donorList) {
    for (int i = 0; i < nnodes; i++) deallocateLinkList(donorList[i]);
    free(donorList);
  }
  return;  /// Why is this here??
  if (interpListCart) delete [] interpListCart;
  if (obb) free(obb);
  if (interp2donor) free(interp2donor);
  if (cancelList) deallocateLinkList2(cancelList);
  if (ctag) free(ctag);
  if (ftag) free(ftag);
  if (!iartbnd) free(iblank_cell);
  if (pointsPerCell) free(pointsPerCell);
  if (pointsPerFace) free(pointsPerFace);
  if (picked) free(picked);
  if (rxyzCart) free(rxyzCart);
  if (donorIdCart) free(donorIdCart);
  if (pickedCart) free(pickedCart);
  if (ctag_cart) free(ctag_cart);

  // need to add code here for other objects as and
  // when they become part of MeshBlock object  
}

//
// set user specified node and cell resolutions
//
void MeshBlock::setResolutions(double *nres,double *cres)
{
  userSpecifiedNodeRes=nres;
  userSpecifiedCellRes=cres;
}

void MeshBlock::setTransform(double* mat, double* off, int ndim)
{
  if (ndim != nDims)
    ThrowException("MeshBlock::set_transform: input ndim != nDims");

  rrot = true;
  Rmat = mat;
  offset = off;

  if (adt)
    adt->setTransform(mat,off,ndim);
}
