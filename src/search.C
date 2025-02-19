// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"
#include "MeshBlock.h"
#include <unordered_map>
#include <iostream>

extern "C" {
  void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes);
  void writebbox(OBB *obb,int bid);
  void writePoints(double *x,int nsearch,int bid);
  void uniquenodes(double *x,int *meshtag,double *rtag,int *itag,int *nn);
  void uniquenodes_octree(double *x,int *meshtag,double *rtag,int *itag,int *nn);
}

namespace {

/** Determine the unique nodes by a global identifier
 *
 *  The function will create a mapping such that all duplicate nodes will point
 *  to the original node (as determined by a global identifier) in the `itag`
 *  array. It will update the nodal resolutions the shared nodes such that the
 *  resolutions upon exit will be the maximum resolution amongst all the
 *  duplicate nodes.
 *
 *  \param[in] node_ids Global IDs for the nodes across all MPI ranks
 *  \param[inout] node_res The nodal resolutions
 *  \param[out] itag The local index of the original node (duplicate to original mapping)
 *  \param[in] nnodes The size of the arrays
 */
void uniquenode_map(uint64_t* node_ids, double* node_res, int* itag, int nnodes)
{
    std::unordered_map<uint64_t, int> lookup;

    for (int i=0; i < nnodes; i++) {
        auto found = lookup.find(node_ids[i]);
        if (found != lookup.end()) {
            // This is a duplicate node, store the index to the original node
            // found previously
            itag[i] = found->second;

            // Update the original node's resolution to be the max of either
            // node resolution
            node_res[found->second] = std::max(node_res[found->second], node_res[i]);
        } else {
            // This is the first appearance of the unique ID, stash it in the
            // lookup table
            lookup[node_ids[i]] = i;
            itag[i] = i;
        }
    }

    // The max node resolution was stored off in the original node, propagate
    // this to all the duplicates
    for (int i=0; i < nnodes; i++)
        node_res[i] = node_res[itag[i]];
}
}


void MeshBlock::search(void)
{
  int i,j,k,l,m,n,p,i3;
  int ndim;
  int iptr,isum,nvert;
  OBB *obq;
  int *icell;
  int *itag;
  int cell_count; 
  int cellindex;
  double xd[3];
  double dxc[3];
  double xmin[3];
  double xmax[3];
  int *dId;
  //
  // form the bounding box of the 
  // query points
  //
  if (nsearch == 0) {
    donorCount=0;
    return;
  }
 
  if (uniform_hex) {
    search_uniform_hex();
    return;
  }

  obq=(OBB *) malloc(sizeof(OBB));
  
findOBB(xsearch,obq->xc,obq->dxc,obq->vec,nsearch);


  //writebbox(obq,4);
  //writePoints(xsearch,nsearch,4);
  //
  // find all the cells that may have intersections with
  // the OBB
  //
  icell=(int *)malloc(sizeof(int)*ncells);
  for(i=0;i<ncells;i++) icell[i]=-1;
  iptr=-1;
  cell_count=0;
  p=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  //
	  // find each cell that has
	  // overlap with the bounding box
	  //
	  xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
	  xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
	  for(m=0;m<nvert;m++)
	    {
	      i3=3*(vconn[n][nvert*i+m]-BASE);	      
	      for(j=0;j<3;j++)
		{
		  xd[j]=0;
		  for(k=0;k<3;k++)
		    xd[j]+=(x[i3+k]-obq->xc[k])*obq->vec[j][k];
		  xmin[j]=TIOGA_MIN(xmin[j],xd[j]);
		  xmax[j]=TIOGA_MAX(xmax[j],xd[j]);
		}
	      for(j=0;j<3;j++)
		{
		  xd[j]=(xmax[j]+xmin[j])*0.5;
		  dxc[j]=(xmax[j]-xmin[j])*0.5;
		}
	    }
	  if (fabs(xd[0]) <= (dxc[0]+obq->dxc[0]) &&
	      fabs(xd[1]) <= (dxc[1]+obq->dxc[1]) &&
	      fabs(xd[2]) <= (dxc[2]+obq->dxc[2])) 
	    {
	      //
	      // create a LIFO stack
	      // with all the cells that 
	      // have bounding box intersection with
	      // the QP bounding box
	      //
	      icell[p]=iptr;
	      iptr=p;
	      cell_count++;
	    }
	  p++;
	}
    }
  //
  // now find the axis aligned bounding box
  // of each cell in the LIFO stack to build the
  // ADT
  //

  if (elementBbox) TIOGA_FREE(elementBbox);
  if (elementList) TIOGA_FREE(elementList);
  elementBbox=(double *)malloc(sizeof(double)*cell_count*6);
  elementList=(int *)malloc(sizeof(int)*cell_count);
  //
  k=iptr;
  l=0;
  p=0;
  //for(k=0;k<ncells;k++)
  while(k!=-1)
    {
      cellindex=k;
      isum=0;
      for(n=0;n<ntypes;n++) 
	{
	  isum+=nc[n];
	  if (cellindex < isum)
	    {
	      i=cellindex-(isum-nc[n]);
	      break;
	    }
	}
      nvert=nv[n];
      xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
      xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
      for(m=0;m<nvert;m++)
	{
	  i3=3*(vconn[n][nvert*i+m]-BASE);
	  for(j=0;j<3;j++)
	    {
	      xmin[j]=TIOGA_MIN(xmin[j],x[i3+j]);
	      xmax[j]=TIOGA_MAX(xmax[j],x[i3+j]);
	    }
	}
      //
      elementBbox[l++]=xmin[0];
      elementBbox[l++]=xmin[1];
      elementBbox[l++]=xmin[2];
      elementBbox[l++]=xmax[0];
      elementBbox[l++]=xmax[1];
      elementBbox[l++]=xmax[2];
      //
      elementList[p++]=k;
      //
      k=icell[k];
    }
  //
  // build the ADT now
  //
  if (adt) 
   {
    adt->clearData();
   }
  else
   {
    adt=new ADT[1];
   }
  ndim=6;
  //
  adt->buildADT(ndim,cell_count,elementBbox);
  //
  if (donorId) TIOGA_FREE(donorId);
  donorId=(int*)malloc(sizeof(int)*nsearch);
  if (xtag) TIOGA_FREE(xtag);
  xtag=(int *)malloc(sizeof(int)*nsearch);
  //
  // create a unique hash
  //
#ifdef TIOGA_HAS_NODEGID
  uniquenode_map(gid_search.data(), res_search, xtag, nsearch);
#else
  uniquenodes_octree(xsearch,tagsearch,res_search,xtag,&nsearch);
#endif
  //
  donorCount=0;
  ipoint=0; 
  dId=(int *) malloc(sizeof(int) *2);
  for(i=0;i<nsearch;i++)
    {
     if (xtag[i]==i) {
	//adt->searchADT(this,&(donorId[i]),&(xsearch[3*i]));
	adt->searchADT(this,dId,&(xsearch[3*i]));
        donorId[i]=dId[0];
      }
      else {
	donorId[i]=donorId[xtag[i]];
      }
      if (donorId[i] > -1) {
	  donorCount++;
	}
       ipoint+=3;
     }
  TIOGA_FREE(dId);
  TIOGA_FREE(icell);
  TIOGA_FREE(obq);
}

void MeshBlock::search_uniform_hex(void)
{
  if (donorId) free(donorId);
  donorId=(int*)malloc(sizeof(int)*nsearch);
  if (xtag) free(xtag);
  xtag=(int *)malloc(sizeof(int)*nsearch);
  //
#ifdef TIOGA_HAS_NODEGID
  uniquenode_map(gid_search.data(), res_search, xtag, nsearch);
#else
  uniquenodes_octree(xsearch,tagsearch,res_search,xtag,&nsearch);
#endif
  //
  int donorCount=0;
  int *dId=(int *) malloc(sizeof(int) *2);
  double xvec[8][3];
  //
  // corners of a cube with of side 4*TOL
  // with origin as the center
  // 
  for(int jj=0;jj<8;jj++)
    for(int k=0;k<3;k++) xvec[jj][k]=(2*((jj & (1 << k)) >> k)-1)*2*TOL;
  //
  double xd[3];
  int dID[2];
  for(int i=0;i<nsearch;i++)
    {
      int idx[3];
      if (xtag[i]==i) {
	for(int j=0;j<3;j++)
	  {
	    xd[j]=0;
	    for(int k=0;k<3;k++)
	      xd[j]+=(xsearch[3*i+k]-xlow[k])*obh->vec[j][k];
            idx[j]=xd[j]/dx[j];
	  }
        if (xd[0] > -TOL && xd[0] < idims[0]*dx[0]+TOL &&
            xd[1] > -TOL && xd[1] < idims[1]*dx[1]+TOL &&
            xd[2] > -TOL && xd[2] < idims[2]*dx[2]+TOL) 
	   {
            for(int k=0;k<3;k++) if (idx[k]==idims[k]) idx[k]--;
            dID[0]=uindx[idx[2]*idims[1]*idims[0]+idx[1]*idims[0]+idx[0]];
            dID[1]=(dID[0] > -1) ? (cellRes[dID[0]]==BIGVALUE) : 1; 
            for(int jj=0;jj<8 && (dId[0]==-1 || dID[1]) ;jj++)
             {
              for(int k=0;k<3;k++)
	       {
                idx[k]=(xd[k]+xvec[jj][k])/dx[k];
                if (idx[k]==idims[k]) idx[k]--;
               }
	        int dtest=uindx[idx[2]*idims[1]*idims[0]+idx[1]*idims[0]+idx[0]];
                dID[1]=(dtest > -1) ? (cellRes[dtest]==BIGVALUE) : 1; 
                dID[0]=(dID[0] == -1) ? dtest : (!dID[1] ? dtest : dID[0]);
              }
             donorId[i]=dID[0]; 
            }
       else
          {
           donorId[i]=-1;
          }
      }
     else {
       donorId[i]=donorId[xtag[i]];
     }
      if (donorId[i] > -1) {
	donorCount++;
      }
      ipoint+=3;
    }
  free(dId);
}
