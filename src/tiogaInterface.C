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
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "tioga.h"
#include "globals.h"
#include <string.h>
//
// All the interfaces that are 
// accessible to third party f90 and C
// flow solvers
//
//
// Jay Sitaraman
// 02/24/2014
//
extern "C" {
  void tioga_init_(int *scomm)
  {
    int id_proc,nprocs;
    MPI_Comm tcomm;
    //tcomm=(MPI_Comm) (*scomm);
    tcomm=MPI_Comm_f2c(*scomm);
    //
    tg=new tioga[1];
    //
    //MPI_Comm_rank(MPI_COMM_WORLD,&id_proc);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(tcomm,&id_proc);
    MPI_Comm_size(tcomm,&nprocs);
    //
    tg->setCommunicator(tcomm,id_proc,nprocs);
    nc=NULL;
    nv=NULL;
    vconn=NULL;
  }
  

  void tioga_registergrid_data_(int *btag,int *nnodes,double *xyz,int *ibl,int *nwbc, int *nobc,int *wbcnode, 
			       int *obcnode,int *ntypes,...)
  {
    va_list arguments;
    int i;

    va_start(arguments, ntypes);

    if(nv) free(nv);
    if(nc) free(nc);
    if(vconn) free(vconn);

    nv=(int *) malloc(sizeof(int)*(*ntypes));    
    nc=(int *) malloc(sizeof(int)*(*ntypes));
    vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0;i<*ntypes;i++)
     {
      nv[i]=*(va_arg(arguments, int *));
      nc[i]=*(va_arg(arguments, int *));
      vconn[i]=va_arg(arguments, int *);
     }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,nv,nc,vconn);
  }

  void tioga_register_amr_global_data_(int *nf, int *qstride, double *qnodein,
				      int *idata,double *rdata,
				      int *ngridsin,int *qnodesize)
  {
    tg->register_amr_global_data(*nf,*qstride,qnodein,idata,rdata,*ngridsin,*qnodesize);
  }


  void tioga_register_amr_patch_count_(int *npatches)
  {
    tg->set_amr_patch_count(*npatches);
  }

  void tioga_register_amr_local_data_(int *ipatch,int *global_id,int *iblank,double *q)
  {
    tg->register_amr_local_data(*ipatch,*global_id,iblank,q);
  }

  void tioga_preprocess_grids_(void)
  {
    tg->profile();
  }

  void tioga_performconnectivity_(void)
  {
    tg->performConnectivity();
  }

  void tioga_performconnectivity_highorder_(void)
  {
    tg->performConnectivityHighOrder();
  }

  void tioga_performconnectivity_amr_(void)
  {
   tg->performConnectivityAMR();
  }

  void tioga_dataupdate_(double *q,int *nvar,char *itype)
  {
    int interptype;
    if (strstr(itype,"row")) 
      {
	interptype=0;
      }
    else if (strstr(itype,"column")) 
      {
	interptype=1;
      }
    else
      {
	printf("#tiogaInterface.C:dataupdate_:unknown data orientation\n");
	return;
      }
    if (tg->ihighGlobal==0) 
      {
	tg->dataUpdate(*nvar,q,interptype);
      }
    else
      {
	tg->dataUpdate_highorder(*nvar,q,interptype);
      }
  }
  void tioga_writeoutputfiles_(double *q,int *nvar,char *itype)
  {
    int interptype;
    if (strstr(itype,"row")) 
      {
	interptype=0;
      }
    else if (strstr(itype,"column")) 
      {
	interptype=1;
      }
    else
      {
	printf("#tiogaInterface.C:dataupdate_:unknown data orientation\n");
	return;
      }
    tg->writeData(*nvar,q,interptype);
  }    
  void tioga_getdonorcount_(int *dcount,int *fcount)
  {
    tg->getDonorCount(dcount,fcount);
  }
  void tioga_getdonorinfo_(int *receptors,int *indices,double *frac,int *dcount)
  {
    tg->getDonorInfo(receptors,indices,frac,dcount);
  }

  void tioga_setsymmetry_(int *isym)
  {
    tg->setSymmetry(*isym);
  }

  void tioga_setresolutions_(double *nres,double *cres)
  {
    tg->setResolutions(nres,cres);
  }
  
  void tioga_setcelliblank_(int *iblank_cell)
  {
    tg->set_cell_iblank(iblank_cell);
  }

  void tioga_set_highorder_callback_(void (*f1)(int*, int*),
				    void (*f2)(int *,int *,double *),
				    void (*f3)(int *,double *,int *,double *),
				    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
				     void (*f5)(int *,int *,double *,int *,int *,double *))
  {
    tg->setcallback(f1,f2,f3,f4,f5);
    //    get_nodes_per_cell=f1;
    //get_receptor_nodes=f2;
    //donor_inclusion_test=f3;
    //donor_frac=f4;
    //convert_to_modal=f5;
  }
  void tioga_set_amr_callback_(void (*f1)(int *,double *,int *,double *))
  {
    tg->set_amr_callback(f1);
  }
  
  void tioga_delete_(void)
   {
    delete [] tg;
    if (nc) free(nc);
    if (nv) free(nv);
    if (vconn) free(vconn);
   }
}
