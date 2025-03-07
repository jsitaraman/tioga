// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "codetypes.h"
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

  void tioga_init_f90_(int *scomm)
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
    for(int i=0;i<MAXBLOCKS;i++)
     {
      idata[i].nc=NULL;
      idata[i].nv=NULL;
      idata[i].vconn=NULL;
     }
  }

  void tioga_init_(MPI_Comm tcomm)
  {
    int id_proc,nprocs;
    //MPI_Comm tcomm;
    //tcomm=(MPI_Comm) (*scomm);
    //tcomm=MPI_Comm_f2c(*scomm);
    //
    tg=new tioga[1];
    //
    //MPI_Comm_rank(MPI_COMM_WORLD,&id_proc);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(tcomm,&id_proc);
    MPI_Comm_size(tcomm,&nprocs);
    //
    tg->setCommunicator(tcomm,id_proc,nprocs);
    for (int i=0;i<MAXBLOCKS;i++)
     {
      idata[i].nc=NULL;
      idata[i].nv=NULL;
      idata[i].vconn=NULL;
     }
  }
  
  void tioga_registergrid_data_(int *btag,int *nnodes,double *xyz,int *ibl,int *nwbc, int *nobc,int *wbcnode,
                               int *obcnode,int *ntypes,...)
   {
    va_list arguments;
    int i;
    int iblk=0;
    va_start(arguments, ntypes);

    if(idata[iblk].nv) TIOGA_FREE(idata[iblk].nv);
    if(idata[iblk].nc) TIOGA_FREE(idata[iblk].nc);
    if(idata[iblk].vconn) TIOGA_FREE(idata[iblk].vconn);
    idata[iblk].nv=(int *) malloc(sizeof(int)*(*ntypes));    
    idata[iblk].nc=(int *) malloc(sizeof(int)*(*ntypes));
    idata[iblk].vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0;i<*ntypes;i++)
     {
      idata[iblk].nv[i]=*(va_arg(arguments, int *));
      idata[iblk].nc[i]=*(va_arg(arguments, int *));
      idata[iblk].vconn[i]=va_arg(arguments, int *);
     }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,idata[iblk].nv,idata[iblk].nc,idata[iblk].vconn);
  }


  void tioga_registergrid_data_mb_(int *bid, int *btag,int *nnodes,double *xyz,int *ibl,int *nwbc, int *nobc,int *wbcnode, 
			       int *obcnode,int *ntypes,...)
  {
    va_list arguments;
    int i;
    int iblk=*bid-BASE;

    va_start(arguments, ntypes);

    if(idata[iblk].nv) TIOGA_FREE(idata[iblk].nv);
    if(idata[iblk].nc) TIOGA_FREE(idata[iblk].nc);
    if(idata[iblk].vconn) TIOGA_FREE(idata[iblk].vconn);
    idata[iblk].nv=(int *) malloc(sizeof(int)*(*ntypes));    
    idata[iblk].nc=(int *) malloc(sizeof(int)*(*ntypes));
    idata[iblk].vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0;i<*ntypes;i++)
     {
      idata[iblk].nv[i]=*(va_arg(arguments, int *));
      idata[iblk].nc[i]=*(va_arg(arguments, int *));
      idata[iblk].vconn[i]=va_arg(arguments, int *);
     }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,idata[iblk].nv,idata[iblk].nc,idata[iblk].vconn);
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

  void tioga_registersolution_(int *bid,double *q)
  {
    tg->registerSolution(*bid,q);
  }

  void tioga_dataupdate_mb_(int *nvar,char *itype)
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

    //tg->dataUpdate(*nvar,interptype);

    if (tg->ihighGlobal==0) 
    {
      if (tg->iamrGlobal==0) 
      	{
     	    tg->dataUpdate(*nvar,interptype);
        }
     	else
     	  {
     	    tg->dataUpdate_AMR(*nvar,interptype);
     	  }
    }
    else
    {
     if (tg->iamrGlobal==0) 
       {
     	    tg->dataUpdate(*nvar,interptype,1);
       }
      else
      {
        printf("Data udpate between high-order near-body and AMR cartesian Not implemented yet\n");
      }
    }
  }

  void tioga_dataupdate_(double *q,int *nvar,char *itype)
  {
    int interptype;
    int bid=0;
    tg->registerSolution(bid,q);
    tioga_dataupdate_mb_(nvar,itype);
  }

  void tioga_writeoutputfiles_(int *nvar,char *itype)
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
    tg->writeData(*nvar,interptype);
  }    
  void tioga_getdonorcount_(int *btag,int *dcount,int *fcount)
  {
    tg->getDonorCount(*btag,dcount,fcount);
  }
  void tioga_getdonorinfo_(int *btag,int *receptors,int *indices,double *frac,int *dcount)
  {
    tg->getDonorInfo(*btag,receptors,indices,frac,dcount);
  }

  void tioga_setsymmetry_(int *isym)
  {
    tg->setSymmetry(*isym);
  }

  void tioga_setresolutions_(double *nres,double *cres)
  {
    tg->setResolutions(nres,cres);
  }

  void tioga_setresolutions_multi_(int *btag,double *nres,double *cres)
  {
    tg->setResolutions(*btag,nres,cres);
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
  
  void tioga_set_p4est_(void)
  {
    tg->set_p4est();
  }
  void tioga_set_amr_callback_(void (*f1)(int *,double *,int *,double *))
  {
    tg->set_amr_callback(f1);
  }
  void tioga_set_p4est_search_callback_(void (*f1)(double *xsearch,int *process_id,int *cell_id,int *npts),
					void (*f2)(int *pid,int *iflag))
  {
    tg->setp4estcallback(f1,f2);
  //jayfixme  tg->set_p4est_search_callback(f1);
  }  

  void tioga_reduce_fringes_(void)
  {
    tg->reduce_fringes();
  }

  void tioga_setnfringe_(int *nfringe)
  {
    tg->setNfringe(nfringe);
  }

  void tioga_setmexclude_(int *mexclude)
  {
   tg->setMexclude(mexclude);
  }

  void tioga_delete_(void)
   {
    delete [] tg;
    for(int i=0;i<MAXBLOCKS;i++)
     {
       if (idata[i].nc) TIOGA_FREE(idata[i].nc);
       if (idata[i].nv) TIOGA_FREE(idata[i].nv);
       if (idata[i].vconn) TIOGA_FREE(idata[i].vconn);
     }
   }

}
