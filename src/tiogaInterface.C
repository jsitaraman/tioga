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

  void tioga_init_f90_(int *scomm)
  {
    int id_proc,nprocs;
    MPI_Comm tcomm = MPI_Comm_f2c(*scomm);
    //
    tg=new tioga();
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

  void tioga_init_(MPI_Comm tcomm)
  {
    int id_proc,nprocs;

    tg = new tioga();

    MPI_Comm_rank(tcomm,&id_proc);
    MPI_Comm_size(tcomm,&nprocs);

    tg->setCommunicator(tcomm,id_proc,nprocs);
    nc=NULL;
    nv=NULL;
    vconn=NULL;
  }
  
  void tioga_registergrid_data_(int btag, int nnodes, double *xyz, int *ibl,
                                int nwbc, int nobc, int *wbcnode, int *obcnode,
                                int ntypes, int _nv, int _nc, int *_vconn)
  {
    free(nv);
    free(nc);
    free(vconn);

    // NOTE: due to SWIG/Python wrapping, removing va_args stuff for now
    nv    = (int *) malloc(sizeof(int)*ntypes);
    nc    = (int *) malloc(sizeof(int)*ntypes);
    vconn = (int **)malloc(sizeof(int *)*ntypes);
    nv[0] = _nv;
    nc[0] = _nc;
    vconn[0] = _vconn;

    tg->registerGridData(btag,nnodes,xyz,ibl,nwbc,nobc,wbcnode,obcnode,ntypes,nv,nc,vconn);
  }

  void tioga_register_face_data_(int *f2c, int *c2f, int *fibl, int nOverFaces,
                                 int nMpiFaces, int *overFaces, int *mpiFaces,
                                 int* mpiProcR, int* mpiFidR, int nftype,
                                 int _nfv, int _nf, int *_fconn)
  {
    free(nfv);
    free(nf);
    free(fconn);

    nfv   = (int *) malloc(sizeof(int)*nftype);
    nf    = (int *) malloc(sizeof(int)*nftype);
    fconn = (int **)malloc(sizeof(int *)*nftype);
    nfv[0] = _nfv;
    nf[0] = _nf;
    fconn[0] = _fconn;

    tg->registerFaceConnectivity(nftype, nf, nfv, fconn, f2c, c2f, fibl,
                                 nOverFaces, nMpiFaces, overFaces, mpiFaces,
                                 mpiProcR, mpiFidR);
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
	if (tg->iamrGlobal==0) 
	  {
	    tg->dataUpdate(*nvar,q,interptype);
	  }
	else
	  {
	    tg->dataUpdate_AMR(*nvar,q,interptype);
	  }
      }
    else
      {
	if (tg->iamrGlobal==0) 
	  {
	    tg->dataUpdate_highorder(*nvar,q,interptype);
	  }
	else
	  {
	    printf("Data udpate between high-order near-body and AMR cartesian Not implemented yet\n");
	  }
      }
  }

  void tioga_dataupdate_ab(int nvar, double *q_spts, int gradFlag)
  {
    tg->dataUpdate_artBnd(nvar, q_spts, gradFlag);
  }

  void tioga_dataupdate_ab_send(int nvar, int gradFlag)
  {
    tg->dataUpdate_artBnd_send(nvar, gradFlag);
  }

  void tioga_dataupdate_ab_recv(int nvar, int gradFlag)
  {
    tg->dataUpdate_artBnd_recv(nvar, gradFlag);
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
    //get_nodes_per_cell=f1;
    //get_receptor_nodes=f2;
    //donor_inclusion_test=f3;
    //donor_frac=f4;
    //convert_to_modal=f5;
  }

  void tioga_set_ab_callback_(void (*gnf)(int* id, int* npf),
                              void (*gfn)(int* id, int* npf, double* xyz),
                              double (*gqs)(int ic, int spt, int var),
                              double& (*gqf)(int ff, int fpt, int var),
                              double (*ggs)(int ic, int spt, int dim, int var),
                              double& (*ggf)(int, int, int, int),
                              double* (*gqss)(int& es, int& ss, int& vs),
                              double* (*gdqs)(int& es, int& ss, int& vs, int& ds))
  {
    tg->set_ab_callback(gnf, gfn, gqs, gqf, ggs, ggf, gqss, gdqs);
  }

  void tioga_set_ab_callback_gpu_(void (*d2h)(int* ids, int nd, int grad),
                                  void (*h2df)(int* ids, int nf, int grad, double *data),
                                  void (*h2dc)(int* ids, int nc, int grad, double *data),
                                  double* (*gqd)(int&, int&, int&),
                                  double* (*gdqd)(int&, int&, int&, int&))
  {
    tg->set_ab_callback_gpu(d2h,h2df,h2dc,gqd,gdqd);
  }

  void tioga_register_moving_grid_data(double* grid_vel)
  {
    tg->registerMovingGridData(grid_vel);
  }

  void tioga_set_amr_callback_(void (*f1)(int *,double *,int *,double *))
  {
    tg->set_amr_callback(f1);
  }

  void tioga_set_transform(double *rmat, double *offset, int ndim)
  {
    tg->setTransform(rmat, offset, ndim);
  }
  
  void tioga_do_point_connectivity(void)
  {
    tg->doPointConnectivity();
  }

  void tioga_set_iter_iblanks(double dt, int nvar)
  {
    tg->setIterIblanks(dt, nvar);
  }

  void tioga_delete_(void)
   {
    delete tg;
    free(nc);
    free(nv);
    free(vconn);
    free(nf);
    free(nfv);
    free(fconn);
   }

  void tioga_set_stream_handle(cudaStream_t stream, cudaEvent_t event)
  {
    tg->set_stream_handle(stream, event);
  }
}
