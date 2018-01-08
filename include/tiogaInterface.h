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
#include <mpi.h>
#ifdef _GPU
#include <cuda_runtime.h>
#endif
//
// All the interfaces that are 
// accessible to third party f90 and C
// flow solvers
//
// This file is intended to be process by SWIG to create a Python interface
//
// Jay Sitaraman
// 02/24/2014
//
extern "C" {

void tioga_init_f90_(int *scomm);
void tioga_init_(MPI_Comm tcomm);

void tioga_registergrid_data_(int btag, int nnodes, double *xyz, int *ibl,
                              int nwbc, int nobc, int *wbcnode, int *obcnode,
                              int ntypes, int *_nv, int *_ncf, int* _nc, int **_vconn);

void tioga_register_face_data_(int gtype, int *f2c, int **c2f, int *fibl, int nOverFaces,
    int nWallFaces, int nMpiFaces, int *overFaces, int *wallFaces, int *mpiFaces,
    int* mpiProcR, int* mpiFidR, int nftype, int *_nfv, int *_nf, int **_fconn);

void tioga_register_amr_global_data_(int nf, int qstride, double *qnodein,
				     int *idata,double *rdata,
				     int ngridsin,int qnodesize);

void tioga_register_amr_patch_count_(int npatches);

void tioga_register_amr_local_data_(int ipatch,int global_id,int *iblank,double *q);

void tioga_preprocess_grids_(void);

void tioga_performconnectivity_(void);

void tioga_performconnectivity_highorder_(void);

void tioga_performconnectivity_amr_(void);

void tioga_dataupdate_(double *q,int *nvar,char *itype);

void tioga_dataupdate_amr(double *q, int nvar, int interptype);

void tioga_dataupdate_ab(int nvar, int gradFlag = 0);

void tioga_dataupdate_ab_send(int nvar, int gradFlag = 0);
void tioga_dataupdate_ab_recv(int nvar, int gradFlag = 0);

void tioga_writeoutputfiles_(double *q,int *nvar,char *itype);

void tioga_getdonorcount_(int *dcount,int *fcount);

void tioga_getdonorinfo_(int *receptors,int *indices,double *frac,int *dcount);

void tioga_setsymmetry_(int *isym);

void tioga_setresolutions_(double *nres,double *cres);

void tioga_setcelliblank_(int *iblank_cell);

void tioga_set_highorder_callback_(void (*f1)(int*, int*),
                                   void (*f2)(int *,int *,double *),
                                   void (*f3)(int *,double *,int *,double *),
                                   void (*f4)(int *,double *,int *,int *,double *,double *,int *),
                                   void (*f5)(int *,int *,double *,int *,int *,double *));

void tioga_set_ab_callback_(void (*gnf)(int* id, int* npf),
                            void (*gfn)(int* id, int* npf, double* xyz),
                            double (*gqs)(int ic, int spt, int var),
                            double& (*gqf)(int ff, int fpt, int var),
                            double (*ggs)(int ic, int spt, int dim, int var),
                            double& (*ggf)(int ff, int fpt, int dim, int var),
                            double* (*gqss)(int& es, int& ss, int& vs, int etype),
                            double* (*gdqs)(int& es, int& ss, int& vs, int& ds, int etype));

void tioga_set_ab_callback_gpu_(void (*d2h)(int* ids, int nd, int grad),
                                void (*h2df)(int* ids, int nf, int grad, double *data),
                                void (*h2dc)(int* ids, int nc, int grad, double *data),
                                double* (*gqd)(int& es, int& ss, int& vs, int etype),
                                double* (*gdqd)(int& es, int& ss, int& vs, int& ds, int etype),
                                void (*gfng)(int* ids, int nf, int* nptf, double* xyz),
                                void (*gcng)(int* ids, int nf, int* nptf, double* xyz),
                                int (*gnw)(int),
                                void (*dfg)(int*, int, double*, double*));


void tioga_register_moving_grid_data(double* grid_vel, double* offset, double* Rmat);

void tioga_set_transform(double *mat, double *off, int ndim);

void tioga_do_point_connectivity(void);

void tioga_set_iter_iblanks(double dt, int nvar);

void tioga_unblank_part_1(void);
void tioga_unblank_part_2(int nvar);

void tioga_set_amr_callback_(void (*f1)(int *,double *,int *,double *));

void tioga_delete_(void);

//! For GPU-based interpolation, pointers to the CUDA stream & event to use
void tioga_set_stream_handle(void* stream, void* event);

void tioga_set_device_geo_data(double* xyz, double* coord, int* ibc, int* ibf);
} /* extern "C" */
