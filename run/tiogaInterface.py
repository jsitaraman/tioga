#! /usr/bin/env python
#
# File: tiogaInterface.py
# Authors: Jacob Crabill
# Last Modified: 9/03/2017

__version__ = "1.00"


# ==================================================================
# Standard Python modules

import sys
import os
import copy
import string
import types
import pickle
#Extension modules
#sys.path.append(os.path.abspath(os.getenv('PYTHONPATH')))
#sys.path.append(os.path.abspath(os.getenv('PYTHONPATH')+'/numpy'))

import numpy as np

CONVERT_DIR = '/home/jcrabill/zefr/external/'
sys.path.append(CONVERT_DIR)

from convert import *

# Try to import the MPI.COMM_WORLD.module
try:
    from mpi4py import MPI
    _parallel = True

except ImportError:
    _parallel = False
    raise


try:
    import tioga as tg
except ImportError:
    print("Import Error: TIOGA connectivity module")
    raise

class Tioga:

    def __init__(self,gridID,nGrids):

        self.nGrids = nGrids
        self.ID = gridID

        Comm = MPI.COMM_WORLD
        rank = Comm.Get_rank()
        nproc = Comm.Get_size()

        self.gridComm = Comm.Split(gridID,rank)
        self.gridRank = self.gridComm.Get_rank()
        self.gridSize = self.gridComm.Get_size()

        tg.tioga_init_(Comm)

        self.name = 'tioga'

    # Initial grid preprocessing
    def preprocess(self):
        tg.tioga_preprocess_grids_()

    # Perform the full domain connectivity
    def performConnectivity(self):
        tg.tioga_performconnectivity_()

    # Perform just the donor/fringe point connecitivity (using current blanking)
    def performPointConnectivity(self):
        tg.tioga_do_point_connectivity()

    # For high-order codes: First part of unblank procedure (t^{n+1} blanking)
    def unblankPart1(self):
        tg.tioga_unblank_part_1()

    # For high-order codes: Second part of unblank procedure (t^n blanking + union)
    def unblankPart2(self):
        tg.tioga_unblank_part_2(self.nfields)

    # Interpolate solution and send/receive all data
    def exchangeSolution(self):
        tg.tioga_dataupdate_ab(self.nfields, 0)

    # Interpolate solution gradient and send/receive all data
    def exchangeGradient(self):
        tg.tioga_dataupdate_ab(self.nfields, 1)

    def sifInitialize(self,properties,conditions):
        self.ndims = int(properties['ndims'])
        self.nfields = int(properties['nfields'])
        self.motion = int(properties['moving-grid'])
        self.gridScale = conditions['meshRefLength']

        try:
            self.useGpu = int(properties['use-gpu'])
        except:
            self.useGpu = False

        # Have TIOGA do any additional setup on input parameters
        #tg.sifinit(dt,Re,mach_fs)

    # sifInit called first to setup simulation input
    def initData(self,gridData,callbacks):
        self.gridData = gridData
        self.callbacks = callbacks

        # Get pointers to grid data
        btag = gridData['bodyTag'][0]
        xyz = arrayToDblPtr(gridData['grid-coordinates'][0])
        c2v = arrayToIntPtr(gridData['hexaConn'][0])
        iblank = arrayToIntPtr(gridData['iblanking'][0])
	#pickle.dump(gridData['grid-coordinates'],open('xyz'+str(MPI.COMM_WORLD.Get_rank()),'wb'))
        overNodes = arrayToIntPtr(gridData['obcnode'][0])
        wallNodes = arrayToIntPtr(gridData['wallnode'][0])

        c2f = arrayToIntPtr(gridData['cell2face'][0])
        f2c = arrayToIntPtr(gridData['face2cell'][0])

        iblank = arrayToIntPtr(gridData['iblanking'][0])
        iblank_face = arrayToIntPtr(gridData['iblank-face'][0])
        iblank_cell = arrayToIntPtr(gridData['iblank-cell'][0])

        overFaces = arrayToIntPtr(gridData['overset-faces'][0])
        wallFaces = arrayToIntPtr(gridData['wall-faces'][0])

        mpiFaces = arrayToIntPtr(gridData['mpi-faces'][0])
        mpiProcR = arrayToIntPtr(gridData['mpi-right-proc'][0])
        mpiFidR = arrayToIntPtr(gridData['mpi-right-id'][0])

        f2v = arrayToIntPtr(gridData['faceConn'][0])

        # Extract metadata
        nCellTypes = 1
        nFaceTypes = 1
        nnodes = gridData['grid-coordinates'][0].shape[0]
        ncells = gridData['hexaConn'][0].shape[0]
        nvert = gridData['hexaConn'][0].shape[1]
        nfaces = gridData['faceConn'][0].shape[0]
        nvertf = gridData['faceConn'][0].shape[1]

        nover = gridData['obcnode'][0].shape[0]
        nwall = gridData['wallnode'][0].shape[0]

        nOverFace = gridData['overset-faces'][0].shape[0]
        nWallFace = gridData['wall-faces'][0].shape[0]
        nMpiFace = gridData['mpi-faces'][0].shape[0]

        gridType = gridData['gridCutType']

        tg.tioga_registergrid_data_(btag, nnodes, xyz, iblank,
            nwall, nover, wallNodes, overNodes, nCellTypes, nvert,
            ncells, c2v)

        tg.tioga_setcelliblank_(iblank_cell)

        tg.tioga_register_face_data_(gridType, f2c, c2f, iblank_face,
            nOverFace, nWallFace, nMpiFace, overFaces, wallFaces,
            mpiFaces, mpiProcR, mpiFidR, nFaceTypes, nvertf, nfaces, f2v);

        # Get solver callbacks
        get_nodes_per_cell = callbacks['nodesPerCell']
        get_nodes_per_face = callbacks['nodesPerFace']
        get_receptor_nodes = callbacks['receptorNodes']
        get_face_nodes = callbacks['faceNodes']
        donor_inclusion_test = callbacks['donorInclusionTest']
        convert_to_modal = callbacks['convertToModal']
        donor_frac = callbacks['donorFrac']

        get_q_spt = callbacks['get_q_spt']
        get_q_fpt = callbacks['get_q_fpt']
        get_dq_spt = callbacks['get_dq_spt']
        get_dq_fpt = callbacks['get_dq_fpt']
        get_q_spts = callbacks['get_q_spts']
        get_dq_spts = callbacks['get_dq_spts']

        tg.tioga_set_highorder_callback_(get_nodes_per_cell,
            get_receptor_nodes, donor_inclusion_test, donor_frac,
            convert_to_modal)

        tg.tioga_set_ab_callback_(get_nodes_per_face, get_face_nodes,
            get_q_spt, get_q_fpt, get_dq_spt, get_dq_fpt, get_q_spts,
            get_dq_spts)

        if self.motion:
            gridV = arrayToDblPtr(gridData['gridVel'][0])
            offset = arrayToDblPtr(gridData['rigidOffset'][0])
            Rmat = arrayToDblPtr(gridData['rigidRotMat'][0])
            tg.tioga_register_moving_grid_data(gridV,offset,Rmat)

        if self.useGpu:
            donorFromDevice = callbacks['donorDataDevice']
            fringeToDevice = callbacks['fringeDataToDevice']
            unblankToDevice = callbacks['unblankToDevice']
            faceNodesGPU = callbacks['faceNodesGPU']
            cellNodesGPU = callbacks['cellNodesGPU']
            qSpts_d = callbacks['q_spts_d']
            dqSpts_d = callbacks['dq_spts_d']
            nWeightsGPU = callbacks['nWeightsGPU']
            weightsGPU = callbacks['weightsGPU']

            tg.tioga_set_ab_callback_gpu_(donorFromDevice, fringeToDevice,
                unblankToDevice, qSpts_d, dqSpts_d, faceNodesGPU, cellNodesGPU,
                nWeightsGPU, weightsGPU)

            coords_d = arrayToDblPtr(gridData['nodesGPU'][0])
            ecoords_d = arrayToDblPtr(gridData['eleCoordsGPU'][0])
            iblankCell_d = arrayToIntPtr(gridData['iblankCellGPU'][0])
            iblankFace_d = arrayToIntPtr(gridData['iblankFaceGPU'][0])

            tg.tioga_set_device_geo_data(coords_d, ecoords_d, iblankCell_d,
                iblankFace_d)

            tb.tioga_set_stream_handle(gridData['cuStream'],gridData['cuEvent'])
            
    def initAMRData(self,gridData):
        pickle.dump(gridData,open('obeGridData'+str(MPI.COMM_WORLD.Get_rank()),'wb'))
        ngrids=len(gridData['gridParam'])
        local2global_tmp=np.zeros((ngrids,),'i')
        local2global=np.zeros((ngrids,),'i')
        for i in range(len(gridData['qParam'])):
            local2global_tmp[gridData['qParam'][i][0]]=i
        MPI.COMM_WORLD.Allreduce(local2global_tmp,local2global)
        #
        # these variables are from high-order FE off-body implementation
        # have to see if p=0 defaults work
        #
        nf=3
        qstride=1
        qnodein=0.0
        qnodeinC=arrayToDblPtr(np.array([qnodein,0.0],'d'))
        qnodesize=1
        #
        idata=np.zeros((ngrids*11),'i')
        rdata=np.zeros((ngrids*6),'d')
        for i in range(ngrids):
            m=i*11
            idata[m]   =gridData['gridParam'][i][0] # global id of patch
            idata[m+1] =gridData['gridParam'][i][1] # level number of patch 
            idata[m+2] =gridData['gridParam'][i][3] # proc id containing patch
            idata[m+3] =0                        # porder (set to zero for FD) 
            idata[m+4] =local2global[idata[m]]   # local number of patch  
            idata[m+5] =gridData['ilo'][i][0]       # lower left-hand front global numbering (x)
            idata[m+6] =gridData['ilo'][i][1]       # lower left-hand front global numbering (y)
            idata[m+7] =gridData['ilo'][i][2]       # lower left-hand front global numbering (z)
            idata[m+8] =gridData['ihi'][i][0]       # upper right-hand back global numbering (x) 
            idata[m+9] =gridData['ihi'][i][1]       # lower right-hand back global numbering (y)
            idata[m+10] =gridData['ihi'][i][2]      # lower left-hand  back global numbering (z)
            m1=i*6
            rdata[m1]   = gridData['xlo'][i][0]
            rdata[m1+1] = gridData['xlo'][i][1]
            rdata[m1+2] = gridData['xlo'][i][2]
            rdata[m1+3] = gridData['dx'][i]
            rdata[m1+4] = gridData['dx'][i]
            rdata[m1+5] = gridData['dx'][i]
        idataC=arrayToIntPtr(idata)
        rdataC=arrayToDblPtr(rdata)
        tg.tioga_register_amr_global_data_(nf,qstride,qnodeinC,idataC,rdataC,ngrids,qnodesize)
        #
        npatches=len(gridData['qParam'])
        tg.tioga_register_amr_patch_count_(npatches)
        for i in range(npatches):
            global_id=gridData['qParam'][i][0]
            qC=arrayToDblPtr(gridData['q-variables'][i])
            iblankC=arrayToIntPtr(gridData['iblanking'][i])
            tg.tioga_register_amr_local_data_(i,global_id,iblankC,qC)

    def performAMRConnectivity(self):
        tg.tioga_performconnectivity_amr_()
        
    def finish(self,step):
        tg.tioga_delete_()
