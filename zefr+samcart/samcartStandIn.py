# stand in for samcart

import sys
import os
import numpy
from mpi4py import MPI
import pickle
class samcartSolver:

    def  __init__(self,sif=None):
        self.worldComm = MPI.COMM_WORLD
        self.name='samcart'

    def initialize(self,startfile,conditions):
        try:
            self.porder = int(conditions['porder'])
        except:
            self.porder = 1
        pass

    def sifInitialize(self,properties,conditions):
        pass

    def initData(self):
        self.gridData=pickle.load(open('pickle_dumps/obeGridData'+str(MPI.COMM_WORLD.Get_rank()),'rb'))
        self.gridData['porder'] = self.porder;

    def finish(self):
        pass
