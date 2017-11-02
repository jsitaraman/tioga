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
		pass
        def sifInitialize(self,properties,conditions):
		pass

        def initData(self):
		self.gridData=pickle.load(open('pickle_dumps/obeGridData'+str(MPI.COMM_WORLD.Get_rank()),'rb'))

        def finish(self):
               pass
