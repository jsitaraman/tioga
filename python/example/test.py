#!/bin/env python
import sys
import numpy as np
sys.path.append("../../build/python")
import pytioga
from helper import *
from mpi4py import MPI
print "Everything imported correctly"

#                                                                                          
# MPI World Setup                                                                          
#                                                                                          
world  = MPI.COMM_WORLD
rank   = world.Get_rank()
procs  = world.Get_size()

jtot, ktot, ltot = 20, 20, 20

width = 1.0*procs-0.9*rank

if(rank == 0):
    obc    = np.array([], dtype="int32")
else:
    obc    = obcnodes(jtot, ktot, ltot)

coord  = cube(jtot, ktot, ltot, width)
ndc8   = connect(jtot, ktot, ltot)

iblank = np.ones((jtot, ktot, ltot), dtype="int32")

q = np.zeros(5*jtot*ktot*ltot, dtype="double");

# filename = "cube%d.dat"%(rank)
# to_tecplot(filename, coord, iblank)

gridData={'gridtype'        :'unstructured',
          'tetConn'         :[np.array([[],[],[],[]],'i')],
          'pyraConn'        :[np.array([[],[],[],[],[]],'i')],
          'prismConn'       :[np.array([[],[],[],[],[],[]],'i')],
          'hexaConn'        :[ndc8],
          'bodyTag'         :[np.array([int(rank+1)],'i')],
          'wallnode'        :[np.array([], dtype='int')],
          'obcnode'         :[obc],
          'grid-coordinates':[coord],
          'iblanking'       :[iblank],
          'istor'           :'row',
          'q-variables'     :[q]}

t = pytioga.PyTioga(world)

# print "initialized"

t.register_data(gridData)
# print "done registering data!"

t.preprocess()
# print "done preprocessing grids!"

t.connect()
# print "done connecting!"

# try:
#     t.set_data(gridData)
# except:
#     print "error!"

filename = "cube%d.dat"%(rank)
to_tecplot(filename, coord, iblank)
print "wrote data file from python"

