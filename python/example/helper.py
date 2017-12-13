import numpy as np

def cube(jtot, ktot, ltot, width=1.0):
    xyz = np.zeros((ltot,ktot,jtot,3), dtype="double")
    for l in range(ltot):
        for k in range(ktot):
            for j in range(jtot):
                xyz[l,k,j,0] = 1.0*width*j/jtot-0.5*width
                xyz[l,k,j,1] = 1.0*width*k/ktot-0.5*width
                xyz[l,k,j,2] = 1.0*width*l/ltot-0.5*width
    return xyz

def connect(jtot, ktot, ltot):
    ndc8 = np.zeros(((ltot-1)*(ktot-1)*(jtot-1), 8), dtype="int32")
    jstride = 1
    kstride = jtot
    lstride = jtot*ktot
    tidx = 0
    for l in range(ltot-1):
        for k in range(ktot-1):
            for j in range(jtot-1):
                # note +1 since tioga starts at 1
                idx = j*jstride + k*kstride + l*lstride + 1;
                ndc8[tidx,0] = idx
                ndc8[tidx,1] = idx + jstride
                ndc8[tidx,2] = idx + jstride + kstride
                ndc8[tidx,3] = idx + kstride
                ndc8[tidx,4] = idx + lstride
                ndc8[tidx,5] = idx + jstride + lstride
                ndc8[tidx,6] = idx + kstride + jstride + lstride
                ndc8[tidx,7] = idx + kstride + lstride
                tidx += 1
    return ndc8

def obcnodes(jtot, ktot, ltot):
    # num = 2*((jtot-1)*(ktot-1)) + 2*((jtot-1)*(ltot-1)) + 2*((ktot-1)*(ltot-1))
    # num = 2*(jtot*ktot) + 2*(jtot*ltot) + 2*(ktot*ltot)
    # print "number of obcnodes is: ", num
    jstride = 1
    kstride = jtot
    lstride = jtot*ktot
    tidx = 0
    num  = 0
    for l in range(ltot):
        for k in range(ktot):
            for j in range(jtot):
                if(j==0 or j==jtot-1 or
                   k==0 or k==ktot-1 or
                   l==0 or l==ltot-1):
                    num+=1
    obc = np.zeros(num, dtype="int32")
    for l in range(ltot):
        for k in range(ktot):
            for j in range(jtot):
                if(j==0 or j==jtot-1 or
                   k==0 or k==ktot-1 or
                   l==0 or l==ltot-1):
                    obc[tidx] = j*jstride + k*kstride + l*lstride + 1
                    tidx+=1
    return obc

def to_tecplot(filename, xyz, iblank):
    ltot, ktot, jtot = iblank.shape
    with open(filename, "w") as f:
        f.write("VARIABLES = X,Y,Z,IBLANK\n")
        f.write("Zone i=%d,j=%d,k=%d,F=POINT\n"%(jtot, ktot, ltot))
        for l in range(ltot):
            for k in range(ktot):
                for j in range(jtot):
                    f.write("%25.16e %25.16e %25.16e %d\n"%(xyz[l,k,j,0],
                                                            xyz[l,k,j,1],
                                                            xyz[l,k,j,2],iblank[l,k,j]))
