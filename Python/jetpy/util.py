import numpy as np
import h5py as h5


def loadCheckpoint(filename):

    f = h5.File(filename, "r")
    Index = f['Grid/Index'][...]
    Nr = f['Grid/Nr'][...]
    T = f['Grid/T'][...][0]
    t_jph = f['Grid/t_jph'][...]

    cells = f['Data/Cells'][...]
    f.close()

    Index = Index[:, 0]
    Nr = Nr[:, 0]

    prim = cells[:, :-1]
    r_iph = cells[:, -1]

    r = np.empty(r_iph.shape)
    th = np.empty(r_iph.shape)
    dV = np.empty(r_iph.shape)

    Nth = Nr.shape[0]

    for i in range(Nth):
        ia = Index[i]
        ib = Index[i] + Nr[i]
        rf = np.empty(Nr[i]+1)
        rf[0] = 0
        rf[1:] = r_iph[ia:ib]
        r[ia:ib] = 0.5*(rf[:-1] + rf[1:])
        th[ia:ib] = 0.5*(t_jph[i]+t_jph[i+1])

        rb = rf[1:]
        ra = rf[:-1]
        dr = rb-ra
        sindth = np.cos(t_jph[i]) - np.cos(t_jph[i+1])

        dV[ia:ib] = (ra*ra+ra*rb+rb*rb)*dr*sindth/3.0

    return T, r, th, prim, dV, (r_iph, t_jph, Index, Nr)


def getTime(filename):

    f = h5.File(filename, "r")
    T = f['Grid/T'][0]
    f.close()
    return T


def getRminmax(filename):

    f = h5.File(filename, "r")
    riph = f['Data/Cells'][:, -1][...]
    f.close()
    return 2*riph.min(), riph.max()
