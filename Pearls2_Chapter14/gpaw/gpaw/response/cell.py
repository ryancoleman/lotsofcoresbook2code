import numpy as np
from math import sqrt, pi

def get_primitive_cell(a, rpad=None):
    """From the unit cell, calculate primitive cell and volume. """

    if rpad is None:
        rpad = np.ones(3, int)

    a = (a.T * rpad).T # ???
    vol = np.abs(np.dot(a[0], np.cross(a[1], a[2])))
    BZvol = (2. * pi)**3 / vol
    b = np.linalg.inv(a.T)
    
    b *= 2 * pi

    assert np.abs((np.dot(b.T, a) - 2.*pi*np.eye(3)).sum()) < 1e-10

    return a, b, vol, BZvol


def set_Gvectors(acell, bcell, nG, Ecut, q=[0., 0., 0.]):
    """Calculate the number of planewaves with a certain cutoff, their reduced coordinates and index."""

    # Refer to R.Martin P85
    Gmax = np.zeros(3, dtype=int)
    Gcut = np.zeros(3, dtype=float)
    for i in range(3):
        a = acell[i]
        Gcut[i] = sqrt(2*Ecut[i])
        if Gcut[i] > 0:
            Gmax[i] = sqrt(a[0]**2 + a[1]**2 + a[2]**2) * Gcut[i] / (2*pi) + 1

    Nmax = np.ones([3], dtype = int)   
    for dim in range(3):
        if Gcut[dim] > 0:
            Nmax[dim] = 2 * Gmax[dim] + 2
    #assert (nG - Nmax >=0).all() # to prevent too many planewaves     
    m = {}
    for dim in range(3):
        m[dim] = np.zeros(Nmax[dim],dtype=int)
        for i in range(Nmax[dim]):
            m[dim][i] = i
            if m[dim][i] > np.int(Gmax[dim]):
                m[dim][i] = i - Nmax[dim]
                
    G = np.zeros((Nmax[0]*Nmax[1]*Nmax[2],3),dtype=int)
    n = 0
    for i in range(Nmax[0]):
        for j in range(Nmax[1]):
            for k in range(Nmax[2]):
                tmp = np.array([m[0][i], m[1][j], m[2][k]])
                tmpG = np.dot(tmp+np.array(q), bcell)
                Gmod = 0
                for dim in range(3):
                    if Gcut[dim] > 0:
                        Gmod += tmpG[dim]**2/Gcut[dim]**2
                if Gmod < 1:
                    G[n] = tmp
                    n += 1
    npw = n
    Gvec = G[:n]

    Gindex = np.zeros(npw, dtype=int) 
    id = np.zeros(3, dtype=int)

    for iG in range(npw):
        G = Gvec[iG]
        for dim in range(3):
            if G[dim] >= 0:
                id[dim] = G[dim]
            else:
                id[dim] = nG[dim] - np.abs(G[dim])
        Gindex[iG] = id[0]*nG[1]*nG[2] + id[1]*nG[2] + id[2] 

    return npw, Gvec, Gindex
