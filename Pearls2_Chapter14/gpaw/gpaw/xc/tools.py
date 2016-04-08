import numpy as np
from ase.units import Hartree

from gpaw.xc import XC
from gpaw.utilities import unpack


def vxc(paw, xc=None, coredensity=True):
    """Calculate XC-contribution to eigenvalues."""
    
    ham = paw.hamiltonian
    dens = paw.density
    wfs = paw.wfs

    if xc is None:
        xc = ham.xc
    elif isinstance(xc, str):
        xc = XC(xc)

    if dens.nt_sg is None:
        dens.interpolate_pseudo_density()

    thisisatest = not True
    
    if xc.orbital_dependent:
        paw.get_xc_difference(xc)

    # Calculate XC-potential:
    vxct_sg = ham.finegd.zeros(wfs.nspins)
    xc.calculate(dens.finegd, dens.nt_sg, vxct_sg)
    vxct_sG = ham.gd.empty(wfs.nspins)
    ham.restrict(vxct_sg, vxct_sG)
    if thisisatest:
        vxct_sG[:] = 1
        
    # ... and PAW corrections:
    dvxc_asii = {}
    for a, D_sp in dens.D_asp.items():
        dvxc_sp = np.zeros_like(D_sp)
        xc.calculate_paw_correction(wfs.setups[a], D_sp, dvxc_sp, a=a, 
                                    addcoredensity=coredensity)
        dvxc_asii[a] = [unpack(dvxc_p) for dvxc_p in dvxc_sp]
        if thisisatest:
            dvxc_asii[a] = [wfs.setups[a].dO_ii]

    vxc_un = np.empty((wfs.kd.mynks, wfs.bd.mynbands))
    for u, vxc_n in enumerate(vxc_un):
        kpt = wfs.kpt_u[u]
        vxct_G = vxct_sG[kpt.s]
        for n in range(wfs.bd.mynbands):
            psit_G = wfs._get_wave_function_array(u, n, realspace=True)
            vxc_n[n] = wfs.gd.integrate((psit_G * psit_G.conj()).real,
                                        vxct_G, global_integral=False)

        for a, dvxc_sii in dvxc_asii.items():
            P_ni = kpt.P_ani[a]
            vxc_n += (np.dot(P_ni, dvxc_sii[kpt.s]) *
                      P_ni.conj()).sum(1).real

    wfs.gd.comm.sum(vxc_un)
    vxc_skn = wfs.kd.collect(vxc_un)

    if xc.orbital_dependent:
        vxc_skn += xc.exx_skn

    return vxc_skn * Hartree
