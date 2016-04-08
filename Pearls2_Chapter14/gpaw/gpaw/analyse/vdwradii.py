import numpy as np
from ase.data import atomic_numbers, chemical_symbols
from ase.units import Bohr

from gpaw.setup import Setups
from gpaw.xc import XC
from gpaw.mpi import world

Bondi64jpc_vdWradii = {  # units Anstrom
    'He': 1.40,
    'Ne': 1.54,
    'Ar': 1.88,
    'Kr': 2.02,
    'Xe': 2.16
}
# Van der Waals Radii after
# Pekka Pyykko, Chem. Rev. 97 (1997) 597-636
# Table 2
Pyykko97cr_vdWradii = {  # units Anstrom
    'Ne': 1.55,
    'Ar': 1.88,
    'Kr': 2.00,
    'Xe': 2.18,
    'Rn': 2.24
}
collected_vdWradii = Bondi64jpc_vdWradii
collected_vdWradii['Rn'] = Pyykko97cr_vdWradii['Rn']


def vdWradii(symbols, xc):
    """Find the elements van der Waals radius.

    Method proposed in:
    Tkatchenko and Scheffler PRL 102 (2009) 073005

    The returned radii are given in Angstroms.
    """
    Z_rare_gas = [atomic_numbers[symbol] for symbol in Bondi64jpc_vdWradii]
    Z_rare_gas.append(atomic_numbers['Rn'])
    Z_rare_gas.sort()

    if isinstance(xc, str):
        xc = XC(xc)

    def get_density(Z):
        """Return density and radial grid from setup."""
        # load setup
        setups = Setups([Z], 'paw', {}, 2,
                        xc, world)
        setup = setups[0].data
        #  create density
        n_g = setup.nc_g.copy()
        for f, phi_g in zip(setup.f_j, setup.phi_jg):
            n_g += f * phi_g ** 2
        return n_g, setup.rgd.r_g

    radii = []
    radius = {}
    for symbol in symbols:
        Z = atomic_numbers[symbol]
        if symbol not in radius:
            # find the rare gas of the elements row
            Zrg = None
            for Zr in Z_rare_gas:
                if Zrg is None and Z <= Zr:
                    Zrg = Zr

            n_g, r_g = get_density(Zrg)
            # find density at R
            R = collected_vdWradii[chemical_symbols[Zrg]] / Bohr
            n = 0
            while r_g[n] < R:
                n += 1
            # linear interpolation
            ncut = (n_g[n - 1] +
                    (n_g[n] - n_g[n - 1]) * (R - r_g[n - 1]) / (r_g[n] - r_g[n - 1]))
#            print "Z, Zrg, ncut", Z, Zrg, ncut

            # find own R at this density
            n_g, r_g = get_density(Z)
            n = 0
            while n_g[n] > ncut:
                n += 1
            # linear interpolation
            R = (r_g[n - 1] +
                 (r_g[n] - r_g[n - 1]) * (ncut - n_g[n - 1]) / (n_g[n] - n_g[n - 1]))
            radius[symbol] = R * Bohr

        radii.append(radius[symbol])

    return radii
