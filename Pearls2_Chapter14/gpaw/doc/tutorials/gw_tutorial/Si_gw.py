import numpy as np
from gpaw.response.gw import GW

gw = GW(
        file='Si_groundstate.gpw',
        nbands=100,                # number of bands for calculation of self-energy
        bands=np.array([2,3,4,5]), # here: two highest valence and two lowest conduction bands
        kpoints=None,              # by default: all points in irreducible Brillouin zone
        ecut=50.,                  # plane wave cutoff for self-energy
        ppa=True,                  # use Plasmon Pole Approximation
        txt='Si-3k_GW.out'
       )

gw.get_exact_exchange()

gw.get_QP_spectrum()
