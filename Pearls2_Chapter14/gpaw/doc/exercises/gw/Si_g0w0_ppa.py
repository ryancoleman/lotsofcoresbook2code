import numpy as np
from ase.parallel import parprint
from gpaw.response.g0w0 import G0W0

# We start by setting up a G0W0 calculator object
gw = G0W0('Si_gs.gpw',             # Path to groundstate gpw file
          filename='Si_g0w0_ppa',  # filename base for output files
          kpts=None,               # List of quasiparticle k-point indices
                                   # or None = all k-points
          bands=(3, 5),            # Range of quasiparticle bands - last
                                   # index is NOT included
          ecut=100.,               # Plane wave basis cut-off energy
          ppa=True)                # Use Plasmon Pole Approximation

# Perform the GW calculation. The results, ie. quasiparticle energies, as
# well as original Kohn-Sham eigenvalues, occupation numbers, DFT XC and
# self-energy contributions and renormalization factors are returned as a
# python dictionary object
result = gw.calculate()

ks_skn = result['eps']               # Get Kohn-Sham eigenvalues
ks_cbmin = np.amin(ks_skn[0, :, 1])  # DFT conduction band minimum
ks_vbmax = np.amax(ks_skn[0, :, 0])  # DFT valence band maximum
ks_gap = ks_cbmin - ks_vbmax         # DFT band gap

qp_skn = result['qp']                # GW quasiparticle energies
qp_cbmin = np.amin(qp_skn[0, :, 1])  # GW conduction band minimum
qp_vbmax = np.amax(qp_skn[0, :, 0])  # GW valence band maximum
qp_gap = qp_cbmin - qp_vbmax         # GW band gap

parprint('Kohn-Sham gap = %.3f' % ks_gap)
parprint('G0W0 gap = %.3f' % qp_gap)
