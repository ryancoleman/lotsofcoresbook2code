from __future__ import print_function
import numpy as np
from ase.transport.calculators import TransportCalculator
import pylab

# onsite energies 0.0, nearest neighbor hopping -1.0, and
# second nearest neighbor hopping 0.2
H_lead = np.array([[0., -1.,  0.2,  0.],
                   [-1.,  0., -1.,  0.2],
                   [0.2, -1.,  0., -1.],
                   [0.,  0.2, -1.,  0.]])

H_scat = np.zeros((6, 6))

# Principal layers on either side of S
H_scat[:2, :2] = H_scat[-2:, -2:] = H_lead[:2, :2]

# Scattering region (hydrogen molecule) - onsite 0.0 and hopping -0.8
H_scat[2:4, 2:4] = [[0.0, -0.8], [-0.8, 0.0]]

# coupling to the leads - nearest neighbor only
H_scat[1, 2] = H_scat[2, 1] = H_scat[3, 4] = H_scat[4, 3] = 0.2

tcalc = TransportCalculator(h=H_scat,  # Scattering Hamiltonian
                            h1=H_lead,  # Lead 1 (left)
                            h2=H_lead,  # Lead 2  (right)
                            energies=np.arange(-3, 3, 0.02))

T_e = tcalc.get_transmission()
pylab.plot(tcalc.energies, T_e)
pylab.title('Transmission function')
pylab.show()

tcalc.set(pdos=[2, 3])
pdos_ne = tcalc.get_pdos()
pylab.plot(tcalc.energies, pdos_ne[0], ':')
pylab.plot(tcalc.energies, pdos_ne[1], '--')
pylab.title('Projected density of states')
pylab.show()

h_rot, s_rot, eps_n, vec_nn = tcalc.subdiagonalize_bfs([2, 3])
tcalc.set(h=h_rot, s=s_rot)  # Set the rotated matrices
for n in range(2):
    print("eigenvalue, eigenvector:", eps_n[n], ',', vec_nn[:, n])

pdos_rot_ne = tcalc.get_pdos()
pylab.plot(tcalc.energies, pdos_rot_ne[0], ':')
pylab.plot(tcalc.energies, pdos_rot_ne[1], '--')
pylab.title('Projected density of states (rotated)')
pylab.show()

h_cut, s_cut = tcalc.cutcoupling_bfs([2])
tcalc.set(h=h_cut, s=s_cut)
T_cut_bonding_e = tcalc.get_transmission()
pylab.plot(tcalc.energies, T_cut_bonding_e)
pylab.title('Transmission (bonding orbital cut)')
pylab.show()
