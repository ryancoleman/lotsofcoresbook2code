from __future__ import print_function
from ase.transport.calculators import TransportCalculator
import numpy as np
import pickle
import pylab

# Principal layer size
# Uncomment this line if going back to gpawtransport again
# pl = 4 * 9 # 9 is the number of bf per Pt atom (basis=szp), see below

# Read in the hamiltonians
h, s = pickle.load(file('scat_hs.pickle'))
# Uncomment this line if going back to gpawtransport again
# h, s = h[pl:-pl, pl:-pl], s[pl:-pl, pl:-pl]
h1, s1 = pickle.load(file('lead1_hs.pickle'))
h2, s2 = pickle.load(file('lead2_hs.pickle'))

tcalc = TransportCalculator(h=h, h1=h1, h2=h2,  # hamiltonian matrices
                            s=s, s1=s1, s2=s2,  # overlap matrices
                            align_bf=1)        # align the Fermi levels

# Calculate the conductance (the energy zero corresponds to the Fermi level)
tcalc.set(energies=[0.0])
G = tcalc.get_transmission()[0]
print('Conductance: %.2f 2e^2/h' % G)

# Determine the basis functions of the two Hydrogen atoms and subdiagonalize
Pt_N = 5    # Number of Pt atoms on each side in the scattering region
Pt_nbf = 9  # number of bf per Pt atom (basis=szp)
H_nbf = 4   # number of bf per H atom (basis=szp)
bf_H1 = Pt_nbf * Pt_N
bfs = range(bf_H1, bf_H1 + 2 * H_nbf)
h_rot, s_rot, eps_n, vec_jn = tcalc.subdiagonalize_bfs(bfs)
for n in range(len(eps_n)):
    print("bf %i corresponds to the eigenvalue %.2f eV" % (bfs[n], eps_n[n]))

# Switch to the rotated basis set
tcalc.set(h=h_rot, s=s_rot)

# plot the transmission function
tcalc.set(energies=np.arange(-8, 4, 0.1))
T = tcalc.get_transmission()
pylab.plot(tcalc.energies, T)
pylab.title('Transmission function')
pylab.show()

# ... and the projected density of states (pdos) of the H2 molecular orbitals
tcalc.set(pdos=bfs)
pdos_ne = tcalc.get_pdos()
pylab.plot(tcalc.energies, pdos_ne[0], label='bonding')
pylab.plot(tcalc.energies, pdos_ne[1], label='anti-bonding')
pylab.title('Projected density of states')
pylab.legend()
pylab.show()

# Cut the coupling to the anti-bonding orbital.
print('Cutting the coupling to the renormalized molecular state at %.2f eV' % (
    eps_n[1]))
h_rot_cut, s_rot_cut = tcalc.cutcoupling_bfs([bfs[1]])
tcalc.set(h=h_rot_cut, s=s_rot_cut)
pylab.plot(tcalc.energies, tcalc.get_transmission())
pylab.title('Transmission without anti-bonding orbital')
pylab.show()

# Cut the coupling to the bonding-orbital.
print('Cutting the coupling to the renormalized molecular state at %.2f eV' % (
    eps_n[0]))
tcalc.set(h=h_rot, s=s_rot)
h_rot_cut, s_rot_cut = tcalc.cutcoupling_bfs([bfs[0]])
tcalc.set(h=h_rot_cut, s=s_rot_cut)
pylab.plot(tcalc.energies, tcalc.get_transmission())
pylab.title('Transmission without bonding orbital')
pylab.show()
