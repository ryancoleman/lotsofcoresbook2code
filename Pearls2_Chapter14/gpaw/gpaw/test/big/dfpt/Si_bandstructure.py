"""Phonon band-structure for silicon using the Appelbaum-Hamann PP.

For comparison, see e.g.:

    * Phys. Rev. B 43, 7231 (1991).

"""

import numpy as np
import pylab as plt

import ase.units as units
from ase.dft.kpoints import ibz_points, get_bandpath

from gpaw.mpi import rank, world
from gpaw.dfpt import PhononCalculator

# Pseudo-potential
PP = 'AH'

# Name of file with ground-state calculation
name = 'Si_%s.gpw' % PP

# Create phonon calculator
ph = PhononCalculator(name,
                      gamma=False,
                      symmetry=False,
                      e_ph=False)

# Run the self-consistent calculation
ph.run()

# Ensure that the master does not enter here before all files have been created
world.barrier()

# Calculate band-structure and plot on master
if rank == 0:

    # High-symmetry points in the Brillouin zone
    points = ibz_points['fcc']
    G = points['Gamma']
    X = points['X']
    W = points['W']
    K = points['K']
    L = points['L']

    atoms = ph.get_atoms()
    path_kc, q, Q = get_bandpath([G, K, X, G, L, X, W, L],
                                 atoms.cell, 100)
    point_names = ['$\Gamma$', 'K', 'X', '$\Gamma$', 'L', 'X', 'W', 'L']
    
    # Calculate band-structure
    omega_kn = ph.band_structure(path_kc)
    
    # Convert from sqrt(Ha / Bohr**2 / amu) -> meV
    s = units.Hartree**0.5 * units._hbar * 1.e10 / \
        (units._e * units._amu)**(0.5) / units.Bohr
    omega_kn *= s * 1000
    
    # Plot the band-structure
    plt.figure(1)
    for n in range(len(omega_kn[0])):
       plt.plot(q, omega_kn[:, n], 'k-', lw=2)
        
    plt.xticks(Q, point_names, fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(q[0], q[-1])
    plt.ylim(0, np.ceil(omega_kn.max() / 10) * 10)
    plt.ylabel("Frequency ($\mathrm{meV}$)", fontsize=22)
    plt.grid('on')
    # plt.show()
    plt.savefig('Si_bandstructure.png')
