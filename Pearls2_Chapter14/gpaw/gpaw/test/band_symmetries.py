from __future__ import print_function
from ase.lattice import bulk
from gpaw import GPAW
from gpaw import restart
from _gpaw import symmetrize
import numpy as np
from gpaw.symmetry import Symmetry
from gpaw.utilities.tools import md5_array
from gpaw.test import equal
import os

# Normal ground state calculation
atoms = bulk('Si')
calc = GPAW(h=0.24, kpts=(4,4,4), 
            convergence={'eigenstates' : 1e-4, 'density' : 1e-3},
    )
atoms.set_calculator(calc)
atoms.get_potential_energy()

# "Bandstructure" calculation (only Gamma point here)
kpts = np.array(((0, 0, 0),))
calc.set(fixdensity=True, kpts=kpts)
atoms.get_potential_energy()
calc.write('Si_gamma.gpw', mode='all')

# Analyse symmetries of wave functions
atoms, calc = restart('Si_gamma.gpw', txt=None)

# Find symmetries (in Gamma-point calculations calculator does not
# use symmetry)
sym = Symmetry(calc.wfs.setups.id_a, atoms.cell, atoms.pbc)
sym.analyze(atoms.get_scaled_positions())

def find_classes(op_all_scc):
    # Find classes of group represented by matrices op_all_scc
    # and return representative operations
    op_scc = [op_all_scc[0]]
    for op1 in op_all_scc[1:]:
        new_class = True
        for op2 in op_all_scc:
            op_tmp = (np.dot(np.dot(op2, op1), np.linalg.inv(op2).astype(int)))
            # Check whether operation op1 belongs to existing class 
            for op in op_scc:
                if np.all((op_tmp - op) == 0):
                    new_class = False
                    break
            if not new_class:
                break
        if new_class:
            op_scc.append(op1)
    return op_scc

op_scc = find_classes(sym.op_scc)
nops = len(op_scc)
nbands = calc.wfs.bd.nbands

characters = np.zeros((nbands, nops))

degeneracy_threshold = 1e-2 # Threshold for considering bands degenerate
symdone = np.zeros(nbands, bool)

eig_n = calc.get_eigenvalues()
for n in range(nbands):
    # For degenerate states analysis might be already done
    if symdone[n]:
        continue
    # Find degenerate bands
    ndeg = 0
    m = n
    while abs(eig_n[m] - eig_n[n]) < degeneracy_threshold:
        ndeg += 1
        m += 1
        if m == nbands:
            break

    # Representation matrix
    representation_nn = np.zeros((ndeg, ndeg, nops))
    for nop, op_cc in enumerate(op_scc):
        # Bands from n to m are (possibly) degenerate
        for n1 in range(n, m):
            for n2 in range(n, m):
                wf1 = calc.get_pseudo_wave_function(band=n1)
                wf2 = calc.get_pseudo_wave_function(band=n2)
                norm1 = np.sqrt(calc.wfs.gd.integrate(wf1, wf1))
                norm2 = np.sqrt(calc.wfs.gd.integrate(wf2, wf2).real)
                wf_rot = np.zeros_like(wf2)
                symmetrize(wf2, wf_rot, op_cc)
                # Indices of representation matrix are from 0 to ndeg
                i1, i2 = n1 - n, n2 - n
                representation_nn[i1, i2, nop] = calc.wfs.gd.integrate(wf1,
                                                                       wf_rot) 
                representation_nn[i1, i2, nop] /= norm1 * norm2

    # Calculate traces of irreducible representations
    # If bands i1 and i2 are accidentally degenerate (i.e. not due to symmetry)
    # they belong to different irreducible representations and the
    # corresponding representation matrix elements are zero  for all 
    # symmetry operations. 
    for i1 in range(ndeg):
        for i2 in range(ndeg):
            if np.any(abs(representation_nn[i1, i2, :]) > 0.01):
                characters[n + i1, :] += representation_nn[i2, i2, :]
        symdone[n + i1] = True

# Use four decimals for characters
characters = np.round(characters, 4)
# Use characters as fingerprints
fingerprints = np.array([md5_array(row) for row in characters])

fmt = "%6.4f " * nops
for i in range(nbands):
     print(fmt % tuple([c for c in characters[i, :nops]]))

# Correct?!? character table
characters_reference = np.array(((1.0, 1.0, 1.0, 1.0, 1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (3.0, 1.0, 0.0, -1.0, -1.0),
                                 (1.0, 1.0, 1.0, 1.0, 1.0)))
assert np.all(np.abs(characters_reference-characters) < 1.0e-4)
os.remove('Si_gamma.gpw')
