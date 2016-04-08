from __future__ import print_function
import pickle
import numpy as np
from gpaw.mpi import world
from parms import V, H, h1, ion_shift
from kgf.GreenFunctions import NonEqNonIntGreenFunction, NonEqIntGreenFunction
from kgf.ScatteringHamiltonianMatrix import ScatteringHamiltonianMatrix
from kgf.LeadHamiltonianMatrix import LeadHamiltonianMatrix
from kgf.CurrentCalculator import CurrentCalculator
from kgf.selfenergies import NonEqConstantSelfEnergy
from kgf.selfenergies.GW2index import Hartree2index, Fock2index, GW2index

hmat = ScatteringHamiltonianMatrix(leftprincipallayer=1,
                                   rightprincipallayer=1,
                                   hamiltonianmatrix=H)

left = LeadHamiltonianMatrix(principallayer=1,
                             hamiltonianmatrix=h1)
    
right = LeadHamiltonianMatrix(principallayer=1,
                             hamiltonianmatrix=h1)

energies = np.arange(-96.0, 96.0, 0.02)
de = energies[1] - energies[0]
assert abs(energies).min() < 1.0e-8
fermi = 0.0
g0 = NonEqNonIntGreenFunction(hmat, left, right, E_Fermi=fermi,
                              energies=energies)
g0.SetInfinitesimal(2*de) 
Ntot = g0.GetTotalParticleNumber()
if world.rank==0: print('Total electron number:', Ntot)
se_null = NonEqConstantSelfEnergy(g0, np.zeros((6, 6), complex), 'null')
se_ion = NonEqConstantSelfEnergy(g0, ion_shift, 'shift')
se_hartree = Hartree2index(g0, V, initialhartree=np.zeros((6, 6), complex))
se_fock = Fock2index(g0, V)
se_gw = GW2index(g0, V)
gf = NonEqIntGreenFunction([se_null, se_ion, se_hartree, se_fock, se_gw])
orbitals = range(1, 5)
pulay = (0.03, 0.25, 1)
# Non interacting calculation
se_ion.SetActive(False)
se_hartree.SetActive(False)
se_fock.SetActive(False)
se_gw.SetActive(False)
gf.WriteSpectralInfoToNetCDFFile('nonint.nc', 
                                 diagonalize=True,
                                 orbitals=orbitals,
                                 spectral='individual')
# HF calculation
se_ion.SetActive(True)
se_hartree.SetActive(True)
se_fock.SetActive(True)
gf.SelfConsistent(log='HF.log', pulay=pulay)
gf.WriteSpectralInfoToNetCDFFile('HF.nc', 
                                 diagonalize=True,
                                 orbitals=orbitals,
                                 spectral='individual')

# GW calculation
se_gw.SetActive(True)
gf.SelfConsistent(log='GW.log', pulay=pulay)
gf.WriteSpectralInfoToNetCDFFile('GW.nc', 
                                 diagonalize=True,
                                 orbitals=orbitals,
                                 spectral='individual')


