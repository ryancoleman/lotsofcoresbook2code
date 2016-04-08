from __future__ import print_function

from gpaw import *
from ase import *
import numpy as np

atom = Atoms(symbols='KTaO3',
             pbc=np.array([ True,  True,  True], dtype=bool),
             cell=np.array(
                           [[ 4.,  0.,  0.],
                            [ 0.,  4.,  0.],
                            [ 0.,  0.,  4.]]
                           ),
             positions=np.array(
                            [[ 0.,  0.,  0.],
                             [ 2.,  2.,  2.],
                             [ 2.,  2.,  0.],
                             [ 0.,  2.,  2.],
                             [ 2.,  0.,  2.]]
                                ),
            )


calc = GPAW(h=0.16,
            kpts=(10,10,10),
            xc='GLLBSC',
            txt='KTaO3.out',
            occupations=FermiDirac(width=0.05),
            )

atom.set_calculator(calc)
atom.get_potential_energy()

#Important commands for calculating the response and the
#derivatice discontinuity
response = calc.hamiltonian.xc.xcs['RESPONSE']
response.calculate_delta_xc()
EKs, Dxc = response.calculate_delta_xc_perturbation()

# fundamental band gap
# EKs = kohn-sham bandgap
# Dxc = derivative discontinuity
Gap = EKs+Dxc

print("Calculated band gap:", Gap)
