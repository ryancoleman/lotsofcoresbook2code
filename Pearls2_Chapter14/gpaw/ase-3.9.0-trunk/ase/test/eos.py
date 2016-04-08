# -*- coding: utf-8 -*-
# test ASE3 eos vs ASE2' on EMT Al bulk

import numpy as np

from ase.test import NotAvailable

try:
    import scipy
except ImportError:
    raise NotAvailable('This needs scipy module.')

from ase.lattice import bulk

from ase.io.trajectory import PickleTrajectory

from ase.calculators.emt import EMT

# old ASE2 conversion factor
eVA3ToGPA = 160.21773

ref = {
    'volumes': [29.205536, 30.581492, 32.000000, 33.461708, 34.967264],
    'energies': [0.0190898, -0.0031172, -0.0096925, -0.0004014, 0.0235753],
    # name: (V0 A**3, E0 eV, B eV/A**3)
    # ASE2: ScientificPython 2.6.2/Numeric 24.2
    'Taylor': (31.896496488942326, -0.0096090164907389405, 0.23802461480382878),
    'Murnaghan': (31.866877784374836, -0.0096119194044206324, 0.24202636566649313),
    'Birch': (31.866809942501359, -0.0096161509968013953, 0.24231157506701367),
    'BirchMurnaghan': (31.867394584147391, -0.009609309015137282, 0.23891301754324207),
    'PourierTarantola': (31.866473067615818, -0.009599545236557528, 0.24120474301680481),
    'Vinet': (31.866741599224699, -0.0096110298949974356, 0.24196956466978184),
    'AntonSchmidt': (31.745672779210317, 0.012772723347888704, 0.19905185689855259),
    # ASE3: scipy 0.7.0/numpy 1.3.0
    'sjeos': (31.867118229937798, -0.0096410046694188622, 0.23984474782755572),
    #
    'taylor': (31.867114798134253, -0.0096606904384420791, 0.24112293515031302),
    'murnaghan': (31.866729811658402, -0.0096340233039666941, 0.23937322901028654),
    'birch': (31.867567845123162, -0.0096525305272843597, 0.24062224387079953),
    'birchmurnaghan': (31.866655489091499, -0.0096461024146103497, 0.24016191590912808),
    'pouriertarantola': (31.866750629512403, -0.0096361387118443446, 0.23951298910150925),
    'vinet': (31.866655146818957, -0.0096368465365208426, 0.23955684756879458),
    'p3': (31.867115199307815, -0.0096606897797322233, 0.24112291100256208),
}

# original ASE2 methods

eos_strl = [
    'Taylor',
    'Murnaghan',
    'Birch',
    'BirchMurnaghan',
    'PourierTarantola',
    'Vinet',
    'AntonSchmidt',
    ]

# AntonSchmidt fails with ASE3!
# RuntimeError: Optimal parameters not found:
# Number of calls to function has reached maxfev = 1000.
eos_strl3 = [m for m in eos_strl]
eos_strl3.remove('AntonSchmidt')

results = {}

# prepare energies and volumes

b = bulk('Al', 'fcc', a=4.0, orthorhombic=True)
b.set_calculator(EMT())
cell = b.get_cell()

volumes = []
energies = []
traj = PickleTrajectory('eos.traj', 'w')
for x in np.linspace(0.97, 1.03, 5):
    b.set_cell(cell * x, scale_atoms=True)
    volumes.append(b.get_volume())
    energies.append(b.get_potential_energy())
    traj.write(b)

for n, (v, e) in enumerate(zip(volumes, energies)):
    vref = ref['volumes'][n]
    eref = ref['energies'][n]
    vabserr = abs((v - vref) / vref)
    vstrerr = str(n) + ' volume: ' + str(v) + ': ' + str(vref) + ': ' + str(vabserr)
    assert vabserr < 1.e-6, vstrerr
    eabserr = abs((e - eref) / eref)
    estrerr = str(n) + ' energy: ' + str(e) + ': ' + str(eref) + ': ' + str(eabserr)
    assert eabserr < 1.e-4, estrerr

# ASE2

try:
    from ASE.Utilities.EquationOfState import EquationOfState as eos2

    import sys

    from ase.utils import devnull

    for e in eos_strl:
        eos = eos2(e, volumes, energies)
        stdout_save = sys.stdout
        sys.stdout = devnull
        print eos
        results[e] = (eos.GetV0(), eos.GetEnergy(), eos.GetBulkModulus()/eVA3ToGPA)
        sys.stdout = stdout_save
except (ImportError, TypeError, ValueError):
    # TypeError, ValueError may appear with weird numeric/numpy/scientificpython
    pass

# ASE3

from ase.utils.eos import EquationOfState as eos3

for e in eos_strl3 + ['sjeos', 'p3']:
    eos = eos3(volumes, energies, eos=e.lower())
    v0, e0, B = eos.fit()
    results[e.lower()] = (v0, e0, B)

# test ASE2 vs ASE2 regression (if available)

for e in eos_strl:
    for n, v2 in enumerate(ref[e]):
        if n in [0, 2]: # only test volume and bulk modulus
            try:
                v3 = results[e][n]
                abserr = abs((v3 - v2) / v2)
                #print e, abserr
                strerr = e + ' 2 vs 2: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
                assert abserr < 1.e-6, strerr
            except KeyError:
                pass

# test ASE3 vs ASE3 regression

for e in eos_strl3:
    for n, v2 in enumerate(ref[e.lower()]):
        if n in [0, 2]: # only test volume and bulk modulus
            v3 = results[e.lower()][n]
            abserr = abs((v3 - v2) / v2)
            #print e, abserr
            strerr = e + ' 3 vs 3: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
            assert abserr < 5.e-6, strerr

# test ASE3 vs ASE2 reference

for e in eos_strl3:
    for n, v2 in enumerate(ref[e]):
        if n in [0, 2]: # only test volume and bulk modulus
            v3 = results[e.lower()][n]
            abserr = abs((v3 - v2) / v2)
            #print e, abserr
            strerr = e + ' 3 vs 2: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
            if n == 0: # volume
                assert abserr < 1.e-3, strerr
            else:
                # ASE2/ScientificPython/Numeric vs ASE2 methods/scipy/numpy error ~1% for bulk modulus!
                assert abserr < 2.e-2, strerr

# test ASE3: various eos between each other

for e1 in eos_strl3:
    for e2 in eos_strl3:
        for n, v2 in enumerate(ref[e1.lower()]):
            if n in [0, 2]: # only test volume and bulk modulus
                v3 = ref[e2.lower()][n]
                abserr = abs((v3 - v2) / v2)
                #print e, abserr
                strerr = '3: '+ e1 + ' vs ' + e2 + ' ' + str(v3) + ' : ' + str(v2) + ': ' + str(abserr)
                if n == 0: # volume
                    assert abserr < 5.e-5, strerr
                else:
                    # different eos disagree by ~1% for bulk modulus!
                    # The differences depend mostly on the sampling interval
                    # and less on the number of sampling points.
                    # Typical 5% of lattice constant is way too large for Al (~2% needed)!
                    assert abserr < 1.e-2, strerr
