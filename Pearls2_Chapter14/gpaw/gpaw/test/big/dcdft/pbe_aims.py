import os
import sys
import time

import numpy as np

import ase.db
from ase.utils import opencew
from ase.calculators.calculator import kpts2mp
from ase.io.trajectory import PickleTrajectory
from ase.calculators.aims import Aims
from ase.test.tasks.dcdft import DeltaCodesDFTCollection as Collection

collection = Collection()

if len(sys.argv) == 1:
    names = collection.names
else:
    names = [sys.argv[1]]

c = ase.db.connect('dcdft_aims.db')

# select the basis set
basis = 'light'
#basis = 'tight'
#basis = 'really_tight'
#basis = 'tier2'

kptdensity = 16.0
width = 0.01

basis_threshold = 0.00001
relativistic = 'none'
relativistic = 1.e-12
relativistic = 'scalar'

linspace = (0.98, 1.02, 5)  # eos numpy's linspace
linspacestr = ''.join([str(t) + 'x' for t in linspace])[:-1]

code = 'aims' + '-' + basis + '_e' + linspacestr
code = code + '_k' + str(kptdensity) + '_w' + str(width)
code = code + '_t' + str(basis_threshold) + '_r' + str(relativistic)

collection = Collection()

for name in names:
    # save all steps in one traj file in addition to the database
    # we should only used the database c.reserve, but here
    # traj file is used as another lock ...
    fd = opencew(name + '_' + code + '.traj')
    if fd is None:
        continue
    traj = PickleTrajectory(name + '_' + code + '.traj', 'w')
    atoms = collection[name]
    cell = atoms.get_cell()
    kpts = tuple(kpts2mp(atoms, kptdensity, even=True))
    kwargs = {}
    if relativistic == 'scalar':
        kwargs.update({'relativistic': ['atomic_zora', relativistic]})
    elif relativistic == 'none':
        kwargs.update({'relativistic': 'none'})
    else:  # e.g. 1.0e-12
        kwargs.update({'relativistic': ['zora', relativistic]})
    if atoms.get_initial_magnetic_moments().any():  # spin-polarization
        magmom = atoms.get_initial_magnetic_moments().sum() / len(atoms)
        kwargs.update({'default_initial_moment': magmom, 'spin': 'collinear'})
    # loop over EOS linspace
    for n, x in enumerate(np.linspace(linspace[0], linspace[1], linspace[2])):
        id = c.reserve(name=name, basis=basis, linspacestr=linspacestr,
                       kptdensity=kptdensity, width=width,
                       basis_threshold=basis_threshold,
                       relativistic=relativistic,
                       x=x)
        if id is None:
            continue
        # perform EOS step
        atoms.set_cell(cell * x, scale_atoms=True)
        # set calculator
        atoms.calc = Aims(
            label=name + '_' + code + '_' + str(n),
            species_dir=os.path.join(os.environ['AIMS_SPECIES_DIR'], basis),
            xc='PBE',
            kpts=kpts,
            KS_method='elpa',
            sc_accuracy_rho=1.e-4,
            sc_accuracy_eev=5.e-3,
            occupation_type=['gaussian', width],
            override_relativity=True,
            override_illconditioning=True,
            basis_threshold=basis_threshold,
            charge_mix_param=0.01,
            )
        atoms.calc.set(**kwargs)  # remaining calc keywords
        t = time.time()
        atoms.get_potential_energy()
        c.write(atoms,
                name=name, basis=basis, linspacestr=linspacestr,
                kptdensity=kptdensity, width=width,
                basis_threshold=basis_threshold,
                relativistic=relativistic,
                x=x,
                time=time.time()-t)
        traj.write(atoms)
        del c[id]
