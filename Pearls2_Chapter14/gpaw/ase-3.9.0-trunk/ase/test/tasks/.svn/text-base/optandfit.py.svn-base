import numpy as np
from ase.tasks.main import run

# molecule

# fit
atoms, task = run('molecule H2 -F 5,2 --atomize -t fitfail')
try:
    atoms, task = run('molecule H2 H -t fitfail -s')
except ValueError:
    pass
# fitting outside of range must fail!
assert task.data == {}

# fit in range
# when only fitting the number of points not must be odd
# in this case data['energy'] is the energy of the middle fit point
#
# test trailing space
atoms, task = run('molecule H2 -F 5,7 --atomize -t fit ')
atoms, task = run('molecule H2 H -t fit -s')
data = task.data['H2']
assert abs(data['energy'] - 1.1589) < 0.0001
assert abs(data['relaxed energy'] - 1.0705) < 0.0001
# note slightly different bondlenght from fitting
assert abs(data['distance'] - 0.77900) < 0.00001
assert abs(data['frequency'] - 0.8676) < 0.0001
assert abs(data['atomic energy'] - data['relaxed energy'] - 5.3495) < 0.0001

# opt then fit
# when fitting after optimization the number of points not need to be odd
# in this case data['energy'] is the original energy before optimization
#
# test leading space
atoms, task = run(' molecule H2 -R 0.001,FIRE -F 6,2 --atomize -t optfit')
atoms, task = run('molecule H2 H -t optfit -s')
data = task.data['H2']
assert abs(data['energy'] - 1.1589) < 0.0001
assert abs(data['relaxed energy'] - 1.0705) < 0.0001
# note slightly different bondlength from the fitting above!
assert abs(data['distance'] - 0.77905) < 0.00001
assert abs(data['frequency'] - 0.8628) < 0.0001
assert abs(data['atomic energy'] - data['relaxed energy'] - 5.3495) < 0.0001

# opt
atoms, task = run('molecule H2 -R 0.001,BFGSLineSearch --atomize -t opt')
atoms, task = run('molecule H2 H -t opt -s')
data = task.data['H2']
assert data['optimizer steps'] == 4
assert data['optimizer force calls'] == 5
assert abs(data['relaxed energy'] - 1.0705) < 0.0001
assert abs(data['distance'] - 0.77905) < 0.00001
assert abs(data['atomic energy'] - data['relaxed energy'] - 5.3495) < 0.0001

# optimization with unsufficient number of steps (ASE does not fail, simpy stops)
atoms, task = run('molecule H2 -R 0.001 --relaxsteps 1 --atomize -t opt')
atoms, task = run('molecule H2 H -t opt -s')
data = task.data['H2']
assert data['optimizer steps'] == 1
assert abs(data['relaxed energy'] - 1.1294) < 0.0001
assert abs(data['distance'] - 0.81717) < 0.00001
assert abs(data['atomic energy'] - data['relaxed energy'] - 5.2906) < 0.0001

# bulk

# fit (the system slightly distorted)
# when only fitting the number of points not must be odd
# in this case data['energy'] is the energy of the middle fit point
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -F 5,-0.5 --modify=system.positions[0,2]+=0.1 -t fit')
# non-default lattice constant must be passed to the analysis part
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -t fit -s')
data = task.data['NiO']
assert abs(data['fitted energy'] - 1.1455) < 0.0001
assert abs(data['volume'] - 20.2594) < 0.0001
assert abs(data['B'] - 0.9317) < 0.0001

# fit sensitivity to sampling (same initial structure)
#
# test mid space
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -F 5,1  --modify=system.positions[0,2]+=0.1 -t fit')
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -t fit -s')
data = task.data['NiO']
assert abs(data['fitted energy'] - 1.1455) < 0.0001
assert abs(data['volume'] - 20.2595) < 0.0001
assert abs(data['B'] - 0.9303) < 0.0001

# fit sensitivity to equation of state (same data)
try:
    import scipy
    atoms, task = run('bulk NiO -x rocksalt -a 4.32 --eos murnaghan -t fit -s')
    data = task.data['NiO']
    assert abs(data['fitted energy'] - 1.1455) < 0.0001
    assert abs(data['volume'] - 20.2595) < 0.0001
    assert abs(data['B'] - 0.9301) < 0.0001
except ImportError:
    pass

# opt and fit (same initial structure)
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -R 0.01,BFGS -F 5,-0.5 --modify=system.positions[0,2]+=0.1 -t optfit5')
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -t optfit5 -s')
data = task.data['NiO']
assert data['optimizer force calls'] == data['optimizer steps'] == 3
assert abs(data['energy'] - 1.1458) < 0.0001
assert abs(data['fitted energy'] - 1.1254) < 0.0001
assert abs(data['volume'] - 20.1513) < 0.0001
assert abs(data['B'] - 0.9407) < 0.0001

# opt and fit (different initial structure)
# when fitting after optimization the number of points not need to be odd
# in this case data['energy'] is the original energy before optimization
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -R 0.01 -F 6,-0.5 -t optfit6')
atoms, task = run('bulk NiO -x rocksalt -a 4.32 -t optfit6 -s')
data = task.data['NiO']
assert abs(data['energy'] - 1.1254) < 0.0001
assert abs(data['fitted energy'] - 1.1254) < 0.0001
assert abs(data['volume'] - 20.1513) < 0.0001
assert abs(data['B'] - 0.9407) < 0.0001
