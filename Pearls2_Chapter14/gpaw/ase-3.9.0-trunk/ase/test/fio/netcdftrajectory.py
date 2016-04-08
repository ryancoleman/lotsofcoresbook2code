from ase.test import NotAvailable

import ase.io.netcdftrajectory as netcdftrajectory

if not netcdftrajectory.have_nc:
    raise NotAvailable('No NetCDF module available (netCDF4-python, '
                       'scipy.io.netcdf or ase.io.pupynere)')

import os
from ase import Atom, Atoms
from ase.io import NetCDFTrajectory

co = Atoms([Atom('C', (0, 0, 0)),
            Atom('O', (0, 0, 1.2))], pbc=True)
traj = NetCDFTrajectory('1.nc', 'w', co)
for i in range(5):
    co.positions[:, 2] += 0.1
    traj.write()
del traj
if netcdftrajectory.have_nc == netcdftrajectory.NC_IS_NETCDF4:
    traj = NetCDFTrajectory('1.nc', 'a')
    co = traj[-1]
    #print co.positions
    co.positions[:] += 1
    traj.write(co)
    del traj
    t = NetCDFTrajectory('1.nc', 'a')
else:
    t = NetCDFTrajectory('1.nc', 'r')
#print t[-1].positions
#print '.--------'
for i, a in enumerate(t):
    if i < 4:
        #print 1, a.positions[-1, 2], 1.3 + i * 0.1
        assert abs(a.positions[-1, 2] - 1.3 - i * 0.1) < 1e-6
        assert a.pbc.all()
    else:
        #print 1, a.positions[-1, 2], 1.7 + i - 4
        assert abs(a.positions[-1, 2] - 1.7 - i + 4) < 1e-6
        assert a.pbc.all()
if netcdftrajectory.have_nc == netcdftrajectory.NC_IS_NETCDF4:
    co.positions[:] += 1
    t.write(co)
    for i, a in enumerate(t):
        if i < 4:
            #print 2, a.positions[-1, 2], 1.3 + i * 0.1
            assert abs(a.positions[-1, 2] - 1.3 - i * 0.1) < 1e-6
        else:
            #print 2, a.positions[-1, 2], 1.7 + i - 4
            assert abs(a.positions[-1, 2] - 1.7 - i + 4) < 1e-6
    assert len(t) == 7
else:
    assert len(t) == 5

co[0].number = 1
try:
    t.write(co)
except ValueError:
    pass
else:
    assert False

if netcdftrajectory.have_nc == netcdftrajectory.NC_IS_NETCDF4:
    co[0].number = 6
    co.pbc = True
    t.write(co)

    co.pbc = False
    o = co.pop(1)
    try:
        t.write(co)
    except ValueError:
        pass
    else:
        assert False

    co.append(o)
    co.pbc = True
    t.write(co)
del t

# append to a nonexisting file
if netcdftrajectory.have_nc == netcdftrajectory.NC_IS_NETCDF4:
    fname = '2.nc'
    if os.path.isfile(fname):
        os.remove(fname)
    t = NetCDFTrajectory(fname, 'a', co)
    del t

fname = '3.nc'
t = NetCDFTrajectory(fname, 'w', co)
# File is not created before first write
co.set_pbc([True, False, False])
d = co.get_distance(0, 1)
t.write(co)
del t
# Check pbc
t = NetCDFTrajectory(fname)
a = t[-1]
assert a.pbc[0] and not a.pbc[1] and not a.pbc[2]
assert abs(a.get_distance(0, 1) - d) < 1e-6
del t 
os.remove(fname)
