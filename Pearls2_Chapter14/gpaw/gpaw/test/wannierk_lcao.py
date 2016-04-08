'Test ase.dft.wannier module with k-points.'
from __future__ import print_function

from ase.lattice import bulk
from ase.dft.wannier import Wannier

from gpaw import GPAW
from gpaw.mpi import world, serial_comm

si = bulk('Si', 'diamond', a=5.43)

k = 2
if 1:
    si.calc = GPAW(kpts=(k, k, k), mode='lcao', txt='Si-ibz.txt')
    #si.calc = GPAW(kpts=(k, k, k), mode='lcao')
    e1 = si.get_potential_energy()
    si.calc.write('Si-ibz', mode='all')
    #si.calc.set(symmetry='off', txt='Si-bz.txt')
    #si.calc.set(symmetry='off')
    #e2 = si.get_potential_energy()
    #si.calc.write('Si-bz', mode='all')
    #print e1, e2

def wan(calc):
    centers = [([0.125, 0.125, 0.125], 0, 1.5),
               ([0.125, 0.625, 0.125], 0, 1.5),
               ([0.125, 0.125, 0.625], 0, 1.5),
               ([0.625, 0.125, 0.125], 0, 1.5)]
    w = Wannier(4, calc,
                nbands=4,
                verbose=False,
                initialwannier=centers)
    w.localize()
    x = w.get_functional_value()
    centers = (w.get_centers(1) * k) % 1
    c = (centers - 0.125) * 2
    #print w.get_radii()  # broken! XXX
    assert abs(c.round() - c).max() < 0.03
    c = c.round().astype(int).tolist()
    c.sort()
    assert c == [[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]]
    if 0:
        from ase.visualize import view
        from ase import Atoms
        watoms = calc.atoms + Atoms(symbols='X4',
                                    scaled_positions=centers,
                                    cell=calc.atoms.cell)
        view(watoms)
    return x

calc1 = GPAW('Si-bz.gpw', txt=None, communicator=serial_comm)
# do full lcao initialization
#print('Re-Init LCAO calc')
_s = calc1.get_atoms()
calc1.set_positions(_s)
calc1.wfs.eigensolver.iterate(calc1.hamiltonian, calc1.wfs)
x1 = wan(calc1)
    
#print x1
x2 = 5.76962685681      # from FD mode
assert abs(x1 - x2) < 0.001
#assert abs(x1 - 9.71) < 0.01

world.barrier()
