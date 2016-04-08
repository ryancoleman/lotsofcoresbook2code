from ase import Atoms
from gpaw import GPAW
from gpaw.mpi import world, serial_comm

a = Atoms('H',
          cell=(1, 3, 3),
          pbc=1)

a.calc = GPAW(mode='pw',
              h=0.15,
              kpts=(4, 1, 1),
              basis='dzp',
              nbands=4,
              eigensolver='rmm-diis',
              txt=None)

a.get_potential_energy()
w1 = a.calc.get_pseudo_wave_function(0, 1)
e1 = a.calc.get_eigenvalues(1)

a.calc.write('H')

if world.size <= 2:
    scalapack = None
else:
    mb = world.size // 4
    scalapack = (2, world.size // 4, 32)

a.calc.diagonalize_full_hamiltonian(nbands=100, scalapack=scalapack)
w2 = a.calc.get_pseudo_wave_function(0, 1)
e2 = a.calc.get_eigenvalues(1)

calc = GPAW('H', txt=None)
calc.diagonalize_full_hamiltonian(nbands=100, scalapack=scalapack)
w3 = calc.get_pseudo_wave_function(0, 1)
e3 = calc.get_eigenvalues(1)

calc.write('Hwf', 'all')

calc = GPAW('Hwf', txt=None, communicator=serial_comm)
w4 = calc.get_pseudo_wave_function(0, 1)
e4 = calc.get_eigenvalues(1)

for w in [w2, w3, w4]:
    assert abs(abs(w[1, 2, 3]) - abs(w1[1, 2, 3])) < 1e-7

for e in [e2, e3, e4]:
    assert abs(e[0] - e1[0]) < 2e-9, abs(e[0] - e1[0])
    assert abs(e[-1] - e2[-1]) < 1e-10
