import numpy as np
# mathtext fails to create Angstrom with matplotlib 0.99 on el6
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
import ase.units as units
from ase.io import write
from gpaw import GPAW
from gpaw.poisson import PoissonSolver
from gpaw.dipole_correction import DipoleCorrection

# this test requires OpenEXR-libs

for name in ['zero', 'periodic', 'corrected']:
    if name == 'corrected':
        calc = GPAW(name, txt=None,
                    poissonsolver=DipoleCorrection(PoissonSolver(), 2))
    else:
        calc = GPAW(name, txt=None)

    efermi = calc.get_fermi_level()

    calc.restore_state()
    v = (calc.hamiltonian.vHt_g * units.Hartree).mean(0).mean(0)
    #v = (calc.hamiltonian.vt_sG[0] * units.Hartree).mean(0).mean(0)
    
    z = np.linspace(0, calc.atoms.cell[2, 2], len(v), endpoint=False)
    if not calc.atoms.pbc[2]:
        z += z[1]

    plt.figure(figsize=(6.5, 4.5))
    plt.plot(z, v, label='xy-averaged potential')
    plt.plot([0, z[-1]], [efermi, efermi], label='Fermi level')
    
    if name == 'corrected':
        plt.plot([0.2, 0.2], [efermi, v[0]], 'r:')
        plt.text(0.23, (efermi + v[0]) / 2,
                 '$\phi$ = %.2f eV' % (v[0] - efermi), va='center')
        plt.plot([z[-1] - 0.2, z[-1] - 0.2], [efermi, v[-1]], 'r:')
        plt.text(z[-1] - 0.23, (efermi + v[-1]) / 2,
                 '$\phi$ = %.2f eV' % (v[-1] - efermi), va='center', ha='right')
    
    plt.xlabel('$z$, $\AA$')
    plt.ylabel('(Pseudo) electrostatic potential, V')
    plt.xlim([0., z[-1]])
    plt.title(name.title() + ' boundary conditions')
    plt.savefig(name + '.png')

write('slab.pov', calc.atoms,
      rotation='-90x',
      show_unit_cell=2,
      transparent=False,
      display=False,
      run_povray=True)
