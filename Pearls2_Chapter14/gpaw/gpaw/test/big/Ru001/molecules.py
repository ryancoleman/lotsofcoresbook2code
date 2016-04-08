from ase.structure import molecule
from ase.optimize import BFGS
from gpaw import GPAW
from gpaw.mixer import MixerDif

for name in ['H2', 'N2', 'O2', 'NO']:
    mol = molecule(name)
    mol.center(vacuum=5.0)
    calc = GPAW(xc='PBE',
                h=0.2,
                txt=name + '.txt')
    if name == 'NO':
        mol.translate((0, 0.1, 0))
        calc.set(mixer=MixerDif(0.05,5))
    mol.set_calculator(calc)
  
    opt = BFGS(mol, logfile=name + '.log', trajectory=name + '.traj')
    opt.run(fmax=0.05)
    calc.write(name)
