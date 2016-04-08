from gpaw.upf import UPFSetupData, upfplot
from gpaw.pseudopotential import pseudoplot, PseudoPotential
from gpaw.eigensolvers.davidson import Davidson

def get(fname):
    s = UPFSetupData(fname)
    return s

s = get('O.pz-hgh.UPF')
#upfplot(s.data, show=True)
bfs = s.create_basis_functions()
pp = PseudoPotential(s, bfs)

import numpy as np
x = np.linspace(0.0, 5.0, 1000)
dr = x[1] - x[0]

psi = pp.phit_j[0].map(x)
p = pp.pt_j[0].map(x)

from gpaw.atom.atompaw import AtomPAW

if 0:
    f = 1.0 #1e-12
    c = AtomPAW('H', [[[2.0], [4.0]]],
                #charge=1 - f,
                h=0.04,
                #setups='paw',
                #setups='hgh',
                setups={'H': s}
                )

    import pylab as pl
    pl.plot(c.wfs.gd.r_g, c.hamiltonian.vt_sg[0])
    pl.show()
    raise SystemExit

#print 'test v201 Au.pz-d-hgh.UPF'
#s = UPFSetupData('Au.pz-d-hgh.UPF')
#print 'v201 ok'

#print 'test horrible version O.pz-mt.UPF'
#s = UPFSetupData('O.pz-mt.UPF')
#print 'horrible version ok, relatively speaking'

if 1:
    from ase import Atoms
    from gpaw import GPAW, PoissonSolver, FermiDirac, Mixer, MixerSum
    from gpaw.utilities import h2gpts
    from ase.structure import molecule

    #s = UPFSetupData('/home/askhl/parse-upf/h_lda_v1.uspp.F.UPF')

    upfsetups = {'H': UPFSetupData('H.pz-hgh.UPF'),
                 'O': UPFSetupData('O.pz-hgh.UPF')}

    system = molecule('H2O')
    system.center(vacuum=3.5)
    calc = GPAW(txt='-',
                nbands=6,
                setups=upfsetups,
                #setups='paw',
                #hund=True,
                #mixer=MixerSum(0.1, 5, 20.0),
                #eigensolver='cg',
                #occupations=FermiDirac(0.1),
                #charge=1-1e-12,
                eigensolver=Davidson(2),
                #eigensolver='rmm-diis',
                gpts=h2gpts(0.12, system.get_cell(), idiv=8),
                poissonsolver=PoissonSolver(relax='GS', eps=1e-7),
                xc='oldLDA',
                #nbands=4
                )

    system.set_calculator(calc)
    system.get_potential_energy()
