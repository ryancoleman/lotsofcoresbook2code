from ase.structure import molecule
from ase.parallel import parprint

from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.analyse.hirshfeld import HirshfeldDensity, HirshfeldPartitioning
from gpaw.analyse.wignerseitz import WignerSeitz
from gpaw.test import equal
from gpaw.mpi import rank

h = 0.4
gpwname = 'H2O' + str(h) + '.gpw'

def run(lastres=[]):
    results = []

    # Hirshfeld ----------------------------------------

    if 1:

        hd = HirshfeldDensity(calc)

        # check for the number of electrons
        expected = [[None, 10],
                    [[0, 1, 2], 10],
                    [[1, 2], 2],
                    [[0], 8],
                    ]
        #expected = [[[0, 2], 9], ]
        #expected = [[None, 10], ]
        for gridrefinement in [1, 2, 4]:
            #Test for all gridrefinements for get_all_electron_density
            parprint('grid refinement', gridrefinement)
            for result in expected:
                indicees, result = result
                full, gd = hd.get_density(indicees, gridrefinement)
                parprint('indicees', indicees, end=': ') 
                parprint('result, expected:', gd.integrate(full), result)
                if gridrefinement < 4:
                    #The highest level of gridrefinement gets wrong electron numbers
                    equal(gd.integrate(full), result, 1.e-8)
                else:
                    equal(gd.integrate(full), result, 1.e-4)

        hp = HirshfeldPartitioning(calc)
        vr = hp.get_effective_volume_ratios()
        parprint('Hirshfeld:', vr)
        if len(lastres):
            equal(vr, lastres.pop(0), 1.e-10)
        results.append(vr)

    # Wigner-Seitz ----------------------------------------

    if 1:
        ws = WignerSeitz(calc.density.finegd, mol, calc)

        vr = ws.get_effective_volume_ratios()
        parprint('Wigner-Seitz:', vr)
        if len(lastres):
            equal(vr, lastres.pop(0), 1.e-10)
        results.append(vr)

    return results

mol = Cluster(molecule('H2O'))
mol.minimal_box(2.5, h=h)

# calculate
if 1:
    parprint('### fresh:')
    calc = GPAW(nbands=6,
                h = h, 
                txt=None)
if 1:
    calc.calculate(mol)
    calc.write(gpwname)
    lastres = run()

# load previous calculation
if 1:
    parprint('### reloaded:')
    calc = GPAW(gpwname, txt=None)
    mol = calc.get_atoms()
    run(lastres)

# periodic modulo test
parprint('### periodic:')
mol.set_pbc(True)
mol.translate(-mol[0].position)
mol.translate([-1.e-24, 0, 0])
calc.calculate(mol)
run()


