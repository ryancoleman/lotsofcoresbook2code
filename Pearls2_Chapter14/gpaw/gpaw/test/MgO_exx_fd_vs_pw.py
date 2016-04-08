from __future__ import print_function
from ase.dft.kpoints import monkhorst_pack
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW


nk = 4
for mode in ('pw', 'fd'):
    kpts = monkhorst_pack((nk, nk, nk))
    kshift = 1.0 / (2 * nk)
    kpts += kshift

    atoms = bulk('MgO', 'rocksalt', a=4.212)

    if mode == 'pw':
        calc = GPAW(mode=PW(800), basis='dzp', xc='PBE', maxiter=300,
                    kpts=kpts, parallel={'band': 1, 'domain': 1},
                    setups={'Mg': '2'},
                    occupations=FermiDirac(0.01))
        atoms.set_calculator(calc)
        E1 = atoms.get_potential_energy()
        from gpaw.xc.hybridg import HybridXC
        exx = HybridXC('EXX', method='acdf')
        E_hf1 = E1 + calc.get_xc_difference(exx)
        
    else:
        calc = GPAW(h=0.12,
                    basis='dzp', kpts=kpts, xc='PBE',
                    setups={'Mg': '2'},
                    parallel={'domain': 1, 'band': 1},
                    occupations=FermiDirac(0.01))
        atoms.set_calculator(calc)
        E2 = atoms.get_potential_energy()
        from gpaw.xc.hybridk import HybridXC
        exx = HybridXC('EXX', acdf=True)
        E_hf2 = E2 + calc.get_xc_difference(exx)

print((E1, E2, E_hf1, E_hf2))
assert abs(E1 - E2) < 0.05
assert abs(E_hf1 - E_hf2) < 0.05
