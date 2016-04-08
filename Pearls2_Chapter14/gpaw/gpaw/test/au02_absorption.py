import numpy as np
from ase import Atoms
from gpaw import GPAW, FermiDirac
from gpaw.response.df import DielectricFunction
from gpaw.test import equal, findpeak

GS = 1
ABS = 1
if GS:
    cluster = Atoms('Au2', [(0, 0, 0), (0, 0, 2.564)])
    cluster.set_cell((6, 6, 6), scale_atoms=False)
    cluster.center()
    calc = GPAW(mode='pw',
                dtype=complex,
                xc='RPBE',
                nbands=16,
                eigensolver='rmm-diis',
                occupations=FermiDirac(0.01))
    
    cluster.set_calculator(calc)
    cluster.get_potential_energy()
    calc.diagonalize_full_hamiltonian(nbands=24, scalapack=True)
    calc.write('Au2.gpw', 'all')

if ABS:
    df = DielectricFunction('Au2.gpw',
                            frequencies=np.linspace(0, 14, 141),
                            hilbert=not True,
                            eta=0.1,
                            ecut=10)

    b0, b = df.get_dielectric_function(filename=None,
                                       direction='z')
    a0, a = df.get_polarizability(filename=None,
                                  direction='z')
    df_ws = DielectricFunction('Au2.gpw',
                               frequencies=np.linspace(0, 14, 141),
                               hilbert=not True,
                               eta=0.1,
                               ecut=10,
                               truncation='wigner-seitz')
    
    a0_ws, a_ws = df_ws.get_polarizability(filename=None,
                                           direction='z')

    w0_ = 5.60491055
    I0_ = 244.693028
    w_ = 5.696528390
    I_ = 207.8
    
    w, I = findpeak(np.linspace(0, 14., 141), b0.imag)
    equal(w, w0_, 0.05)
    equal(6**3 * I / (4 * np.pi), I0_, 0.5)
    w, I = findpeak(np.linspace(0, 14., 141), a0.imag)
    equal(w, w0_, 0.05)
    equal(I, I0_, 0.5)
    w, I = findpeak(np.linspace(0, 14., 141), a0_ws.imag)
    equal(w, w0_, 0.05)
    equal(I, I0_, 0.5)
    w, I = findpeak(np.linspace(0, 14., 141), b.imag)
    equal(w, w_, 0.05)
    equal(6**3 * I / (4 * np.pi), I_, 0.5)
    w, I = findpeak(np.linspace(0, 14., 141), a.imag)
    equal(w, w_, 0.05)
    equal(I, I_, 0.5)
    # The Wigner-Seitz truncation does not give exactly the same for small cell
    w, I = findpeak(np.linspace(0, 14., 141), a_ws.imag)
    equal(w, w_, 0.2)
    equal(I, I_, 8.0)
