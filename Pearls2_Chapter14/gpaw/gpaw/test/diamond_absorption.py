import numpy as np
from ase.units import Bohr
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.response.df import DielectricFunction
from gpaw.test import equal, findpeak

a = 6.75 * Bohr
atoms = bulk('C', 'diamond', a=a)

calc = GPAW(mode='pw',
            kpts=(3, 3, 3),
            eigensolver='rmm-diis',
            occupations=FermiDirac(0.001))

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('C.gpw', 'all')

# Macroscopic dielectric constant calculation
df = DielectricFunction('C.gpw', frequencies=(0.,), eta=0.001, ecut=200,
                        hilbert=False)
eM1, eM2 = df.get_macroscopic_dielectric_constant()

eM1_ = 9.725
eM2_ = 9.068

equal(eM1, eM1_, 0.01)
equal(eM2, eM2_, 0.01)

# Absorption spectrum calculation
df = DielectricFunction('C.gpw', eta=0.25, ecut=200,
                        frequencies=np.linspace(0, 24., 241), hilbert=False)
b0, b = df.get_dielectric_function(filename=None)
df.check_sum_rule(b.imag)

equal(b0[0].real, eM1_, 0.01)
equal(b[0].real, eM2_, 0.01)

a0, a = df.get_polarizability(filename=None)

df_ws = DielectricFunction('C.gpw', eta=0.25, ecut=200,
                           frequencies=np.linspace(0, 24., 241), hilbert=False,
                           truncation='wigner-seitz')

a0_ws, a_ws = df_ws.get_polarizability(filename=None)

w0_ = 10.778232265664668
I0_ = 5.5467658790816268
w_ = 10.9530497246
I_ = 6.09704008088

w, I = findpeak(np.linspace(0, 24., 241), b0.imag)
equal(w, w0_, 0.05)
equal(I / (4 * np.pi), I0_, 0.05)
w, I = findpeak(np.linspace(0, 24., 241), a0.imag)
equal(w, w0_, 0.05)
equal(I, I0_, 0.05)
w, I = findpeak(np.linspace(0, 24., 241), a0_ws.imag)
equal(w, w0_, 0.05)
equal(I, I0_, 0.05)
w, I = findpeak(np.linspace(0, 24., 241), b.imag)
equal(w, w_, 0.05)
equal(I / (4 * np.pi), I_, 0.05)
w, I = findpeak(np.linspace(0, 24., 241), a.imag)
equal(w, w_, 0.05)
equal(I, I_, 0.05)
# The Wigner-Seitz truncation does not give exactly the same for kpts=(3,3,3)
w, I = findpeak(np.linspace(0, 24., 241), a_ws.imag)
equal(w, w_, 0.1)
equal(I, I_, 0.1)
