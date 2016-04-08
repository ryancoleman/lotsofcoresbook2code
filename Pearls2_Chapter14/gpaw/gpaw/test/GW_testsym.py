from __future__ import print_function
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from ase.units import Hartree
import numpy as np
from gpaw.response.gw import GW

a = 5.431
atoms = bulk('Si', 'diamond', a=a)

kpts =(3,3,3)

calc = GPAW(
            h=0.24,
            kpts=kpts,
            xc='LDA',
            txt='Si_gs.txt',
            nbands=20,
            convergence={'bands':8},
            occupations=FermiDirac(0.001)
           )

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Si_kpt3.gpw','all')

    

file='Si_kpt3.gpw'

gw = GW(
        file=file,
        nbands=8,
        bands=np.array([2,3]),
        kpoints=np.arange(27),
        w=np.linspace(0., 60., 601),
        ecut=100.,
        eta=0.1
       )

gw.get_QP_spectrum()

gw.Qp_kn *= Hartree

checkQp_kn = np.zeros((gw.kd.nibzkpts, gw.gwnband))
nn = np.zeros(gw.kd.nibzkpts)
for k in range(gw.nkpt):
    ibzk = gw.kd.bz2ibz_k[k]
    checkQp_kn[ibzk,:] += gw.Qp_kn[k,:]
    nn[ibzk] += 1

for k in range(gw.kd.nibzkpts):
    checkQp_kn[k] /= nn[k]

for k in range(gw.nkpt):
    ibzk = gw.kd.bz2ibz_k[k]
    print(np.abs(checkQp_kn[ibzk] - gw.Qp_kn[k]))

# gw.Qp_kn
#[[ 6.19522578  6.19534951]
# [ 3.99972211  3.98923808]
# [ 1.27697908  4.28309876]
# [ 3.99980109  3.98910082]
# [ 6.18609008  6.16523424]
# [ 1.27796725  4.28057501]
# [ 1.27694066  4.28322734]
# [ 1.27793528  4.28067604]
# [ 1.27709399  4.27762828]
# [ 4.00034637  3.98848107]
# [ 6.19054241  6.16985403]
# [ 1.27781506  4.28614094]
# [ 6.19044301  6.16981151]
# [ 7.1159977   7.11541927]
# [ 6.19469536  6.17384384]
# [ 1.27780361  4.28614089]
# [ 6.19473472  6.17384953]
# [ 4.00680613  3.98247912]
# [ 1.27709219  4.27768906]
# [ 1.27792268  4.28066888]
# [ 1.27693712  4.28320036]
# [ 1.27796113  4.28054092]
# [ 6.19021475  6.16917482]
# [ 4.0073578   3.98186593]
# [ 1.27696901  4.28309883]
# [ 4.00743753  3.98172983]
# [ 6.19094467  6.19115549]]

