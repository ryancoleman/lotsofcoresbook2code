"""creates: bigpicture.svg bigpicture.png"""

import os
from math import pi, cos, sin

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


class Box:
    def __init__(self, name, description=(), attributes=(), color='grey'):
        self.name = name
        if isinstance(description, str):
            description = [description]
        self.description = description
        self.attributes = attributes
        self.color = color

        self.owns = []
        self.position = None
            
    def set_position(self, position):
        self.position = np.asarray(position)

    def has(self, other, name, angle=None, distance=None, x=0.4, style='<-'):
        self.owns.append((other, name, x, style))
        if angle is not None:
            angle *= pi / 180
            other.set_position(self.position +
                               [cos(angle) * distance, sin(angle) * distance])


def cut(size, dir):
    if abs(size[0] * dir[1]) < abs(size[1] * dir[0]):
        x = min(max(-size[0] / 2, dir[0]), size[0] / 2)
        y = x * dir[1] / dir[0]
    else:
        y = min(max(-size[1] / 2, dir[1]), size[1] / 2)
        x = y * dir[0] / dir[1]
    return x, y


class MPL:
    def __init__(self, boxes):
        self.boxes = boxes

    def plot(self):
        a4 = 100 * np.array([2**-1.75, 2**-2.25])
        inch = 2.54
        
        self.fig = plt.figure(1, a4 / inch)
        self.ax = ax = self.fig.add_axes([0, 0, 1, 1], frameon=False)
        ax.set_xlim(0, a4[0])
        ax.set_ylim(0, a4[1])
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        ax.add_patch(mpatches.Rectangle((22.5, 16), 6, 4, fc='orange'))
        ax.text(22.7, 19.5, 'ASE package')

        for b in boxes:
            x, y = b.position

            text = b.name
            for txt in b.description:
                text += '\n' + txt
            for txt in b.attributes:
                text += '\n' + txt

            b.text = ax.text(x, y,
                             text,
                             fontsize=9,
                             ha='center',
                             va='center',
                             bbox=dict(boxstyle='round',
                                       facecolor=b.color,
                                       alpha=0.75))

        self.fig.canvas.mpl_connect('draw_event', self.on_draw)
        plt.savefig('bigpicture.png', dpi=50)
        plt.savefig('bigpicture.svg')
        os.system('cp bigpicture.svg ../_build')

    def on_draw(self, event):
        for b in self.boxes:
            bbox = b.text.get_window_extent()
            t = b.text.get_transform()
            b.size = t.inverted().transform(bbox.size)
        
        for b in self.boxes:
            for other, name, s, style in b.owns:
                d = other.position - b.position
                p1 = b.position + cut(b.size, d)
                p2 = other.position + cut(other.size, -d)
                if style == '-|>':
                    arrowprops = dict(arrowstyle=style, fc='white')
                else:
                    arrowprops = dict(arrowstyle=style)

                self.ax.annotate('', p1, p2,
                                 arrowprops=arrowprops)
                if name:
                    p = (1 - s) * p1 + s * p2
                    self.ax.text(p[0], p[1], name, fontsize=7,
                                 ha='center', va='center',
                                 bbox=dict(facecolor='white', ec='white'))

        self.fig.canvas.callbacks.callbacks[event.name] = {}
        self.fig.canvas.draw()
        return False


boxes = []

def box(*args, **kwargs):
    b = Box(*args, **kwargs)
    boxes.append(b)
    return b

atoms = box('Atoms', [''], ['positions, numbers, cell, pbc'],
            color='white')
paw = box('PAW', [], [], 'green')
scf = box('SCFLoop', [])
density = box('Density', 
              [r'$\tilde{n}_\sigma = \sum_{\mathbf{k}n}' +
               r'|\tilde{\psi}_{\sigma\mathbf{k}n}|^2' +
               r'+\frac{1}{2}\sum_a \tilde{n}_c^a$',
               r'$\tilde{\rho}(\mathbf{r}) = ' +
               r'\sum_\sigma\tilde{n}_\sigma + \sum_{aL}Q_L^a \hat{g}_L^a$'],
              ['nspins, nt_sG, nt_sg,', 'rhot_g, Q_aL, D_asp'])
mixer = box('Mixer')#, color='blue')
hamiltonian = box('Hamiltonian',
                  [r'$-\frac{1}{2}\nabla^2 + \tilde{v} + ' +
                   r'\sum_a \sum_{i_1i_2} |\tilde{p}_{i_1}^a \rangle ' +
                   r'\Delta H_{i_1i_2} \langle \tilde{p}_{i_2}^a|$'],
                  ['nspins, vt_sG, vt_sg, vHt_g, dH_asp',
                   'Etot, Ekin, Exc, Epot, Ebar'])
wfs = box('WaveFunctions',
          [r'$\tilde{\psi}_{\sigma\mathbf{k}n}(\mathbf{r})$'],
          ['nspins, ibzk_qc, mynbands',
           'kpt_comm, band_comm'], color='magenta')
gd = box('GridDescriptor', ['(coarse grid)'],
         ['cell_cv, N_c,', 'pbc_c, dv, comm'], 'orange')
finegd = box('GridDescriptor', '(fine grid)',
         ['cell_cv, N_c, pbc_c, dv, comm'], 'orange')
rgd = box('RadialGridDescriptor', [], ['r_g, dr_g, rcut'], color='orange')
setups = box('Setups', ['', '', '', ''], ['nvalence, nao, Eref, corecharge'])
xccorrection = box('XCCorrection')
nct = box('LFC', r'$\tilde{n}_c^a(r)$', [], 'red')
vbar = box('LFC', r'$\bar{v}^a(r)$', [], 'red')
ghat = box('LFC', r'$\hat{g}_{\ell m}^a(\mathbf{r})$', [], 'red')
fd = box('FDWaveFunctions',
         r"""$\tilde{\psi}_{\sigma\mathbf{k}n}(ih,jh,kh)$""",
         [], 'magenta')
pt = box('LFC', r'$\tilde{p}_i^a(\mathbf{r})$', [], 'red')
lcao = box('LCAOWaveFunctions',
           r"$\tilde{\psi}_{\sigma\mathbf{k}n}(\mathbf{r})=\sum_{\mu\mathbf{R}} C_{\sigma\mathbf{k}n\mu} \Phi_\mu(\mathbf{r} - \mathbf{R}) \exp(i\mathbf{k}\cdot\mathbf{R})$",
           ['S_qMM, T_qMM, P_aqMi'], 'magenta')
atoms0 = box('Atoms', '(copy)', ['positions, numbers, cell, pbc'],
             color='grey')
parameters = box('InputParameters', [], ['xc, nbands, ...'])
forces = box('ForceCalculator')
occupations = box(
    'OccupationNumbers',
    r'$\epsilon_{\sigma\mathbf{k}n} \rightarrow f_{\sigma\mathbf{k}n}$')
poisson = box('PoissonSolver',
              r'$\nabla^2 \tilde{v}_H(\mathbf{r}) = -4\pi \tilde{\rho}(\mathbf{r})$')
eigensolver = box('EigenSolver')
symmetry = box('Symmetry')
restrictor = box('Transformer', '(fine -> coarse)',
                 color='yellow')
interpolator = box('Transformer', '(coarse -> fine)',
                   color='yellow')
xc = box('XCFunctional')
kin = box('FDOperator', r'$-\frac{1}{2}\nabla^2$')
hsoperator = box('HSOperator',
                 [r"$\langle \psi_n | A | \psi_{n'} \rangle$",
                  r"$\sum_{n'}U_{nn'}|\tilde{\psi}_{n'}\rangle$"])
                 
overlap = box('Overlap')
basisfunctions = box('BasisFunctions', r'$\Phi_\mu(\mathbf{r})$',
                     color='red')
tci = box('TwoCenterIntegrals',
          r'$\langle\Phi_\mu|\Phi_\nu\rangle,'
          r'\langle\Phi_\mu|\hat{T}|\Phi_\nu\rangle,'
          r'\langle\tilde{p}^a_i|\Phi_\mu\rangle$')

atoms.set_position((25, 18.3))
atoms.has(paw, 'calculator', -160, 7.5)
paw.has(scf, 'scf', 160, 4, x=0.48)
paw.has(density, 'density', -150, 14, 0.23)
paw.has(hamiltonian, 'hamiltonian', 180, 10, 0.3)
paw.has(wfs, 'wfs', -65, 5.5, x=0.48)
paw.has(atoms0, 'atoms', 9, 7.5)
paw.has(parameters, 'input_parameters', 90, 4)
paw.has(forces, 'forces', 50, 4)
paw.has(occupations, 'occupations', 136, 4)
density.has(mixer, 'mixer', 130, 3.3)
density.has(gd, 'gd', x=0.33)
density.has(finegd, 'finegd', 76, 3.5)
density.has(setups, 'setups', 0, 7, 0.45)
density.has(nct, 'nct', -90, 3)
density.has(ghat, 'ghat', -130, 3.4)
density.has(interpolator, 'interpolator', -45, 4)
hamiltonian.has(restrictor, 'restrictor', 40, 4)
hamiltonian.has(xc, 'xc', 160, 6, x=0.6)
hamiltonian.has(vbar, 'vbar', 80, 4)
hamiltonian.has(setups, 'setups', x=0.3)
hamiltonian.has(gd, 'gd', x=0.45)
hamiltonian.has(finegd, 'finegd')
hamiltonian.has(poisson, 'poissonsolver', 130, 4)
wfs.has(gd, 'gd', 160, 4.8, x=0.48)
wfs.has(setups, 'setups', x=0.4)
wfs.has(lcao, None, -55, 5.9, style='-|>')
wfs.has(fd, None, -112, 5.0, style='-|>')
wfs.has(eigensolver, 'eigensolver', 30, 5, x=0.6)
wfs.has(symmetry, 'symmetry', 80, 3)
fd.has(pt, 'pt', -45, 3.6)
fd.has(kin, 'kin', -90, 3)
fd.has(overlap, 'overlap', -135, 3.5)
lcao.has(basisfunctions, 'basis_functions', -50, 3.5)
lcao.has(tci, 'tci', -90, 4.2)
overlap.has(setups, 'setups', x=0.4)
overlap.has(hsoperator, 'operator', -115, 2.5, x=0.41)

for i in range(3):
    setup = box('Setup', [],
                ['Z, Nv, Nc, pt_j, nct,', 'vbar, ghat_l, Delta_pl'],
                'blue')
    setup.set_position(setups.position +
                       (0.9 - i * 0.14, 0.3 - i * 0.14))
setup.has(xccorrection, 'xc_correction', -110, 3.7)
xccorrection.has(rgd, 'rgd', -105, 2.4, 0.4)

kpts = [box('KPoint', [], ['psit_nG, C_nM,', 'eps_n, f_n, P_ani'],
            color='cyan') for i in range(3)]
wfs.has(kpts[1], 'kpt_u', 0, 5.4, 0.48)
kpts[0].set_position(kpts[1].position - 0.14)
kpts[2].set_position(kpts[1].position + 0.14)

MPL(boxes).plot()
