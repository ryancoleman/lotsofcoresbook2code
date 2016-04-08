from __future__ import print_function
from ase import Atoms
from ase.structure import molecule
from ase.parallel import paropen
from gpaw import GPAW, Mixer, MixerDif
from gpaw.utilities.tools import split_formula

cell = [14.4, 14.4, 14.4]
data = paropen('data.txt', 'a')

##Reference from J. Chem. Phys. Vol 120 No. 15, 15 April 2004, page 6898
tpss_de = [
('H2' , 112.9),
('LiH', 59.1),
('OH' , 106.8),
('HF' , 139.1),
('Li2', 22.5),
('LiF', 135.7),
('Be2', 8.1),
('CO' , 254.2),
('N2' , 227.7),
('O2' , 126.9),
('F2' , 46.4),
('P2' , 116.1),
('Cl2', 60.8)
]

exp_bonds_dE = [
('H2' , 0.741,109.5),
('LiH', 1.595,57.8),
('OH' , 0.970,106.4),
('HF' , 0.917,140.8),
('Li2', 2.673,24.4),
('LiF', 1.564,138.9),
('Be2', 2.440,3.0),
('CO' , 1.128,259.3),
('N2' , 1.098,228.5),
('O2' , 1.208,120.5),
('F2' , 1.412,38.5),
('P2' , 1.893,117.3),
('Cl2', 1.988,58.0)
]

systems = [ a[0] for a in tpss_de ]
ref = [ a[1] for a in tpss_de ]

# Add atoms
for formula in systems:
    temp = split_formula(formula)
    for atom in temp:
        if atom not in systems:
            systems.append(atom)
energies = {}

# Calculate energies
i = 0
for formula in systems:
    if formula == 'Be2':
        loa = Atoms('Be2', [(0, 0, 0), (0, 0, 2.0212)])
    else:
        loa = molecule(formula)
    loa.set_cell(cell)
    loa.center()
    width = 0.0
    calc = GPAW(h=.18,
                nbands=-5,
                maxiter=333,
                xc='PBE',
                txt=formula + '.txt')
    if len(loa) == 1:
        calc.set(hund=True)
        calc.set(fixmom=True)
        calc.set(mixer=MixerDif())
        calc.set(eigensolver='cg')
    else:
        calc.set(mixer=Mixer())
        pos = loa.get_positions()
        pos[1,:] = pos[0,:] + [exp_bonds_dE[i][1],0.0,0.0]
        loa.set_positions(pos)
        loa.center()
    loa.set_calculator(calc)
    try:
        energy = loa.get_potential_energy()
        difft = calc.get_xc_difference('TPSS')
        diffr = calc.get_xc_difference('revTPSS')
        diffm = calc.get_xc_difference('M06L')
        energies[formula]=(energy, energy+difft, energy+diffr,energy+diffm)
    except:
        print(formula, 'Error', file=data)
    else:
        print(formula, energy, energy+difft, energy+diffr, energy+diffm, file=data)
        data.flush()
    i += 1

#calculate atomization energies
ii =0
file = paropen('atom_en.dat', 'a')
print("# formula \t PBE \t TPSS \t revTPSS \t M06L \t Exp", file=file)
for formula in systems[:13]:
    try:
        atoms_formula = split_formula(formula)
        de_tpss = -1.0 * energies[formula][1]
        de_revtpss = -1.0 * energies[formula][2]
        de_m06l = -1.0 * energies[formula][3]
        de_pbe = -1.0 * energies[formula][0]
        for atom_formula in atoms_formula:
            de_tpss += energies[atom_formula][1]
            de_revtpss += energies[atom_formula][2]
            de_m06l += energies[atom_formula][3]
            de_pbe += energies[atom_formula][0]
    except:
        print(formula, 'Error', file=file)
    else:
        de_tpss *= 627.5/27.211
        de_revtpss *= 627.5/27.211
        de_m06l *= 627.5/27.211
        de_pbe *= 627.5/27.211
        out = "%s\t%.1f \t%.1f \t%.1f \t%.1f \t%.1f" %(formula, de_pbe, de_tpss, de_revtpss, de_m06l ,exp_bonds_dE[ii][2])
        print(out, file=file)
        file.flush()
    ii += 1
