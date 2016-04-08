# -*- coding: utf-8 -*-
# creates: atomization_energies.csv bondlengths.csv bondlengths.png

from ase.atoms import string2symbols
from ase.data.g2_1 import data
from ase.data.g2_1_ref import atomization_vasp, diatomic
from ase.data.molecules import latex
from ase.structure import molecule
from ase.units import kcal, mol
import numpy as np
import pylab as plt


dimers = diatomic.keys()
dimers.remove('FH')
molecules = atomization_vasp.keys()
atoms = set()
for m in molecules:
    atoms.update(molecule(m).get_chemical_symbols())
systems = molecules + list(atoms)

E = {'NH3': -19.889, 'S2': -7.082, 'SiH2_s3B1d': -8.765, 'CH3OH': -30.696,
     'SiH4': -18.877, 'Si2H6': -30.888, 'PH3': -15.567, 'PH2': -10.792,
     'HF': -8.706, 'O2': -10.598, 'SiH3': -13.816, 'NH': -8.361,
     'SH2': -11.166, 'ClO': -6.119, 'H2O2': -18.884, 'NO': -13.042,
     'ClF': -4.948, 'LiH': -3.741, 'HCO': -17.574, 'CH3': -18.262,
     'CH4': -24.157, 'Cl2': -3.609, 'HOCl': -11.314, 'SiH2_s1A1d': -9.483,
     'SiO': -11.503, 'F2': -5.172, 'P2': -8.988, 'Si2': -5.217, 'CH': -6.239,
     'CO': -15.281, 'CN': -13.384, 'LiF': -7.701, 'Na2': -1.194,
     'SO2': -17.548, 'NaCl': -4.699, 'Li2': -1.445, 'NH2': -13.831,
     'CS': -10.285, 'C2H6': -40.737, 'N2': -17.382, 'C2H4': -32.205,
     'HCN': -20.159, 'C2H2': -23.174, 'CH2_s3B1d': -12.125, 'CH3Cl': -22.544,
     'BeH': -3.520, 'CO2': -23.886, 'CH3SH': -27.720, 'OH': -8.089,
     'N2H4': -31.003, 'H2O': -14.579, 'SO': -9.356, 'CH2_s1A1d': -11.451,
     'H2CO': -22.638, 'HCl': -6.110, 'H': -1.120, 'Li': -0.290, 'Be': -0.001,
     'B': -0.470, 'C': -1.444, 'N': -3.405, 'O': -2.220, 'F': -1.430,
     'Na': -0.214, 'Mg': 0.003, 'Al': -0.306, 'Si': -0.855, 'P': -1.878,
     'S': -1.051, 'Cl': -0.384}

dE = [('BeH', (0.0231, 0.0057, -0.0005, 0.0032, 0.0157)),
      ('ClO', (0.0810, 0.0182, -0.0064, 0.0004, 0.0335)),
      ('CO', (0.1673, 0.0454, -0.0025, 0.0122, 0.0800)),
      ('CN', (0.0438, -0.0510, -0.0785, -0.0491, 0.0298)),
      ('ClF', (0.0743, 0.0183, -0.0002, 0.0091, 0.0427)),
      ('LiH', (0.0147, 0.0026, -0.0023, -0.0006, 0.0068)),
      ('F2', (0.0623, 0.0144, -0.0004, 0.0122, 0.0504)),
      ('LiF', (0.0556, 0.0133, -0.0096, -0.0137, -0.0013)),
      ('Na2', (0.0046, -0.0013, -0.0027, -0.0002, 0.0058)),
      ('CH', (0.0254, 0.0036, -0.0028, 0.0042, 0.0225)),
      ('HCl', (0.0470, 0.0112, -0.0012, 0.0071, 0.0326)),
      ('Li2', (0.0081, 0.0010, -0.0012, 0.0008, 0.0066)),
      ('N2', (0.0857, -0.0260, -0.0566, -0.0173, 0.0816)),
      ('O2', (0.1818, 0.0625, 0.0046, -0.0027, 0.0323)),
      ('Cl2', (0.0826, 0.0228, -0.0000, 0.0088, 0.0430))]


def atomization_energies(E):
    """Write given atomization energies to file atomization_energies.csv."""
    Ea = {}
    fd = open('atomization_energies.csv', 'w')
    for formula in sorted(molecules):
        ea = -E[formula]
        for a in string2symbols(data[formula]['symbols']):
            ea += E[a]
        eavasp = atomization_vasp[formula][1] * kcal / mol
        Ea[formula] = (ea, eavasp)
        name = latex(data[formula]['name'])
        fd.write('`%s`, %.3f, %.3f, %+.3f\n' %
                 (name[1:-1], ea, eavasp, ea - eavasp))

    return Ea


def bondlengths(Ea, dE):
    """Calculate bond lengths and write to bondlengths.csv file"""
    B = []
    E0 = []
    csv = open('bondlengths.csv', 'w')
    for formula, energies in dE:
        bref = diatomic[formula][1]
        b = np.linspace(0.96 * bref, 1.04 * bref, 5)
        e = np.polyfit(b, energies, 3)
        if not formula in Ea:
            continue
        ea, eavasp = Ea[formula]
        dedb = np.polyder(e, 1)
        b0 = np.roots(dedb)[1]
        assert abs(b0 - bref) < 0.1
        b = np.linspace(0.96 * bref, 1.04 * bref, 20)
        e = np.polyval(e, b) - ea
        if formula == 'O2':
            plt.plot(b, e, '-', color='0.7', label='GPAW')
        else:
            plt.plot(b, e, '-', color='0.7', label='_nolegend_')
        name = latex(data[formula]['name'])
        plt.text(b[0], e[0] + 0.2, name)
        B.append(bref)
        E0.append(-eavasp)
        csv.write('`%s`, %.3f, %.3f, %+.3f\n' %
                  (name[1:-1], b0, bref, b0 - bref))
        
    plt.plot(B, E0, 'g.', label='reference')
    plt.legend(loc='lower right')
    plt.xlabel('Bond length $\mathrm{\AA}$')
    plt.ylabel('Energy [eV]')
    plt.savefig('bondlengths.png')


Ea = atomization_energies(E)
bondlengths(Ea, dE)
