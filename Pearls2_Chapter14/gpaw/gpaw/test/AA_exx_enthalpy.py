from __future__ import print_function
from ase import Atoms, Atom
from ase.structure import molecule
from ase.parallel import barrier
from ase.units import Hartree, mol, kcal
from gpaw import GPAW, setup_paths
from gpaw.mixer import Mixer, MixerSum
from gpaw.occupations import FermiDirac
from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw.test import equal, gen
from gpaw.mpi import rank

from os import remove
from os.path import exists

data = {}

# data (from tables.pdf of 10.1063/1.1626543)
data['N'] = { # intermolecular distance (A),
              #formation enthalpy(298) (kcal/mol) on B3LYP geometry
    'exp': (1.098, 0.0, 'none', 'none'),
    'HCTH407': (1.097, 7.9, 'HCTH407', 'gga'),
    'PBE': (1.103, -15.1, 'PBE', 'gga'),
    'BLYP': (1.103, -11.7, 'BLYP', 'gga'),
    'BP86': (1.104, -15.6, 'BP86', 'gga'),
    'BPW91': (1.103, -8.5, 'BPW91', 'gga'),
    'B3LYP': (1.092, -1.03, 'BLYP', 'hyb_gga'),
    'B3PW91': (1.091, 2.8, 'PW91', 'hyb_gga'),
    'PBE0': (1.090, 3.1, 'PBE', 'hyb_gga'),
    'PBEH': (1.090, 3.1, 'PBE', 'hyb_gga'),
    'magmom': 3.0,
    # tables.pdf:
    # http://ftp.aip.org/epaps/journ_chem_phys/E-JCPSA6-119-302348/tables.pdf
    'R_AA_B3LYP': 1.092, # (from tables.pdf of 10.1063/1.1626543) (Angstom)
    'ZPE_AA_B3LYP': 0.005457 * Hartree, # (from benchmarks.txt of
                                        # 10.1063/1.1626543) (eV)
    'H_298_H_0_AA_B3LYP': 0.003304 * Hartree, # (from benchmarks.txt of
                                              # 10.1063/1.1626543) (eV)
    'H_298_H_0_A': 1.04 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    'dHf_0_A': 112.53 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    }

data['O'] = { # intermolecular distance (A),
              # formation enthalpy(298) (kcal/mol) on B3LYP geometry
    'exp': (1.208, 0.0, 'none', 'none'),
    'HCTH407': (1.202, -14.5, 'HCTH407', 'gga'),
    'PBE': (1.218, -23.6, 'PBE', 'gga'),
    'BLYP': (1.229, -15.4, 'BLYP', 'gga'),
    'BP86': (1.220, -21.9, 'BP86', 'gga'),
    'BPW91': (1.219, -17.9, 'BPW91', 'gga'),
    'B3LYP': (1.204, -3.7, 'BLYP', 'hyb_gga'),
    'B3PW91': (1.197, -5.1, 'PW91', 'hyb_gga'),
    'PBE0': (1.192, -4.3, 'PBE', 'hyb_gga'),
    'PBEH': (1.192, -4.3, 'PBE', 'hyb_gga'),
    'magmom': 2.0,
    # tables.pdf:
    # http://ftp.aip.org/epaps/journ_chem_phys/E-JCPSA6-119-302348/tables.pdf
    'R_AA_B3LYP': 1.204, # (from tables.pdf of 10.1063/1.1626543) (Angstom)
    'ZPE_AA_B3LYP': 0.003736 * Hartree, # (from benchmarks.txt of
                                        # 10.1063/1.1626543) (eV)
    'H_298_H_0_AA_B3LYP': 0.003307*Hartree, # (from benchmarks.txt of
                                            # 10.1063/1.1626543) (eV)
    'H_298_H_0_A': 1.04 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    'dHf_0_A': 58.99 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    }

data['H'] = { # intermolecular distance (A),
              # formation enthalpy(298) (kcal/mol) on B3LYP geometry
    'exp': (0.741, 0.0, 'none', 'none'),
    'HCTH407': (0.744, 1.8, 'HCTH407', 'gga'),
    'PBE': (0.750, 5.1, 'PBE', 'gga'),
    'BLYP': (0.746, 0.3, 'BLYP', 'gga'),
    'BP86': (0.750, -1.8, 'BP86', 'gga'),
    'BPW91': (0.748, 4.0, 'BPW91', 'gga'),
    'B3LYP': (0.742, -0.5, 'BLYP', 'hyb_gga'),
    'B3PW91': (0.744, 2.4, 'PW91', 'hyb_gga'),
    'PBE0': (0.745, 5.3, 'PBE', 'hyb_gga'),
    'PBEH': (0.745, 5.3, 'PBE', 'hyb_gga'),
    'magmom': 1.0,
    # tables.pdf:
    # http://ftp.aip.org/epaps/journ_chem_phys/E-JCPSA6-119-302348/tables.pdf
    'R_AA_B3LYP': 0.742, # (from tables.pdf of 10.1063/1.1626543) (Angstom)
    'ZPE_AA_B3LYP': 0.010025 * Hartree, # (from benchmarks.txt of
                                        # 10.1063/1.1626543) (eV)
    'H_298_H_0_AA_B3LYP': 0.003305 * Hartree, # (from benchmarks.txt of
                                              # 10.1063/1.1626543) (eV)
    'H_298_H_0_A': 1.01 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    'dHf_0_A': 51.63 / (mol / kcal), # (from 10.1063/1.473182) (eV)
    }

def calculate(element, h, vacuum, xc, magmom):

    atom = Atoms([Atom(element, (0, 0, 0))])
    if magmom > 0.0:
        mms = [magmom for i in range(len(atom))]
        atom.set_initial_magnetic_moments(mms)

    atom.center(vacuum=vacuum)

    mixer = MixerSum(beta=0.2)
    if element == 'O':
        mixer = MixerSum(nmaxold=1, weight=100)
        atom.set_positions(atom.get_positions()+[0.0, 0.0, 0.0001])

    calc_atom = GPAW(h=h, xc=data[element][xc][2],
                     eigensolver='rmm-diis',
                     occupations=FermiDirac(0.0, fixmagmom=True),
                     mixer=mixer,
                     nbands=-2,
                     txt='%s.%s.txt' % (element, xc))
    atom.set_calculator(calc_atom)

    mixer = Mixer(beta=0.2, weight=100)
    compound = molecule(element+'2')
    if compound == 'O2':
        mixer = MixerSum(beta=0.2)
        mms = [1.0 for i in range(len(compound))]
        compound.set_initial_magnetic_moments(mms)

    calc = GPAW(h=h, xc=data[element][xc][2],
                eigensolver='rmm-diis',
                mixer=mixer,
                txt='%s2.%s.txt' % (element, xc))
    compound.set_distance(0,1, data[element]['R_AA_B3LYP'])
    compound.center(vacuum=vacuum)

    compound.set_calculator(calc)

    if data[element][xc][3] == 'hyb_gga': # only for hybrids
        e_atom = atom.get_potential_energy()
        e_compound = compound.get_potential_energy()

        calc_atom.set(xc=xc)
        calc.set(xc=xc)

    if 0:
        qn = QuasiNewton(compound)
        qn.attach(PickleTrajectory(
            element+'2'+'_'+xc+'.traj', 'w', compound).write)
        qn.run(fmax=0.02)
    e_atom = atom.get_potential_energy()
    e_compound = compound.get_potential_energy()

    dHf_0   = (e_compound - 2 * e_atom + data[element]['ZPE_AA_B3LYP'] +
               2 * data[element]['dHf_0_A'])
    dHf_298 = (dHf_0 + data[element]['H_298_H_0_AA_B3LYP'] -
               2 * data[element]['H_298_H_0_A']) * (mol / kcal)
    dist_compound = compound.get_distance(0,1)
    de = dHf_298-data[element][xc][1]
    E[element][xc] = de
    if rank == 0:
        print((xc, h, vacuum, dHf_298, data[element][xc][1], de,
               de/data[element][xc][1]))
        if element == 'H':
            equal(dHf_298, data[element][xc][1], 0.25, msg=xc+': ') # kcal/mol
        elif element == 'O':
            equal(dHf_298, data[element][xc][1], 7.5, msg=xc+': ') # kcal/mol
        else:
            equal(dHf_298, data[element][xc][1], 2.15, msg=xc+': ') # kcal/mol
        equal(de, E_ref[element][xc], 0.06, msg=xc+': ') # kcal/mol

E = {}

E_ref = {'H': {'HCTH407': 0.19286893273630645,
               'B3LYP': -0.11369634560501423,
               'PBE0': -0.21413764474738262,
               'PBEH': -0.14147808591211231},
         'N': {'HCTH407': 2.1354017840869268,
               'B3LYP': 0.63466589919873972,
               'PBE0': -0.33376468078480226,
               'PBEH': -0.30365500626180042}} # svnversion 5599 # -np 4

for element in ['H']:#, 'N']:#, 'O']: # oxygen atom fails to converge
    E[element] = {}
    for xc in ['HCTH407', 'PBE0', 'B3LYP']:
        setup = data[element][xc][2]
        if data[element][xc][3] == 'hyb_gga': # only for hybrids
            exx = True
        else:
            exx = False
        gen(element, exx=exx, xcname=setup)
        for h in [0.20]:
            for vacuum in [4.5]:
                calculate(element, h, vacuum, xc, data[element]['magmom'])
        barrier()
        if rank == 0:
            if exists(element+'.'+setup):
                remove(element+'.'+setup)
