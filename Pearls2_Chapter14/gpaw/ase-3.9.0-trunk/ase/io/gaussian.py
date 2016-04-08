"""
Read/write functions for Gaussian.
Written by:

   Glen R. Jenness
   University of Wisconsin - Madison

See accompanying license files for details.
"""

import numpy as np

import ase.units
from ase.atoms import Atoms
from ase.atom import Atom
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.gaussian_reader import GaussianReader as GR
from ase.calculators.gaussian import Gaussian

# http://www.gaussian.com/g_tech/g_ur/k_dft.htm
allowed_dft_functionals = ['lsda',  # = 'svwn'
                           'svwn',
                           'svwn5',  # != 'svwn'
                           'blyp',
                           'b3lyp',
                           'bp86',
                           'pbepbe',
                           'pbe1pbe',  # pbe0
                           'm06',
                           'm06hf',
                           'm062x',
                           'tpssh',
                           'tpsstpss',
                           'wb97xd']


def read_gaussian_out(filename, index=-1, quantity='atoms'):
    """"Interface to GaussianReader and returns various quantities"""
    energy = 0.0

    data = GR(filename)[index]

    formula = data['Chemical_formula']
    positions = np.array(data['Positions'])
    method = data['Method']
    version = data['Version']
    charge = data['Charge']
    multiplicity = data['Multiplicity']

    if method.lower()[1:] in allowed_dft_functionals:
        method = 'HF'

    atoms = Atoms(formula, positions=positions)

    for key, value in data.items():
        if (key in method):
            energy = value

    try:
# Re-read in the log file
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        forces = list()
        for n, line in enumerate(lines):
            if ('Forces (Hartrees/Bohr)' in line):
                for j in range(len(atoms)):
                    forces += [[float(lines[n + j + 3].split()[2]),
                                float(lines[n + j + 3].split()[3]),
                                float(lines[n + j + 3].split()[4])]]
        convert = ase.units.Hartree / ase.units.Bohr
        forces = np.array(forces) * convert
    except:
        forces = None

    energy *= ase.units.Hartree  # Convert the energy from a.u. to eV
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    atoms.set_calculator(calc)

    if (quantity == 'energy'):
        return energy
    elif (quantity == 'forces'):
        return forces
    elif (quantity == 'dipole'):
        return np.array(data['Dipole'])
    elif (quantity == 'atoms'):
        return atoms
    elif (quantity == 'version'):
        return version
    elif (quantity == 'multiplicity'):
        return multiplicity
    elif (quantity == 'charge'):
        return charge


def read_gaussian(filename):
    """Reads a Gaussian input file"""
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    atoms = Atoms()
    for n, line in enumerate(lines):
        if ('#' in line):
            i = 0
            while (lines[n + i + 5] != '\n'):
                info = lines[n + i + 5].split()
                symbol = info[0]
                position = [float(info[1]), float(info[2]), float(info[3])]
                atoms += Atom(symbol, position=position)
                i += 1
    return atoms


def write_gaussian(filename, atoms):
    """Writes a basic Gaussian input file"""
# Since Gaussian prints the geometry directly into the input file, we'll just
# the write_input method from the Gaussian calculator, and just use the
# default settings
    calc = Gaussian()
    calc.initialize(atoms)
    calc.write_input(filename, atoms)
