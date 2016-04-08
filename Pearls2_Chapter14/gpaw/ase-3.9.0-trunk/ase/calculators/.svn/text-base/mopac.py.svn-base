"""
Version 2012/08/20, Torsten Kerber

Contributors:
  Torsten Kerber, Ecole normale superieure de Lyon:
  Paul Fleurat-Lessard, Ecole normale superieure de Lyon
  based on a script by Rosa Bulo, Ecole normale superieure de Lyon

This work is supported by Award No. UK-C0017, made by King Abdullah
University of Science and Technology (KAUST), Saudi Arabia

See accompanying license files for details.
"""
import os
import numpy as np

from ase.units import kcal, mol
from ase.calculators.general import Calculator

str_keys = ['functional', 'job_type', 'command']
int_keys = ['restart', 'spin']
bool_keys = ['OPT']
float_keys = ['RELSCF']


class Mopac(Calculator):
    name = 'MOPAC'
    def __init__(self, label='ase', **kwargs):
        # define parameter fields
        self.str_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.float_params = {}
        
        # initials parameter fields
        for key in str_keys:
            self.str_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in float_keys:
            self.float_params[key] = None
                        
        # set initial values
        self.set(restart=0,
                 spin=0,
                 OPT=False,
                 functional='PM6',
                 job_type='NOANCI 1SCF GRADIENTS AUX(0,PRECISION=9)',
                 RELSCF=0.0001)
        # set user values
        self.set(**kwargs)

        # save label
        self.label = label
        
        #set atoms
        self.atoms = None
        # initialize the results
        self.version = None
        self.energy_zero = None
        self.energy_free = None
        self.forces = None
        self.stress = None
        
        # initialize the results
        self.occupations = None
        
    def set(self, **kwargs):
        """
        Sets the parameters on the according keywords
        Raises RuntimeError when wrong keyword is provided
        """
        for key in kwargs:
            if key in self.bool_params:
                self.bool_params[key] = kwargs[key]
            elif key in self.int_params:
                self.int_params[key] = kwargs[key]
            elif key in self.str_params:
                self.str_params[key] = kwargs[key]
            elif key in self.float_params:
                self.float_params[key] = kwargs[key]
            else:
                raise RuntimeError('MOPAC calculator: unknown keyword: ' + key)

    def get_version(self):
        return self.version

    def initialize(self, atoms):
        pass

    def write_input(self, fname, atoms):
        """
        Writes the files that have to be written each timestep
        """
        
        # start the input
        mopac_input = ''

        #write functional and job_type
        for key in 'functional', 'job_type':
            if self.str_params[key] != None:
                mopac_input += self.str_params[key] + ' '
                
        if self.float_params['RELSCF'] != None:
            mopac_input += 'RELSCF=' + str(self.float_params['RELSCF']) + ' '
            
        #write charge
        charge = sum(atoms.get_initial_charges())
        if charge != 0:
            mopac_input += 'CHARGE=%i ' % (charge)
        
        #write spin
        spin = self.int_params['spin']
        if spin == 1.:
            mopac_input += 'DOUBLET '
        elif spin == 2.:
            mopac_input += 'TRIPLET '

        #input down
        mopac_input += '\n'
        mopac_input += 'Title: ASE job\n\n'

        f = 1
        # write coordinates
        for iat in xrange(len(atoms)):
            atom = atoms[iat]
            xyz = atom.position
            mopac_input += ' %2s' % atom.symbol
            # write x, y, z
            for idir in xrange(3):
                mopac_input += '    %16.5f %i' % (xyz[idir], f)
            mopac_input += '\n'

        if atoms.pbc.any():
            for v in atoms.get_cell():
                mopac_input += 'Tv %8.3f %8.3f %8.3f\n' % (v[0], v[1], v[2])
        
        # write input
        myfile = open(fname, 'w')
        myfile.write(mopac_input)
        myfile.close()

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if self.str_params['command'] is not None:
            command = self.str_params['command']
        elif ('MOPAC_COMMAND' in os.environ):
            command = os.environ['MOPAC_COMMAND']
        return command

    def run(self):
        """
        Writes input in label.mop
        Runs MOPAC
        Reads Version, Energy and Forces
        """
        # set the input file name
        finput = self.label + '.mop'
        foutput = self.label + '.out'
        
        self.write_input(finput, self.atoms)

        command = self.get_command()
        if command is None:
            raise RuntimeError('MOPAC command not specified')
        
        exitcode = os.system('%s %s' % (command, finput))
        
        if exitcode != 0:
            raise RuntimeError('MOPAC exited with error code')

        self.version = self.read_version(foutput)

        energy = self.read_energy(foutput)
        self.energy_zero = energy
        self.energy_free = energy
        
        self.forces = self.read_forces(foutput)

    def read_version(self, fname):
        """
        Reads the MOPAC version string from the second line
        """
        version = 'unknown'
        lines = open(fname).readlines()
        for line in lines:
            if "  Version" in line:
                version = line.split()[-2]
                break
        return version

    def read_energy(self, fname):
        """
        Reads the ENERGY from the output file (HEAT of FORMATION in kcal / mol)
        Raises RuntimeError if no energy was found
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()

        energy = None
        for line in lines:
            if line.find('HEAT OF FORMATION') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('H.o.F. per unit cell') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('UNABLE TO ACHIEVE SELF-CONSISTENCE') != -1:
                energy = None
        if energy is None:
            raise RuntimeError('MOPAC: could not find total energy')
        
        energy *= (kcal / mol)
        return energy

    def read_forces(self, fname):
        """
        Reads the FORCES from the output file
        search string: (HEAT of FORMATION in kcal / mol / AA)
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()

        nats = len(self.atoms)
        forces = np.zeros((nats, 3), float)
        
        for i, line in enumerate(lines):
            if line.find('GRADIENT\n') != -1:
                for j in range(nats * 3):
                    gline = lines[i + j + 1]
                    forces[j / 3, j % 3] = float(gline[49:62])
                break
        
        forces *= - (kcal / mol)
        return forces
        
    def atoms_are_equal(self, atoms_new):
        ''' (adopted from jacapo.py)
        comparison of atoms to self.atoms using tolerances to account
        for float/double differences and float math.
        '''
    
        TOL = 1.0e-6  # angstroms

        # check for change in cell parameters
        test = len(atoms_new) == len(self.atoms)
        if test is not True:
            return False
        
        # check for change in cell parameters
        test = (abs(self.atoms.get_cell() - atoms_new.get_cell()) <= TOL).all()
        if test is not True:
            return False
        
        old = self.atoms.arrays
        new = atoms_new.arrays
        
        # check for change in atom position
        test = (abs(new['positions'] - old['positions']) <= TOL).all()
        if test is not True:
            return False
        
        # passed all tests
        return True

    def update(self, atoms_new):
        if not self.atoms_are_equal(atoms_new):
            self.atoms = atoms_new.copy()
            self.run()
