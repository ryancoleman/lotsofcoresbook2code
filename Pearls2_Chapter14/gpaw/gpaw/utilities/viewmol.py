"""Trajectory module with viewmol support.""" 
from __future__ import print_function

from math import sqrt
import numpy as np

from ase import Angstrom, Atoms, Atom, Hartree, PickleTrajectory
from ase.io import filetype as ase_filetype
from ase.parallel import paropen
from ase.io.gpawtext import read_gpaw_text

from gpaw.cluster import Cluster
import gpaw.mpi as mpi
MASTER = 0

def filetype_from_ending(filename):
    # estimate file type from name ending
    ft = filename.split('.')[-1]
    return ft.lower()

class Trajectory(PickleTrajectory):
    """Basic trajectory class.

    In the moment, this class inherits from ase.PickleTrajectory.
    Might be other way round in some unknown future."""

    def __init__(self, filename=None, mode='r', atoms=None):
        self.images = []
        if filename is not None:
            self.read(filename)
        else:
            self.set_atoms(atoms)

    def __len__(self):
        if hasattr(self, 'offsets'):
            return PickleTrajectory.__len__(self)
        else:
            return len(self.images)

    def __del__(self):
        pass # overwrite PickleTrajectory.__del__
        
    def __getitem__(self, i=-1):
        if hasattr(self, 'offsets'):
            return PickleTrajectory.__getitem__(self, i)
        else:
            return self.images[i]

    def __iter__(self):
        if hasattr(self, 'offsets'):
            return PickleTrajectory.__iter__(self)
        else:
            return self.images.__iter__()

    def next(self):
        if hasattr(self, 'offsets'):
            return PickleTrajectory.next(self)
        else:
            return self.images.next()

    def append(self, images):
        if type(images) == type(Atoms):
            images = [images]
        for image in images:
            self.images.append(image)
        
    def read(self, filename, filetype=None):
        """Read from a file"""

        ftfe = filetype_from_ending(filename)
        aft = ase_filetype(filename)
        
        if ftfe == 'vmol':
            self.read_viewmol(filename)
        elif aft == 'traj':
            PickleTrajectory.__init__(self, filename)
        elif aft == 'gpaw-text':
            # this is strange code, but I don't know how to create the full
            # slice XXXXX
            last = read_gpaw_text(filename)
            self.__init__(atoms=last)
            s = slice(-1)
            rest = read_gpaw_text(filename, s)
            rest.append(last)
            self.append(rest)
        else:
            raise NotImplementedError('unknown file type "'+ aft +'"')

    def write(self, filename, filetype=None):
        """Write to a file"""

        if filetype is None:
            filetype = filetype_from_ending(filename)
        if filetype == 'vmol':
            self.write_viewmol(filename)
        else:
            raise NotImplementedError
    
    def write_viewmol(self, filename, mode='w'):
        """Write yourself in viewmol style."""

        f = paropen(filename, mode)
        if mode == 'w':
            print(' $coord', 1. / Angstrom, file=f)
            self.write_viewmol_positions(f, self[0])
            print(' $grad', file=f)
        else:
            raise NotImplementedError('append')

        for i, atoms in enumerate(self):
            self.write_viewmol_atoms(f, atoms, i)
            
        f.close()

    def write_viewmol_atoms(self, file, atoms, index):
        """Write current atomic position, energy and forces
        to the output file"""
        print('cycle=', index, end=' ', file=file)
        try:
            E = atoms.get_potential_energy() * Hartree
        except:
            E = 0
        print('SCF energy=', E, end=' ', file=file)
        try:
            forces = atoms.get_forces() * (Hartree / Angstrom)
            max_force = 0
            for f in forces:
                max_force = max(max_force, sqrt(np.sum(f * f)))
        except:
            forces = atoms.get_positions() * 0.
            max_force = 0
        print('|max dE/dxyz|=', max_force, file=file)
        self.write_viewmol_positions(file, atoms)
        for atom, f in zip(atoms, forces):
            print('%10.4g %10.4g %10.4g' % (f[0],f[1],f[2]), file=file)

    def write_viewmol_positions(self, file, atoms):
        for c, s in zip(atoms.get_positions(), atoms.get_chemical_symbols()):
            print('%10.4f %10.4f %10.4f' % (c[0],c[1],c[2]), s, file=file)

    def read_viewmol(self, filename):
        f = open(filename)

        # read the definition first
        definition=False
        for line in f:
            w = line.split()
            if not definition:
                if w[0] == '$coord':
                    definition=True
                    self.scale = float(w[1])
                    loa = Cluster([])
            else:
                if w[0] == '$grad':
                    # definition ends here
                    self.definition = loa
                    break
                else:
                    # we assume this is a coordinate entry
                    coo = (float(w[0]),  float(w[1]), float(w[2]))
                    loa.append(Atom(w[3], coo))
 
        # get the iterations
        cycle = False
        for line in f:
            w = line.split()
            if not cycle:
                # search for the cycle keyword
                if w[0] == 'cycle=':
                    cycle=True
                    n_coo=0
                    n_F=0
                    self.images.append(Cluster([]))
            else:
                if n_coo < len(self.definition):
                    n_coo += 1
                    coo = (float(w[0]),  float(w[1]), float(w[2]))
                    self[-1].append(Atom(w[3], coo))
                elif n_F < len(self.definition):
                    F = (float(w[0]),  float(w[1]), float(w[2]))
                    self[-1][n_F].F = F
                    n_F += 1
                    if n_F == len(self.definition):
                        cycle=False

class ViewmolTrajectory(Trajectory):
    """Write a trajectory for viewmol (http://viewmol.sourceforge.net)

    You can attach the writing to the Calculator:

    from gpaw.utilities.viewmol import ViewmolTrajectory

    c = Calculator()
    H2 = Cluster([Atom('H',(0,0,-.9)),Atom('H',(0,0,.9))])
    H2.SetCalculator(c)

    vmt = ViewmolTrajectory(H2)
    c.attach(vmt.add,100000)
    """
    def __init__(self, atoms, filename='trajectory.vmol',mode='w'):
        self.Ha = 1.0 / Hartree
        self.Ang = Angstrom
        
        self.n = 0
        self.atoms = atoms
        if mpi.rank == MASTER:
            self.file = open(filename, mode)
        else:
            self.file = open('/dev/null', mode)
        print(' $coord', 1. / Angstrom, file=self.file)
        self.write_viewmol_positions(self.file, self.atoms)
        print(' $grad', file=self.file)
    
    def __del__(self):
        print(' $end', file=self.file)
        self.file.close()

    def add(self, atoms=None):
        """Write current atomic position to the output file"""
        if atoms is None:
            atoms = self.atoms
        self.n += 1
        print('cycle=', self.n, end=' ', file=self.file)
        print('SCF energy=', \
              atoms.get_potential_energy() * Hartree, end=' ', file=self.file)
        forces = atoms.get_forces() * (Hartree / Angstrom)
        max_force = 0
        for f in forces:
            max_force = max(max_force, sqrt(np.sum(f * f)))
        print('|max dE/dxyz|=', max_force, file=self.file)
        self.write_viewmol_positions(self.file, atoms)
        for atom, f in zip(self.atoms, forces):
            print('%10.4g %10.4g %10.4g' % (f[0],f[1],f[2]), file=self.file)
        self.file.flush()
       
    def read(self, filename='trajectory.vmol', position=0):
        """Read atom configurations of step position"""
        self.file = None
        f = open(filename)
        # find coordinates
        loa = Cluster([])
        coords = False
        for l in f.readlines():
            if coords:
                w = l.split()
                loa.append(Atom(w[3].replace ("\n", "" ),
                                (float(w[0]), float(w[1]), float(w[2]))))
                
def write_viewmol(pt, filename, mode='w'):
    """Write PickleTrajectory as viewmol file."""
    atoms = pt[0]
    vt = ViewmolTrajectory(atoms, filename, mode)
    for atoms in pt:
        vt.add(atoms)
    del(vt)
    
