import os.path

import numpy as np

from ase.atom import Atom
from ase.atoms import Atoms
from ase.units import Bohr

import _gpaw
from gpaw import debug
from gpaw.external_potential import ElectrostaticPotential

class PointCharges(Atoms, ElectrostaticPotential):
    def __init__(self, *args, **kwargs):
        self.pc_nc = None
        self.charge_n = None
        Atoms.__init__(self, *args, **kwargs)
        
    def charge(self):
        """Return the summed charge of all point charges."""
        charge = 0
        for pc in self:
            charge += pc.charge
        return charge

    def read(self, file, filetype=None):
        """Read point charges from a file."""

        if hasattr(self, 'potential'):
            del(self.potential)
            del(self.gd)

        if filetype is None and isinstance(file, str):
            # estimate file type from name ending
            filetype = os.path.split(file)[-1].split('.')[-1]
        filetype = filetype.lower()

        if filetype == 'pc_info':
            self.read_PC_info(file)
        elif filetype == 'xyz':
            self.read_xyz(file)
        else:
            raise NotImplementedError('unknown file type "'+filetype+'"')
##        print "<PointCharges::read> found %d PC's" % len(self)


    def read_PC_info(self, file):
        """Read point charges from a PC_info file."""

        if isinstance(file, str):
            f = open(file)
        else:
            f = file

        lines = f.readlines()
        L0 = lines[0].split()
        del lines[0]

        for line in lines:
            words = line.split()
            if len(words) > 3:
                q, x, y, z = words[:4]
                self.append(PointCharge(position=(float(x) * Bohr,
                                                  float(y) * Bohr,
                                                  float(z) * Bohr),
                                        charge=float(q) ) )
            else:
                break

    def read_xyz(self, file):
        """Read point charges from a xyz file."""

        if isinstance(file, str):
            f = open(file)
        else:
            f = file

        lines = f.readlines()
        L0 = lines[0].split()
        del lines[0:2]

        n = int(L0[0])
        for i in range(n):
            words = lines[i].split()
            dummy, x, y, z, q = words[:5]
            self.append(PointCharge(position=(float(x),
                                              float(y),
                                              float(z)),
                                    charge=float(q) ) )

    def get_potential(self, gd=None):
        """Create the Coulomb potential on the grid."""

        if hasattr(self, 'potential'):
            if gd == self.gd or gd is None:
                # nothing changed
                return self.potential

        if gd is None:
            gd = self.gd

        assert gd.orthogonal
        potential = gd.empty()

        n = len(self)
        pc_nc = np.empty((n, 3))
        charge_n = np.empty((n))
        for a, pc in enumerate(self):
            pc_nc[a] = pc.position / Bohr
            charge_n[a] = pc.charge
        self.pc_nc = pc_nc
        self.charge_n = charge_n

        _gpaw.pc_potential(potential, pc_nc, charge_n,
                           gd.beg_c, gd.end_c, gd.h_cv.diagonal())

        # save grid descriptor and potential for future use
        self.potential = potential
        self.gd = gd

        return potential

    def get_nuclear_energy(self, nucleus):
        return -1. * nucleus.setup.Z * self.get_value(spos_c = nucleus.spos_c)

    def get_value(self, position=None, spos_c=None):
        """The potential value (as seen by an electron)
        at a certain grid point.

        position [Angstrom]
        spos_c scaled position on the grid"""
        if position is None:
            vr = np.dot(spos_c, self.gd.h_cv * self.gd.N_c)
        else:
            vr = position / Bohr

        if debug:
            v = 0
            for pc in self:
                # use c function XXXXX
                d = np.sqrt(np.sum((vr - pc.position / Bohr)**2))
                v -= pc.charge / d
        else:
            if self.pc_nc is None or self.charge_n is None:
                n = len(self)
                pc_nc = np.empty((n, 3))
                charge_n = np.empty((n))
                for a, pc in enumerate(self):
                    pc_nc[a] = pc.position / Bohr 
                    charge_n[a] = pc.charge
                self.pc_nc = pc_nc
                self.charge_n = charge_n

            v = _gpaw.pc_potential_value(vr, self.pc_nc, self.charge_n)
        return v

    def get_taylor(self, position=None, spos_c=None):
        """Get the Taylor expansion around a point

        position [Angstrom]
        output [Hartree, Hartree/Bohr]
        """
        if position is None:
            gd = self.gd
            pos = np.dot(spos_c, gd.h_cv * gd.N_c) * Bohr
        else:
            pos = position
        vr = pos / Bohr

        nabla = np.zeros((3))
        for pc in self:
            dist = vr - pc.position / Bohr
            d2 = np.sum(dist**2)
            nabla += dist * ( pc.charge / (d2 * np.sqrt(d2)) )
            
        # see spherical_harmonics.py for the assignment
        return [[self.get_value(position=pos)],
                np.array([nabla[1], nabla[2], nabla[0]])]
        
    def write(self, file='PC.xyz', filetype=None):

        if filetype is None and isinstance(file, str):
            # estimate file type from name ending
            filetype = os.path.split(file)[-1].split('.')[-1]
        filetype = filetype.lower()

        if filetype == 'xyz':
            self.write_xyz(file)
        else:
            raise NotImplementedError('unknown file type "'+filetype+'"')
##        print "<PointCharges::read> found %d PC's" % len(self)

    def write_xyz(self, file='PC.xyz'):
        if isinstance(file, str):
            f = open(file, 'w')
        else:
            f = file
        f.write('%d\nPoint charges\n' % len(self))
        for pc in self:
            (x, y, z) = pc.position
            q = pc.charge
            f.write('PC  %12.6g %12.6g %12.6g  %12.6g\n' % (x, y, z, q))
        f.close()

    def __eq__(self, other):
        """
        ASE atoms object does not compare charges. Hence, when calling
        GPAW.set(external=...) two identical PointCharge object with different
        charges won't trigger a reinitialization of the Hamiltionian object.
        """
        try:
            return Atoms.__eq__(self, other) and \
                   np.all(self.get_charges() == other.get_charges())
        except:
            return NotImplemented
                    
       

class PointCharge(Atom):
    def __init__(self, position, charge):
        Atom.__init__(self, position=position, charge=charge)
