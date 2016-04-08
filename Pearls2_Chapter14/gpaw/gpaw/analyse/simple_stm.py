from __future__ import print_function
from math import sqrt

import numpy as np

from ase.atoms import Atoms
from ase.units import Bohr, Hartree
from ase.io.cube import write_cube
from ase.io.plt import write_plt, read_plt

import gpaw.mpi as mpi
from gpaw.mpi import MASTER
from gpaw.grid_descriptor import GridDescriptor


class SimpleStm:

    """Simple STM object to simulate STM pictures.

    The simulation uses either a single pseudo-wavefunction (PWF)
    or the PWFs inside the given bias range."""

    def __init__(self, atoms):

        self.file = None
        self.is_wf = False
        self.bias = None
        self.ldos = None
        self.heights = None

        if isinstance(atoms, str):
            self.read_3D(atoms)
            self.calc = None
        else:
            if isinstance(atoms, Atoms):
                self.calc = atoms.get_calculator()
            else:
                self.calc = atoms
            self.calc.converge_wave_functions()

            self.gd = self.calc.wfs.gd
            self.offset_c = [int(not a) for a in self.gd.pbc_c]

    def calculate_ldos(self, bias):
        """bias is the n, k, s list/tuple."""
        if self.calc is None:
            return

        self.bias = bias
        self.is_wf = True
        self.ldos = self.gd.zeros()

        if hasattr(bias, '__len__') and len(bias) == 3:
            n, k, s = bias
            # only a single wf requested
            u = self.calc.get_myu(k, s)
            if u is not None:
                self.add_wf_to_ldos(n, u, weight=1)

        else:
            # energy bias
            try:
                if self.calc.occupations.fixmagmom is True:
                    efermi_s = self.calc.get_fermi_levels()
                else:
                    efermi_s = np.array([self.calc.get_fermi_level()] * 2)
            except:
                efermi_s = np.array([self.calc.get_homo_lumo().mean()] * 2)

            if isinstance(bias, (int, long, float)):
                # bias given
                if bias > 0:
                    # positive bias = negative tip
                    # -> probe unoccupied states
                    emin_s = efermi_s
                    emax_s = efermi_s + bias
                    occupied = False
                else:
                    # negative bias = positive tip
                    # -> probe occupied states
                    emin_s = efermi_s + bias
                    emax_s = efermi_s
                    occupied = True
            else:
                # emin and emax given
                emin, emax = bias
                if abs(emin) > abs(emax):
                    occupied = True
                else:
                    occupied = False
                emin_s = np.array([emin + efermi_s] * 2)
                emax_s = np.array([emax + efermi_s] * 2)

            emin_s /= Hartree
            emax_s /= Hartree

            for u in range(len(self.calc.wfs.kpt_u)):
                kpt = self.calc.wfs.kpt_u[u]
                emin = emin_s[kpt.s]
                emax = emax_s[kpt.s]
                for n, eps in enumerate(kpt.eps_n):
                    if eps > emin and eps < emax:
                        if occupied:
                            weight = kpt.f_n[n]
                        else:
                            weight = kpt.weight - kpt.f_n[n]
                        self.add_wf_to_ldos(n, u, weight)

    def add_wf_to_ldos(self, n, u, weight=None):
        """Add the wf with given kpoint and spin to the ldos"""
        kpt = self.calc.wfs.kpt_u[u]
        psi = kpt.psit_nG[n]
        w = weight
        if w is None:
            w = kpt.weight
# print "w=", w, kpt.weight
        self.ldos += w * (psi * np.conj(psi)).real

    def write_3D(self, bias, file, filetype=None):
        """Write the density as a 3D file.

        Units: [e/A^3]"""
        self.calculate_ldos(bias)
        self.calc.wfs.kd.comm.sum(self.ldos)
        ldos = self.gd.collect(self.ldos)
# print "write: integrated =", self.gd.integrate(self.ldos)

        if mpi.rank != MASTER:
            return

        if filetype is None:
            # estimate file type from name ending
            filetype = file.split('.')[-1]
        filetype.lower()

        if filetype == 'cube':
            write_cube(file, self.calc.get_atoms(), ldos / Bohr ** 3)
        elif filetype == 'plt':
            write_plt(file, self.calc.get_atoms(), ldos / Bohr ** 3)
        else:
            raise NotImplementedError('unknown file type "' + filetype + '"')

    def read_3D(self, file, filetype=None):
        """Read the density from a 3D file"""

        if filetype is None:
            # estimate file type from name ending
            filetype = file.split('.')[-1]
        filetype.lower()

        if filetype == 'plt':
            data, cell = read_plt(file)

            pbc_c = [True, True, True]
            N_c = np.array(data.shape)
            for c in range(3):
                if N_c[c] % 2 == 1:
                    pbc_c[c] = False
                    N_c[c] += 1
            self.gd = GridDescriptor(N_c, cell.diagonal() / Bohr, pbc_c)
            self.offset_c = [int(not a) for a in self.gd.pbc_c]

        else:
            raise NotImplementedError('unknown file type "' + filetype + '"')

        self.file = file
        self.ldos = np.array(data * Bohr ** 3, np.float)
# print "read: integrated =", self.gd.integrate(self.ldos)

    def current_to_density(self, current):
        """The connection between density n and current I

        n [e/Angstrom^3] = 0.0002 sqrt(I [nA])

        as given in Hofer et al., RevModPhys 75 (2003) 1287
        """
        return 0.0002 * sqrt(current)

    def density_to_current(self, density):
        return 5000. * density ** 2

    def scan_const_current(self, current, bias=None,
                           interpolate=False, hmax=None):
        """Get the height image for constant current I [nA].

        hmax is the maximal height to consider
        """
        return self.scan_const_density(self.current_to_density(current),
                                       bias, interpolate, hmax)

    def scan_const_density(self, density, bias, interpolate=False, hmax=None):
        """Get the height image for constant density [e/Angstrom^3].
        """

        self.calculate_ldos(bias)

        self.density = density

        gd = self.gd
        h_c = [np.linalg.norm(gd.h_cv[c]) for c in range(3)]
        nx, ny = (gd.N_c - self.offset_c)[:2]

        # each cpu will have the full array, but works on its
        # own part only
        heights = np.zeros((nx, ny)) - 1
        if hmax is None:
            hmax = h_c[2] * self.ldos.shape[2] + h_c[2] / 2.
        else:
            hmax /= Bohr
        ihmax = min(gd.end_c[2] - 1, int(hmax / h_c[2]))

        for i in range(gd.beg_c[0], gd.end_c[0]):
            ii = i - gd.beg_c[0]
            for j in range(gd.beg_c[1], gd.end_c[1]):
                jj = j - gd.beg_c[1]

                zline = self.ldos[ii, jj]

                # check from above until you find the required density
                for k in range(ihmax, gd.beg_c[2] - 1, -1):
                    kk = k - gd.beg_c[2]
                    if zline[kk] > density:
                        heights[i - self.offset_c[0],
                                j - self.offset_c[1]] = k
                        break

        # collect the results
        gd.comm.max(heights)

        if interpolate:
            # collect the full grid to enable interpolation
            fullgrid = gd.collect(self.ldos)

            kmax = self.ldos.shape[2] - 1
            for i in range(gd.beg_c[0], gd.end_c[0]):
                ii = i - gd.beg_c[0]
                i -= self.offset_c[0]
                for j in range(gd.beg_c[1], gd.end_c[1]):
                    jj = j - gd.beg_c[1]
                    j -= self.offset_c[1]
                    if heights[i, j] > 0:
                        if heights[i, j] < kmax:
                            c1 = fullgrid[i, j, int(heights[i, j])]
                            c2 = fullgrid[i, j, int(heights[i, j]) + 1]
                            k = heights[i, j] + (density - c1) / (c2 - c1)
                        else:
                            k = kmax

        self.heights = np.where(heights > 0,
                               (heights + self.offset_c[2]) * h_c[2], -1)

        return heights

    def write(self, file=None):
        """Write STM data to a file in gnuplot readable tyle."""

        if mpi.rank != MASTER:
            return

        xvals, yvals, heights = self.pylab_contour()
        nx, ny = heights.shape[:2]

        if file is None:
            n, k, s = self.bias
            fname = 'stm_n%dk%ds%d.dat' % (n, k, s)
        else:
            fname = file
        f = open(fname, 'w')

        try:
            import datetime
            print('#', datetime.datetime.now().ctime(), file=f)
        except:
            pass
        print('# Simulated STM picture', file=f)
        if hasattr(self, 'file'):
            print('# density read from', self.file, file=f)
        else:
            if self.is_wf:
                print('# pseudo-wf n=%d k=%d s=%d' % tuple(self.bias), file=f)
            else:
                print('# bias=', self.bias, '[eV]', file=f)
        print('#', file=f)
        print('# density=', self.density, '[e/Angstrom^3]', end=' ', file=f)
        print(
            '(current=', self.density_to_current(self.density), '[nA])', file=f)
        print('# x[Angs.]   y[Angs.]     h[Angs.] (-1 is not found)', file=f)
        for i in range(nx):
            for j in range(ny):
                if heights[i, j] == -1:
                    height = -1
                else:
                    height = heights[i, j] * Bohr
                print('%10g %10g %12g' % (yvals[j], xvals[i], height), file=f)
            print(file=f)
        f.close()

    def pylab_contour(self):
        """Return the countour to be plotted using pylab."""

        nx, ny = self.heights.shape[:2]
        h_c = np.array([np.linalg.norm(self.gd.h_cv[c])
                        for c in range(3)]) * Bohr
        # the lowest point is not stored for non-periodic BCs
        xvals = [(i + self.offset_c[0]) * h_c[0] for i in range(nx)]
        yvals = [(i + self.offset_c[1]) * h_c[1] for i in range(ny)]
        heights = self.heights * Bohr

        # pylab interprets heights[y_i][x_i]
        return np.array(xvals), np.array(yvals), heights.swapaxes(0, 1)
