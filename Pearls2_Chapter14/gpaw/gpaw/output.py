import os
import sys
import time
import platform
from math import log

import numpy as np
import ase
from ase.units import Bohr, Hartree
from ase.data import chemical_symbols
from ase.version import version as ase_version
from ase.utils import devnull
from ase.parallel import get_txt

import gpaw
import _gpaw
from gpaw.version import version
from gpaw.utilities.memory import maxrss
from gpaw import dry_run, extra_parameters


class PAWTextOutput:
    """Class for handling all text output."""

    def __init__(self):
        self.txt = devnull

    def set_txt(self, txt):
        """Set the stream for text output.

        If `txt` is not a stream-object, then it must be one of:

        * None:  Throw output away.
        * '-':  Use stdout (``sys.stdout``) on master, elsewhere throw away.
        * A filename:  Open a new file on master, elsewhere throw away.
        """

        self.txt = get_txt(txt, self.wfs.world.rank)
        self.print_header()

    def text(self, *args, **kwargs):
        self.txt.write(kwargs.get('sep', ' ').join([str(arg)
                                                    for arg in args]) +
                       kwargs.get('end', '\n'))

    def print_header(self):
        self.text()
        self.text('  ___ ___ ___ _ _ _  ')
        self.text(' |   |   |_  | | | | ')
        self.text(' | | | | | . | | | | ')
        self.text(' |__ |  _|___|_____| ', version)
        self.text(' |___|_|             ')
        self.text()

        uname = platform.uname()
        self.text('User: ', os.getenv('USER', '???') + '@' + uname[1])
        self.text('Date: ', time.asctime())
        self.text('Arch: ', uname[4])
        self.text('Pid:  ', os.getpid())
        self.text('gpaw: ', os.path.dirname(gpaw.__file__))
        
        # Find C-code:
        c = getattr(_gpaw, '__file__', None)
        if not c:
            c = sys.executable
        self.text('_gpaw:', os.path.normpath(c))
                  
        self.text('ase:   %s (version %s)' %
                  (os.path.dirname(ase.__file__), ase_version))
        self.text('numpy: %s (version %s)' %
                  (os.path.dirname(np.__file__), np.version.version))
        try:
            import scipy as sp
            self.text('scipy: %s (version %s)' %
                      (os.path.dirname(sp.__file__), sp.version.version))
            # Explicitly deleting SciPy seems to remove garbage collection
            # problem of unknown cause
            del sp
        except ImportError:
            self.text('scipy: Not available')
        self.text('units: Angstrom and eV')
        self.text('cores:', self.wfs.world.size)

        if gpaw.debug:
            self.text('DEBUG MODE')

        if extra_parameters:
            self.text('Extra parameters:', extra_parameters)

    def print_cell_and_parameters(self):
        self.plot_atoms(self.atoms)
        self.print_unit_cell(self.atoms.get_positions() / Bohr)
        self.print_parameters()

    def print_unit_cell(self, pos_ac):
        self.text()
        self.text('Unit Cell:')
        self.text('           Periodic     X           Y           Z' +
                  '      Points  Spacing')
        self.text('  -----------------------------------------------' +
                  '---------------------')
        gd = self.wfs.gd
        h_c = gd.get_grid_spacings()
        pbc_c = self.atoms.pbc
        for c in range(3):
            self.text('  %d. axis:    %s  %10.6f  %10.6f  %10.6f   %3d   %8.4f'
                      % ((c + 1, ['no ', 'yes'][int(pbc_c[c])]) +
                         tuple(Bohr * gd.cell_cv[c]) +
                         (gd.N_c[c], Bohr * h_c[c])))
        self.text()

    def print_positions(self):
        t = self.text
        t()
        t('Positions:')
        symbols = self.atoms.get_chemical_symbols()
        for a, pos_c in enumerate(self.atoms.get_positions()):
            symbol = symbols[a]
            t('%3d %-2s %9.4f %9.4f %9.4f' % ((a, symbol) + tuple(pos_c)))
        t()

    def print_parameters(self):
        t = self.text
        p = self.input_parameters

        # Write PAW setup information in order of appearance:
        ids = set()
        for id in self.wfs.setups.id_a:
            if id in ids:
                continue
            ids.add(id)
            setup = self.wfs.setups.setups[id]
            setup.print_info(self.text)
            basis_descr = setup.get_basis_description()
            t(basis_descr)
            t()

        t('Using the %s Exchange-Correlation Functional.'
          % self.hamiltonian.xc.name)
        if self.wfs.nspins == 2:
            t('Spin-Polarized Calculation.')
            t('Magnetic Moment:  (%.6f, %.6f, %.6f)' %
              tuple(self.density.magmom_av.sum(0)), end='')
            if self.occupations.fixmagmom:
                t('(fixed)')
            else:
                t()
        else:
            t('Spin-Paired Calculation')
        t('Total Charge: %.6f' % p['charge'])
        t('Fermi Temperature: %.6f' % (self.occupations.width * Hartree))
        self.wfs.summary(self.txt)
        t('Eigensolver: %s' % self.wfs.eigensolver)

        self.hamiltonian.summary(self.txt)

        t('Reference Energy: %.6f' % (self.wfs.setups.Eref * Hartree))
        t()

        nibzkpts = self.wfs.kd.nibzkpts

        # Print parallelization details
        t('Total number of cores used: %d' % self.wfs.world.size)
        if self.wfs.kd.comm.size > 1:  # kpt/spin parallization
            if self.wfs.nspins == 2 and nibzkpts == 1:
                t('Parallelization over spin')
            elif self.wfs.nspins == 2:
                t('Parallelization over k-points and spin: %d' %
                  self.wfs.kd.comm.size)
            else:
                t('Parallelization over k-points: %d' %
                  self.wfs.kd.comm.size)
        if self.wfs.gd.comm.size > 1:  # domain parallelization
            t('Domain Decomposition: %d x %d x %d' %
              tuple(self.wfs.gd.parsize_c))
        if self.wfs.bd.comm.size > 1:  # band parallelization
            t('Parallelization over states: %d'
              % self.wfs.bd.comm.size)

        if self.wfs.mode == 'fd':
            # XXX Why is this not in wfs summary? -askhl
            # Also, why wouldn't PW mode use the diagonalizer?
            if self.wfs.diagksl.buffer_size is not None:
                t('MatrixOperator buffer_size (KiB): %d'
                  % self.wfs.diagksl.buffer_size)
            else:
                t('MatrixOperator buffer_size: default value or \n' +
                  ' %s see value of nblock in input file' % (26 * ' '))
            diagonalizer_layout = self.wfs.diagksl.get_description()
            t('Diagonalizer layout: ' + diagonalizer_layout)
            orthonormalizer_layout = self.wfs.orthoksl.get_description()
            t('Orthonormalizer layout: ' + orthonormalizer_layout)
        t()

        self.wfs.kd.symmetry.print_symmetries(self.txt)

        t(self.wfs.kd.description)
        t(('%d k-point%s in the Irreducible Part of the Brillouin Zone') %
          (nibzkpts, ' s'[1:nibzkpts]))
        
        kd = self.wfs.kd
        w_k = kd.weight_k * kd.nbzkpts
        assert np.allclose(w_k, w_k.round())
        w_k = w_k.round()
        
        t()
        t('          k-points in crystal coordinates                weights')
        for k in range(nibzkpts):
            if k < 10 or k == nibzkpts - 1:
                t('%4d:   %12.8f  %12.8f  %12.8f     %6d/%d' %
                  ((k,) + tuple(kd.ibzk_kc[k]) + (w_k[k], kd.nbzkpts)))
            elif k == 10:
                t('          ...')
        t()

        if self.scf.fixdensity > self.scf.maxiter:
            t('Fixing the initial density')
        else:
            mixer = self.density.mixer
            t('Mixer Type:', mixer.__class__.__name__)
            t('Linear Mixing Parameter: %g' % mixer.beta)
            t('Mixing with %d Old Densities' % mixer.nmaxold)
            if mixer.weight == 1:
                t('No Damping of Long Wave Oscillations')
            else:
                t('Damping of Long Wave Oscillations: %g' % mixer.weight)

        cc = p['convergence']
        t()
        t('Convergence Criteria:')
        t('    Total Energy Change: %g eV / electron' %
          (cc['energy']))
        t('    Integral of Absolute Density Change: %g electrons' %
          cc['density'])
        t('    Integral of Absolute Eigenstate Change: %g eV^2' %
          cc['eigenstates'])
        t('Number of Atoms: %d' % len(self.wfs.setups))
        t('Number of Atomic Orbitals: %d' % self.wfs.setups.nao)
        if self.nbands_parallelization_adjustment != 0:
            t('Adjusting Number of Bands by %+d to Match Parallelization'
              % self.nbands_parallelization_adjustment)
        t('Number of Bands in Calculation: %d' % self.wfs.bd.nbands)
        t('Bands to Converge: ', end='')
        if cc['bands'] == 'occupied':
            t('Occupied States Only')
        elif cc['bands'] == 'all':
            t('All')
        else:
            t('%d Lowest Bands' % cc['bands'])
        t('Number of Valence Electrons: %g'
          % (self.wfs.setups.nvalence - p.charge))

    def print_converged(self, iter):
        t = self.text
        t('------------------------------------')
        t('Converged After %d Iterations.' % iter)

        t()
        self.print_all_information()

    def print_all_information(self):
        t = self.text
        if len(self.atoms) == 1:
            t('Energy Contributions Relative to Reference Atom:', end='')
        else:
            t('Energy Contributions Relative to Reference Atoms:', end='')
        t('(reference = %.6f)' % (self.wfs.setups.Eref * Hartree))

        t('-------------------------')

        energies = [('Kinetic:      ', self.hamiltonian.Ekin),
                    ('Potential:    ', self.hamiltonian.Epot),
                    ('External:     ', self.hamiltonian.Eext),
                    ('XC:           ', self.hamiltonian.Exc),
                    ('Entropy (-ST):', -self.hamiltonian.S),
                    ('Local:        ', self.hamiltonian.Ebar)]

        for name, e in energies:
            t('%-14s %+11.6f' % (name, Hartree * e))

        t('-------------------------')
        t('Free Energy:   %+11.6f' % (Hartree * self.hamiltonian.Etot))
        t('Zero Kelvin:   %+11.6f' % (Hartree * (self.hamiltonian.Etot +
                                                 0.5 * self.hamiltonian.S)))
        t()
        self.occupations.print_fermi_level(self.txt)

        self.print_eigenvalues()

        self.hamiltonian.xc.summary(self.txt)

        t()

        try:
            dipole = self.get_dipole_moment()
        except NotImplementedError:
            pass
        else:
            if self.density.charge == 0:
                t('Dipole Moment: %s' % dipole)
            else:
                t('Center of Charge: %s' % (dipole / abs(self.density.charge)))

        try:
            c = self.hamiltonian.poisson.corrector.c
            epsF = self.occupations.fermilevel
        except AttributeError:
            pass
        else:
            wf_a = -epsF * Hartree - dipole[c]
            wf_b = -epsF * Hartree + dipole[c]
            t('Dipole-corrected work function: %f, %f' % (wf_a, wf_b))

        if self.wfs.nspins == 2:
            t()
            magmom = self.occupations.magmom
            t('Total Magnetic Moment: %f' % magmom)
            try:
                # XXX This doesn't always work, HGH, SIC, ...
                sc = self.density.get_spin_contamination(self.atoms,
                                                         int(magmom < 0))
                t('Spin contamination: %f electrons' % sc)
            except (TypeError, AttributeError):
                pass
            t('Local Magnetic Moments:')
            for a, mom in enumerate(self.get_magnetic_moments()):
                t(a, mom)
            t()
        elif not self.wfs.collinear:
            self.txt.write('Local Magnetic Moments:\n')
            for a, mom_v in enumerate(self.get_magnetic_moments()):
                self.txt.write('%4d  (%.3f, %.3f, %.3f)\n' %
                               (a, mom_v[0], mom_v[1], mom_v[2]))

    def print_iteration(self, iter):
        # Output from each iteration:
        t = self.text

        nvalence = self.wfs.setups.nvalence - self.input_parameters.charge
        if nvalence > 0:
            eigerr = self.scf.eigenstates_error * Hartree**2 / nvalence
        else:
            eigerr = 0.0

        if self.verbose != 0:
            T = time.localtime()
            t()
            t('------------------------------------')
            t('iter: %d %d:%02d:%02d' % (iter, T[3], T[4], T[5]))
            t()
            t('Poisson Solver Converged in %d Iterations' %
              self.hamiltonian.npoisson)
            t('Fermi Level Found  in %d Iterations' % self.occupations.niter)
            t('Error in Wave Functions: %.13f' % eigerr)
            t()
            self.print_all_information()

        else:
            if iter == 1:
                header = """\
                     log10-error:    Total        Iterations:
           Time      WFS    Density  Energy       Fermi  Poisson"""
                if self.wfs.nspins == 2:
                    header += '  MagMom'
                if self.scf.max_force_error is not None:
                    l1 = header.find('Total')
                    header = header[:l1] + '       ' + header[l1:]
                    l2 = header.find('Energy')
                    header = header[:l2] + 'Force  ' + header[l2:]
                t(header)

            T = time.localtime()

            if eigerr == 0.0:
                eigerr = ''
            else:
                eigerr = '%+.2f' % (log(eigerr) / log(10))

            denserr = self.density.mixer.get_charge_sloshing()
            if denserr is None or denserr == 0 or nvalence == 0:
                denserr = ''
            else:
                denserr = '%+.2f' % (log(denserr / nvalence) / log(10))

            niterocc = self.occupations.niter
            if niterocc == -1:
                niterocc = ''
            else:
                niterocc = '%d' % niterocc

            if self.hamiltonian.npoisson == 0:
                niterpoisson = ''
            else:
                niterpoisson = str(self.hamiltonian.npoisson)

            t('iter: %3d  %02d:%02d:%02d %6s %6s  ' %
              (iter,
               T[3], T[4], T[5],
               eigerr,
               denserr), end='')

            if self.scf.max_force_error is not None:
                if self.scf.force_error is not None:
                    t('  %+.2f' %
                      (log(self.scf.force_error) / log(10)), end='')
                else:
                    t('       ', end='')

            t('%11.6f    %-5s  %-7s' %
              (Hartree * (self.hamiltonian.Etot + 0.5 * self.hamiltonian.S),
               niterocc,
               niterpoisson), end='')

            if self.wfs.nspins == 2:
                t('  %+.4f' % self.occupations.magmom, end='')

            t()

        self.txt.flush()

    def print_forces(self):
        if self.forces.F_av is None:
            return
        t = self.text
        t()
        t('Forces in eV/Ang:')
        c = Hartree / Bohr
        symbols = self.atoms.get_chemical_symbols()
        for a, symbol in enumerate(symbols):
            t('%3d %-2s %10.5f %10.5f %10.5f' %
              ((a, symbol) + tuple(self.forces.F_av[a] * c)))

    def print_eigenvalues(self):
        """Print eigenvalues and occupation numbers."""
        self.text(eigenvalue_string(self))

    def plot_atoms(self, atoms):
        self.text(plot(atoms))

    def __del__(self):
        """Destructor:  Write timing output before closing."""
        if not dry_run:
            mr = maxrss()
            if mr > 0:
                if mr < 1024.0**3:
                    self.text('Memory usage: %.2f MiB' % (mr / 1024.0**2))
                else:
                    self.text('Memory usage: %.2f GiB' % (mr / 1024.0**3))

            self.timer.write(self.txt)


def eigenvalue_string(paw, comment=None):
    """
    Write eigenvalues and occupation numbers into a string.
    The parameter comment can be used to comment out non-numers,
    for example to escape it for gnuplot.
    """

    if not comment:
        comment = ' '

    if len(paw.wfs.kd.ibzk_kc) > 1:
        # not implemented yet:
        return ''

    s = ''
    if paw.wfs.nspins == 1:
        s += comment + 'Band   Eigenvalues  Occupancy\n'
        eps_n = paw.get_eigenvalues(kpt=0, spin=0)
        f_n = paw.get_occupation_numbers(kpt=0, spin=0)
        if paw.wfs.world.rank == 0:
            for n in range(paw.wfs.bd.nbands):
                s += ('%4d   %10.5f  %10.5f\n' % (n, eps_n[n], f_n[n]))
    else:
        s += comment + '                 Up                     Down\n'
        s += comment + 'Band  Eigenvalues  Occupancy  Eigenvalues  Occupancy\n'
        epsa_n = paw.get_eigenvalues(kpt=0, spin=0, broadcast=False)
        epsb_n = paw.get_eigenvalues(kpt=0, spin=1, broadcast=False)
        fa_n = paw.get_occupation_numbers(kpt=0, spin=0, broadcast=False)
        fb_n = paw.get_occupation_numbers(kpt=0, spin=1, broadcast=False)
        if paw.wfs.world.rank == 0:
            for n in range(paw.wfs.bd.nbands):
                s += (' %4d  %11.5f  %9.5f  %11.5f  %9.5f\n' %
                      (n, epsa_n[n], fa_n[n], epsb_n[n], fb_n[n]))
    return s


def plot(atoms):
    """Ascii-art plot of the atoms."""

##   y
##   |
##   .-- x
##  /
## z

    cell_cv = atoms.get_cell()
    if (cell_cv - np.diag(cell_cv.diagonal())).any():
        atoms = atoms.copy()
        atoms.cell = [1, 1, 1]
        atoms.center(vacuum=2.0)
        cell_cv = atoms.get_cell()
        plot_box = False
    else:
        plot_box = True

    cell = np.diagonal(cell_cv) / Bohr
    positions = atoms.get_positions() / Bohr
    numbers = atoms.get_atomic_numbers()

    s = 1.3
    nx, ny, nz = n = (s * cell * (1.0, 0.25, 0.5) + 0.5).astype(int)
    sx, sy, sz = n / cell
    grid = Grid(nx + ny + 4, nz + ny + 1)
    positions = (positions % cell + cell) % cell
    ij = np.dot(positions, [(sx, 0), (sy, sy), (0, sz)])
    ij = np.around(ij).astype(int)
    for a, Z in enumerate(numbers):
        symbol = chemical_symbols[Z]
        i, j = ij[a]
        depth = positions[a, 1]
        for n, c in enumerate(symbol):
            grid.put(c, i + n + 1, j, depth)
    if plot_box:
        k = 0
        for i, j in [(1, 0), (1 + nx, 0)]:
            grid.put('*', i, j)
            grid.put('.', i + ny, j + ny)
            if k == 0:
                grid.put('*', i, j + nz)
            grid.put('.', i + ny, j + nz + ny)
            for y in range(1, ny):
                grid.put('/', i + y, j + y, y / sy)
                if k == 0:
                    grid.put('/', i + y, j + y + nz, y / sy)
            for z in range(1, nz):
                if k == 0:
                    grid.put('|', i, j + z)
                grid.put('|', i + ny, j + z + ny)
            k = 1
        for i, j in [(1, 0), (1, nz)]:
            for x in range(1, nx):
                if k == 1:
                    grid.put('-', i + x, j)
                grid.put('-', i + x + ny, j + ny)
            k = 0
    return '\n'.join([''.join([chr(x) for x in line])
                      for line in np.transpose(grid.grid)[::-1]])


class Grid:
    def __init__(self, i, j):
        self.grid = np.zeros((i, j), np.int8)
        self.grid[:] = ord(' ')
        self.depth = np.zeros((i, j))
        self.depth[:] = 1e10

    def put(self, c, i, j, depth=1e9):
        if depth < self.depth[i, j]:
            self.grid[i, j] = ord(c)
            self.depth[i, j] = depth
