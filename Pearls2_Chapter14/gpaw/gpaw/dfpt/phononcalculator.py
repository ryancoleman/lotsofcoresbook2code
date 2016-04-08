"""This module provides an interface class for phonon calculations."""

__all__ = ["PhononCalculator"]

import time
import os
from math import pi, sqrt, sin

import numpy as np
import numpy.linalg as la

import ase.units as units
from ase.io.trajectory import PickleTrajectory

from gpaw import GPAW
from gpaw.mpi import serial_comm, rank
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.dfpt.poisson import PoissonSolver, FFTPoissonSolver
from gpaw.dfpt.responsecalculator import ResponseCalculator
from gpaw.dfpt.phononperturbation import PhononPerturbation
from gpaw.dfpt.wavefunctions import WaveFunctions
from gpaw.dfpt.dynamicalmatrix import DynamicalMatrix
from gpaw.dfpt.electronphononcoupling import ElectronPhononCoupling

from gpaw.symmetry import Symmetry

class PhononCalculator:
    """This class defines the interface for phonon calculations."""

    def __init__(self, calc, gamma=True, symmetry=False, e_ph=False,
                 communicator=serial_comm):
        """Inititialize class with a list of atoms.

        The atoms object must contain a converged ground-state calculation.

        The set of q-vectors in which the dynamical matrix will be calculated
        is determined from the ``symmetry`` kwarg. For now, only time-reversal
        symmetry is used to generate the irrecducible BZ.

        Add a little note on parallelization strategy here.

        Parameters
        ----------
        calc: str or Calculator
            Calculator containing a ground-state calculation.
        gamma: bool
            Gamma-point calculation with respect to the q-vector of the
            dynamical matrix. When ``False``, the Monkhorst-Pack grid from the
            ground-state calculation is used.
        symmetry: bool
            Use symmetries to reduce the q-vectors of the dynamcial matrix
            (None, False or True). The different options are equivalent to the
            old style options in a ground-state calculation (see usesymm).
        e_ph: bool
            Save the derivative of the effective potential.
        communicator: Communicator
            Communicator for parallelization over k-points and real-space
            domain.
        """

        # XXX
        assert symmetry in [None, False], "Spatial symmetries not allowed yet"

        if isinstance(calc, str):
            self.calc = GPAW(calc, communicator=serial_comm, txt=None)
        else:
            self.calc = calc

        cell_cv = self.calc.atoms.get_cell()
        setups = self.calc.wfs.setups
        # XXX - no clue how to get magmom - ignore it for the moment
        # m_av = magmom_av.round(decimals=3)  # round off
        # id_a = zip(setups.id_a, *m_av.T)
        id_a = setups.id_a

        if symmetry is None:
            self.symmetry = Symmetry(id_a, cell_cv, point_group=False,
                                     time_reversal=False)
        else:
            self.symmetry = Symmetry(id_a, cell_cv, point_group=False,
                                     time_reversal=True)

        # Make sure localized functions are initialized
        self.calc.set_positions()
        # Note that this under some circumstances (e.g. when called twice)
        # allocates a new array for the P_ani coefficients !!

        # Store useful objects
        self.atoms = self.calc.get_atoms()
        # Get rid of ``calc`` attribute
        self.atoms.calc = None

        # Boundary conditions
        pbc_c = self.calc.atoms.get_pbc()

        if np.all(pbc_c == False):
            self.gamma = True
            self.dtype = float
            kpts = None
            # Multigrid Poisson solver
            poisson_solver = PoissonSolver()
        else:
            if gamma:
                self.gamma = True
                self.dtype = float
                kpts = None
            else:
                self.gamma = False
                self.dtype = complex
                # Get k-points from ground-state calculation
                kpts = self.calc.input_parameters.kpts

            # FFT Poisson solver
            poisson_solver = FFTPoissonSolver(dtype=self.dtype)

        # K-point descriptor for the q-vectors of the dynamical matrix
        # Note, no explicit parallelization here.
        self.kd = KPointDescriptor(kpts, 1)
        self.kd.set_symmetry(self.atoms, self.symmetry)
        self.kd.set_communicator(serial_comm)

        # Number of occupied bands
        nvalence = self.calc.wfs.nvalence
        nbands = nvalence / 2 + nvalence % 2
        assert nbands <= self.calc.wfs.bd.nbands

        # Extract other useful objects
        # Ground-state k-point descriptor - used for the k-points in the
        # ResponseCalculator
        # XXX replace communicators when ready to parallelize
        kd_gs = self.calc.wfs.kd
        gd = self.calc.density.gd
        kpt_u = self.calc.wfs.kpt_u
        setups = self.calc.wfs.setups
        dtype_gs = self.calc.wfs.dtype

        # WaveFunctions
        wfs = WaveFunctions(nbands, kpt_u, setups, kd_gs, gd, dtype=dtype_gs)

        # Linear response calculator
        self.response_calc = ResponseCalculator(self.calc, wfs,
                                                dtype=self.dtype)

        # Phonon perturbation
        self.perturbation = PhononPerturbation(self.calc, self.kd,
                                               poisson_solver,
                                               dtype=self.dtype)

        # Dynamical matrix
        self.dyn = DynamicalMatrix(self.atoms, self.kd, dtype=self.dtype)

        # Electron-phonon couplings
        if e_ph:
            self.e_ph = ElectronPhononCoupling(self.atoms, gd, self.kd,
                                               dtype=self.dtype)
        else:
            self.e_ph = None

        # Initialization flag
        self.initialized = False

        # Parallel communicator for parallelization over kpts and domain
        self.comm = communicator

    def initialize(self):
        """Initialize response calculator and perturbation."""

        # Get scaled atomic positions
        spos_ac = self.atoms.get_scaled_positions()

        self.perturbation.initialize(spos_ac)
        self.response_calc.initialize(spos_ac)

        self.initialized = True

    def __getstate__(self):
        """Method used when pickling.

        Bound method attributes cannot be pickled and must therefore be deleted
        before an instance is dumped to file.

        """

        # Get state of object and take care of troublesome attributes
        state = dict(self.__dict__)
        state['kd'].__dict__['comm'] = serial_comm
        state.pop('calc')
        state.pop('perturbation')
        state.pop('response_calc')

        return state

    def run(self, qpts_q=None, clean=False, name=None, path=None):
        """Run calculation for atomic displacements and update matrix.

        Parameters
        ----------
        qpts: List
            List of q-points indices for which the dynamical matrix will be
            calculated (only temporary).

        """

        if not self.initialized:
            self.initialize()

        if self.gamma:
            qpts_q = [0]
        elif qpts_q is None:
            qpts_q = range(self.kd.nibzkpts)
        else:
            assert isinstance(qpts_q, list)

        # Update name and path attributes
        self.set_name_and_path(name=name, path=path)
        # Get string template for filenames
        filename_str = self.get_filename_string()

        # Delay the ranks belonging to the same k-point/domain decomposition
        # equally
        time.sleep(rank // self.comm.size)

        # XXX Make a single ground_state_contributions member function
        # Ground-state contributions to the force constants
        self.dyn.density_ground_state(self.calc)
        # self.dyn.wfs_ground_state(self.calc, self.response_calc)

        # Calculate linear response wrt q-vectors and displacements of atoms
        for q in qpts_q:

            if not self.gamma:
                self.perturbation.set_q(q)

            # First-order contributions to the force constants
            for a in self.dyn.indices:
                for v in [0, 1, 2]:

                    # Check if the calculation has already been done
                    filename = filename_str % (q, a, v)
                    # Wait for all sub-ranks to enter
                    self.comm.barrier()

                    if os.path.isfile(os.path.join(self.path, filename)):
                        continue

                    if self.comm.rank == 0:
                        fd = open(os.path.join(self.path, filename), 'w')

                    # Wait for all sub-ranks here
                    self.comm.barrier()

                    components = ['x', 'y', 'z']
                    symbols = self.atoms.get_chemical_symbols()
                    print("q-vector index: %i" % q)
                    print("Atom index: %i" % a)
                    print("Atomic symbol: %s" % symbols[a])
                    print("Component: %s" % components[v])

                    # Set atom and cartesian component of perturbation
                    self.perturbation.set_av(a, v)
                    # Calculate linear response
                    self.response_calc(self.perturbation)

                    # Calculate row of the matrix of force constants
                    self.dyn.calculate_row(self.perturbation,
                                           self.response_calc)

                    # Write force constants to file
                    if self.comm.rank == 0:
                        self.dyn.write(fd, q, a, v)
                        fd.close()

                    # Store effective potential derivative
                    if self.e_ph is not None:
                        v1_eff_G = self.perturbation.v1_G + \
                            self.response_calc.vHXC1_G
                        self.e_ph.v1_eff_qavG.append(v1_eff_G)

                    # Wait for the file-writing rank here
                    self.comm.barrier()

        # XXX
        # Check that all files are valid and collect in a single file
        # Remove the files
        if clean:
            self.clean()

    def get_atoms(self):
        """Return atoms."""

        return self.atoms

    def get_dynamical_matrix(self):
        """Return reference to ``dyn`` attribute."""

        return self.dyn

    def get_filename_string(self):
        """Return string template for force constant filenames."""

        name_str = (self.name + '.' + 'q_%%0%ii_' % len(str(self.kd.nibzkpts)) +
                    'a_%%0%ii_' % len(str(len(self.atoms))) + 'v_%i' + '.pckl')

        return name_str

    def set_atoms(self, atoms):
        """Set atoms to be included in the calculation.

        Parameters
        ----------
        atoms: list
            Can be either a list of strings, ints or ...
        """

        assert isinstance(atoms, list)

        if isinstance(atoms[0], str):
            assert np.all([isinstance(atom, str) for atom in atoms])
            sym_a = self.atoms.get_chemical_symbols()
            # List for atomic indices
            indices = []
            for type in atoms:
                indices.extend([a for a, atom in enumerate(sym_a)
                                if atom == type])
        else:
            assert np.all([isinstance(atom, int) for atom in atoms])
            indices = atoms

        self.dyn.set_indices(indices)

    def set_name_and_path(self, name=None, path=None):
        """Set name and path of the force constant files.

        name: str
            Base name for the files which the elements of the matrix of force
            constants will be written to.
        path: str
            Path specifying the directory where the files will be dumped.
        """

        if name is None:
            self.name = 'phonon.' + self.atoms.get_chemical_formula()
        else:
            self.name = name
        # self.name += '.nibzkpts_%i' % self.kd.nibzkpts

        if path is None:
            self.path = '.'
        else:
            self.path = path

        # Set corresponding attributes in the ``dyn`` attribute
        filename_str = self.get_filename_string()
        self.dyn.set_name_and_path(filename_str, self.path)

    def clean(self):
        """Delete generated files."""

        filename_str = self.get_filename_string()

        for q in range(self.kd.nibzkpts):
            for a in range(len(self.atoms)):
                for v in [0, 1, 2]:
                    filename = filename_str % (q, a, v)
                    if os.path.isfile(os.path.join(self.path, filename)):
                        os.remove(filename)

    def band_structure(self, path_kc, modes=False, acoustic=True):
        """Calculate phonon dispersion along a path in the Brillouin zone.

        The dynamical matrix at arbitrary q-vectors is obtained by Fourier
        transforming the real-space matrix. In case of negative eigenvalues
        (squared frequency), the corresponding negative frequency is returned.

        Parameters
        ----------
        path_kc: ndarray
            List of k-point coordinates (in units of the reciprocal lattice
            vectors) specifying the path in the Brillouin zone for which the
            dynamical matrix will be calculated.
        modes: bool
            Returns both frequencies and modes (mass scaled) when True.
        acoustic: bool
            Restore the acoustic sum-rule in the calculated force constants.
        """

        for k_c in path_kc:
            assert np.all(np.asarray(k_c) <= 1.0), \
                "Scaled coordinates must be given"

        # Assemble the dynanical matrix from calculated force constants
        self.dyn.assemble(acoustic=acoustic)
        # Get the dynamical matrix in real-space
        DR_lmn, R_clmn = self.dyn.real_space()

        # Reshape for the evaluation of the fourier sums
        shape = DR_lmn.shape
        DR_m = DR_lmn.reshape((-1,) + shape[-2:])
        R_cm = R_clmn.reshape((3, -1))

        # Lists for frequencies and modes along path
        omega_kn = []
        u_kn = []
        # Number of atoms included
        N = len(self.dyn.get_indices())

        # Mass prefactor for the normal modes
        m_inv_av = self.dyn.get_mass_array()

        for q_c in path_kc:

            # Evaluate fourier transform
            phase_m = np.exp(-2.j * pi * np.dot(q_c, R_cm))
            # Dynamical matrix in unit of Ha / Bohr**2 / amu
            D_q = np.sum(phase_m[:, np.newaxis, np.newaxis] * DR_m, axis=0)

            if modes:
                omega2_n, u_avn = la.eigh(D_q, UPLO='L')
                # Sort eigenmodes according to eigenvalues (see below) and
                # multiply with mass prefactor
                u_nav = u_avn[:, omega2_n.argsort()].T.copy() * m_inv_av
                # Multiply with mass prefactor
                u_kn.append(u_nav.reshape((3*N, -1, 3)))
            else:
                omega2_n = la.eigvalsh(D_q, UPLO='L')

            # Sort eigenvalues in increasing order
            omega2_n.sort()
            # Use dtype=complex to handle negative eigenvalues
            omega_n = np.sqrt(omega2_n.astype(complex))

            # Take care of imaginary frequencies
            if not np.all(omega2_n >= 0.):
                indices = np.where(omega2_n < 0)[0]
                print(("WARNING, %i imaginary frequencies at "
                       "q = (% 5.2f, % 5.2f, % 5.2f) ; (omega_q =% 5.3e*i)"
                       % (len(indices), q_c[0], q_c[1], q_c[2],
                          omega_n[indices][0].imag)))

                omega_n[indices] = -1 * np.sqrt(np.abs(omega2_n[indices].real))

            omega_kn.append(omega_n.real)

        # Conversion factor from sqrt(Ha / Bohr**2 / amu) -> eV
        s = units.Hartree**0.5 * units._hbar * 1.e10 / \
            (units._e * units._amu)**(0.5) / units.Bohr
        # Convert to eV and Ang
        omega_kn = s * np.asarray(omega_kn)
        if modes:
            u_kn = np.asarray(u_kn) * units.Bohr
            return omega_kn, u_kn

        return omega_kn

    def write_modes(self, q_c, branches=0, kT=units.kB*300, repeat=(1, 1, 1),
                    nimages=30, acoustic=True):
        """Write mode to trajectory file.

        The classical equipartioning theorem states that each normal mode has
        an average energy::

            <E> = 1/2 * k_B * T = 1/2 * omega^2 * Q^2

                =>

              Q = sqrt(k_B*T) / omega

        at temperature T. Here, Q denotes the normal coordinate of the mode.

        Parameters
        ----------
        q_c: ndarray
            q-vector of the modes.
        branches: int or list
            Branch index of calculated modes.
        kT: float
            Temperature in units of eV. Determines the amplitude of the atomic
            displacements in the modes.
        repeat: tuple
            Repeat atoms (l, m, n) times in the directions of the lattice
            vectors. Displacements of atoms in repeated cells carry a Bloch
            phase factor given by the q-vector and the cell lattice vector R_m.
        nimages: int
            Number of images in an oscillation.

        """

        if isinstance(branches, int):
            branch_n = [branches]
        else:
            branch_n = list(branches)

        # Calculate modes
        omega_n, u_n = self.band_structure([q_c], modes=True,
                                           acoustic=acoustic)

        # Repeat atoms
        atoms = self.atoms * repeat
        pos_mav = atoms.positions.copy()
        # Total number of unit cells
        M = np.prod(repeat)

        # Corresponding lattice vectors R_m
        R_cm = np.indices(repeat[::-1]).reshape(3, -1)[::-1]
        # Bloch phase
        phase_m = np.exp(2.j * pi * np.dot(q_c, R_cm))
        phase_ma = phase_m.repeat(len(self.atoms))

        for n in branch_n:
            omega = omega_n[0, n]
            u_av = u_n[0, n]  # .reshape((-1, 3))
            # Mean displacement at high T ?
            u_av *= sqrt(kT / abs(omega))

            mode_av = np.zeros((len(self.atoms), 3), dtype=self.dtype)
            indices = self.dyn.get_indices()
            mode_av[indices] = u_av
            mode_mav = (np.vstack([mode_av]*M) * phase_ma[:, np.newaxis]).real

            traj = PickleTrajectory('%s.mode.%d.traj' % (self.name, n), 'w')

            for x in np.linspace(0, 2*pi, nimages, endpoint=False):
                # XXX Is it correct to take out the sine component here ?
                atoms.set_positions(pos_mav + sin(x) * mode_mav)
                traj.write(atoms)

            traj.close()
