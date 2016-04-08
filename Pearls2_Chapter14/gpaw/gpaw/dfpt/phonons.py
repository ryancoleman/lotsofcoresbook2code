"""Post processing module for dfpt force constants."""

import pickle

import numpy as np
import numpy.fft as fft

import ase.units as units
import ase.phonons as phonons

from gpaw.kpt_descriptor import KPointDescriptor

class Phonons(phonons.Phonons):
    """DFPT version of the ``Phonons`` class from ase.

    This class recycle the functionality from the ASE class.

    """

    def __init__(self, atoms, kpts, symmetry):
        """Initialize base class and attributes.

        Parameters
        ----------
        atoms: Atoms
            ASE atoms.
        kpts: tuple or list of tuples
            Shape of Monkhorst-Pack grid or list of k-points used in the dfpt
            calculation.
        symmetry: bool or None
            Symmetry parameter used in dfpt calculation.
            
        """

        # Init base class with ``Atoms`` object
        phonons.Phonons.__init__(atoms)

        # Create k-point descriptor
        self.kd = KPointDescriptor(kpts, 1)
        self.kd.set_symmetry(self.atoms, symmetry)

        # Overwrite ``N_c`` attribute
        self.N_c = tuple(self.kd.N_c)

        # Index of the gamma point -- for the acoustic sum-rule
        self.gamma_index = None
        if self.kd.gamma:
            self.gamma_index = 0
            self.dtype == float
        else:
            self.dtype == comples
            for k, k_c in enumerate(self.kd.ibzk_kc):
                if np.all(k_c == 0.):
                    self.gamma_index = k

        assert self.gamma_index is not None
        
    def run(self):
        """Overwrite base class member function."""

        raise RuntimeError, "Use only this class for post-processing"
    
    def read(self):
        """Read force constants from files."""

        # Data structure for force constants
        self.C_qaavv = [dict([(a, dict([(a_, np.zeros((3, 3), dtype=dtype))
                                        for a_ in self.indices]))
                              for a in self.indices])
                        for q in range(self.kd.nibzkpts)]
        
        assert self.name is not None
        assert self.path is not None
        
        for q in range(self.kd.nibzkpts):
            for a in self.indices:
                for v in [0, 1, 2]:
                    filename = self.name % (q, a, v)
                    try:
                        fd = open(os.path.join(self.path, filename))
                    except EOFError:
                        print(("Redo file %s "
                               % os.path.join(self.path, filename)))
                    C_qav_a = pickle.load(fd)
                    fd.close()
                    for a_ in self.indices:
                        self.C_qaavv[q][a][a_][v] = C_qav_a[a_]

    def assemble(self, acoustic=True):
        """Assemble dynamical matrix from the force constants attribute.

        The elements of the dynamical matrix are given by::

            D_ij(q) = 1/(M_i + M_j) * C_ij(q) ,
                      
        where i and j are collective atomic and cartesian indices.

        During the assembly, various symmetries of the dynamical matrix are
        enforced::

            1) Hermiticity
            2) Acoustic sum-rule
            3) D(q) = D*(-q)

        Parameters
        ----------
        acoustic: bool
            When True, the diagonal of the matrix of force constants is
            corrected to ensure that the acoustic sum-rule is fulfilled.
            
        """

        # Read force constants from files
        self.read()

        # Number of atoms included
        N = len(self.indices)
        
        # Assemble matrix of force constants
        self.C_q = []
        for q, C_aavv in enumerate(self.C_qaavv):

            C_avav = np.zeros((3*N, 3*N), dtype=self.dtype)
    
            for i, a in enumerate(self.indices):
                for j, a_ in enumerate(self.indices):
                    C_avav[3*i : 3*i + 3, 3*j : 3*j + 3] += C_aavv[a][a_]

            self.C_q.append(C_avav)

        # XXX Figure out in which order the corrections should be done
        # Make C(q) Hermitian
        for C in self.C_q:
            C *= 0.5
            C += C.T.conj()

        # Get matrix of force constants in the Gamma-point (is real!)
        C_gamma = self.C_q[self.gamma_index].real
        # Make Gamma-component real
        self.C_q[self.gamma_index] = C_gamma.copy()
            
        # Apply acoustic sum-rule if requested
        if acoustic:
            # Correct atomic diagonal for each q-vector
            for C in self.C_q:
                for a in range(N):
                    for a_ in range(N):
                        C[3*a : 3*a + 3, 3*a : 3*a + 3] -= \
                              C_gamma[3*a: 3*a+3, 3*a_: 3*a_+3]

            # Check sum-rule for Gamma-component in debug mode
            if debug:
                C = self.C_q[self.gamma_index]
                assert np.all(np.sum(C.reshape((3*N, N, 3)), axis=1) < 1e-15)

        
        # Move this bit to an ``unfold`` member function
        # XXX Time-reversal symmetry
        C_q = np.asarray(self.C_q)
        if self.kd.nibzkpts != self.kd.nbzkpts:
            self.D_k = np.concatenate((C_q[:0:-1].conjugate(), C_q))
        else:
            self.D_k = 0.5 * C_q
            self.D_k += self.D_k[::-1].conjugate()
            
        # Mass prefactor for the dynamical matrix
        m = self.atoms.get_masses()
        self.m_inv_av = np.repeat(m[self.indices]**-0.5, 3)
        M_avav = self.m_inv_av[:, np.newaxis] * self.m_inv_av

        for C in self.D_k:
            C *= M_avav

        self.assembled = True
       
    def real_space(self):
        """Fourier transform the dynamical matrix to real-space."""

        if not self.assembled:
            self.assemble()

        # Shape of q-point grid
        N_c = self.N_c

        # Reshape before Fourier transforming
        shape = self.D_k.shape
        Dq_lmn = self.D_k.reshape(N_c + shape[1:])
        DR_lmn = fft.ifftn(fft.ifftshift(Dq_lmn, axes=(0, 1, 2)), axes=(0, 1, 2))

        if debug:
            # Check that D_R is real enough
            assert np.all(DR_lmn.imag < 1e-8)
            
        DR_lmn = DR_lmn.real

        # Corresponding R_m vectors in units of the basis vectors
        R_cm = np.indices(N_c).reshape(3, -1)
        N1_c = np.array(N_c)[:, np.newaxis]        
        R_cm += N1_c // 2
        R_cm %= N1_c
        R_cm -= N1_c // 2
        R_clmn = R_cm.reshape((3,) + N_c)

        return DR_lmn, R_clmn
    
    def fourier_interpolate(self, N_c):
        """Fourier interpolate dynamical matrix onto a finer q-vector grid."""

        # Move this member function to the ASE class
        raise NotImplementedError    
        





