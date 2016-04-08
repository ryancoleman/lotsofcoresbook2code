"""This module provides a class for assembling the dynamical matrix."""

__all__ = ["DynamicalMatrix"]

from math import sqrt, pi
import os
import pickle

import numpy as np
import numpy.fft as fft

from gpaw import debug
from gpaw.mpi import serial_comm
from gpaw.utilities import unpack, unpack2

class DynamicalMatrix:
    """Class for assembling the dynamical matrix from first-order responses.

    The second order derivative of the total energy with respect to atomic
    displacements (for periodic systems collective atomic displacemnts
    characterized by a q-vector) can be obtained from an expression involving
    the first-order derivatives of the density and the wave-functions.
    
    Each of the various contributions to the second order derivative of the
    total energy are implemented in separate functions.
    
    """
    
    def __init__(self, atoms, kd, dtype=float):
        """Inititialize class with a list of atoms.

        Parameters
        ----------
        atoms: Atoms
            List of atoms for the system.
        kd: KPointDescriptor
            Descriptor for the q-vector grid.

        """

        # Store useful objects
        self.atoms = atoms
        self.kd = kd
        self.dtype = dtype
        self.masses = atoms.get_masses()
        # Array with inverse sqrt of masses repeated to match shape of mode
        # arrays        
        self.m_inv_av = None
        
        # String template for files with force constants
        self.name = None
        # Path to directory with force constant files
        self.path = None
        
        # List of atomic indices to be included (default is all atoms)
        self.indices = range(len(self.atoms))
        
        # Index of the gamma point -- for the acoustic sum-rule
        self.gamma_index = None
        
        if self.kd.gamma:
            self.gamma_index = 0
            assert dtype == float
        else:
            for k, k_c in enumerate(self.kd.ibzk_kc):
                if np.all(k_c == 0.):
                    self.gamma_index = k

        assert self.gamma_index is not None
        
        # Matrix of force constants -- dict of dicts in atomic indices
        # Only q-vectors in the irreducible BZ stored here
        self.C_qaavv = [dict([(a, dict([(a_, np.zeros((3, 3), dtype=dtype))
                                        for a_ in self.indices]))
                              for a in self.indices])
                        for q in range(self.kd.nibzkpts)]
        
        # Force constants and dynamical matrix attributes
        # Irreducible zone -- list with array entrances C_avav
        self.C_q = None
        # Full BZ (ndarray)
        self.D_k = None
        
        self.assembled = False

    def __getstate__(self): 
        """Method used to pickle an instance of ``DynamicalMatrix``.

        Bound method attributes cannot be pickled and must therefore be deleted
        before an instance is dumped to file.

        """

        # Get state of object and take care of troublesome attributes
        state = dict(self.__dict__)
        state['kd'].__dict__['comm'] = serial_comm

        return state

    def __setstate__(self, state):
        """Method used to unpickle an instance of ``DynamicalMatrix``."""

        self.__dict__.update(state)

    def get_mass_array(self):
        """Return inverse sqrt of masses (matches shape of mode array)."""

        assert self.m_inv_av is not None
        return self.m_inv_av

    def get_indices(self):
        """Return indices of included atoms."""

        return self.indices
    
    def set_indices(self, indices):
        """Set indices of atoms to be included in the calculation."""
        
        self.indices = indices
        self.C_qaavv = [dict([(a, dict([(a_, np.zeros((3, 3), dtype=self.dtype))
                                        for a_ in self.indices]))
                              for a in self.indices])
                        for q in range(self.kd.nibzkpts)]

    def set_name_and_path(self, name, path):
        """Set name and path of the force constant files.

        name: str
            Base name template for the files which the elements of the matrix
            of force constants are written to.
        path: str
            Path specifying the directory with the files.
            
        """
        
        self.name = name
        self.path = path
        
    def write(self, fd, q, a, v):
        """Write force constants for specified indices to file."""
        
        C_qav_a = [self.C_qaavv[q][a][a_][v] for a_ in self.indices]
        pickle.dump(C_qav_a, fd)
        
    def read(self):
        """Read the force constants from files."""

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
            C += C.conj().T

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
        self.m_inv_av = np.repeat(self.masses[self.indices]**-0.5, 3)
        M_avav = self.m_inv_av[:, np.newaxis] * self.m_inv_av

        for C in self.D_k:
            C *= M_avav

        self.assembled = True
       
    def real_space(self):
        """Fourier transform the dynamical matrix to real-space."""

        if not self.assembled:
            self.assemble()

        # Shape of q-point grid
        N_c = tuple(self.kd.N_c)

        # Reshape before Fourier transforming
        shape = self.D_k.shape
        Dq_lmn = self.D_k.reshape(N_c + shape[1:])
        DR_lmn = fft.ifftn(fft.ifftshift(Dq_lmn, axes=(0, 1, 2)), axes=(0, 1, 2))

        if debug:
            # Check that D_R is real enough
            assert np.all(DR_lmn.imag < 1e-8)
            
        DR_lmn = DR_lmn.real

        # Corresponding R_m vectors in units of the lattice vectors
        R_cm = np.indices(N_c).reshape(3, -1)
        N1_c = np.array(N_c)[:, np.newaxis]        
        R_cm += N1_c // 2
        R_cm %= N1_c
        R_cm -= N1_c // 2
        R_clmn = R_cm.reshape((3,) + N_c)

        return DR_lmn, R_clmn

    def fourier_interpolate(self, N_c):
        """Fourier interpolate dynamical matrix onto a finer q-vector grid."""

        raise NotImplementedError

    def calculate_row(self, perturbation, response_calc):
        """Calculate row of force constant matrix from first-order derivatives.

        Parameters
        ----------
        perturbation: PhononPerturbation
            The perturbation which holds the derivative of the
            pseudo-potential.
        response_calc: ResponseCalculator
            Calculator with the corresponding derivatives of the density and
            the wave-functions.
            
        """

        self.density_derivative(perturbation, response_calc)
        self.wfs_derivative(perturbation, response_calc)
        
    def density_ground_state(self, calc):
        """Contributions involving ground-state density.

        These terms contains second-order derivaties of the localized functions
        ghat and vbar. They are therefore diagonal in the atomic indices.

        """

        # Use the GS LFC's to integrate with the ground-state quantities !
        ghat = calc.density.ghat
        vbar = calc.hamiltonian.vbar
        # Compensation charge coefficients
        Q_aL = calc.density.Q_aL
        
        # Integral of Hartree potential times the second derivative of ghat
        vH_g = calc.hamiltonian.vHt_g
        d2ghat_aLvv = dict([(atom.index, np.zeros((3, 3)))
                            for atom in self.atoms])
        ghat.second_derivative(vH_g, d2ghat_aLvv)

        # Integral of electron density times the second derivative of vbar
        nt_g = calc.density.nt_g
        d2vbar_avv = dict([(atom.index, np.zeros((3, 3)))
                           for atom in self.atoms ])
        vbar.second_derivative(nt_g, d2vbar_avv)

        # Matrix of force constants to be updated; q=-1 for Gamma calculation!
        for C_aavv in self.C_qaavv:

            for a in self.indices:
                # XXX: HGH has only one ghat pr atoms -> generalize when
                # implementing PAW
                C_aavv[a][a] += d2ghat_aLvv[a] * Q_aL[a]
                C_aavv[a][a] += d2vbar_avv[a]

    def wfs_ground_state(self, calc, response_calc):
        """Ground state contributions from the non-local potential."""

        # Projector functions
        pt = calc.wfs.pt
        # Projector coefficients
        dH_asp = calc.hamiltonian.dH_asp
      
        # K-point
        kpt_u = response_calc.wfs.kpt_u
        nbands = response_calc.nbands
        
        for kpt in kpt_u:

            # Index of k
            k = kpt.k
            P_ani = kpt.P_ani
            dP_aniv = kpt.dP_aniv
            
            # Wave functions
            psit_nG = kpt.psit_nG
            psit1_nG = kpt.psit1_nG

            # Calculate d2P_anivv coefficients
            # d2P_anivv = self.calculate_d2P_anivv()
            d2P_anivv = dict([(atom.index,
                               np.zeros(
                (nbands, pt.get_function_count(atom.index), 3, 3)
                )) for atom in self.atoms])
            #XXX Temp dict, second_derivative method only takes a_G array
            # -- no extra dims
            d2P_avv = dict([(atom.index, np.zeros((3, 3)))
                            for atom in self.atoms])
         
            for n in range(nbands):
                pt.second_derivative(psit_nG[n], d2P_avv)
                # Insert in other dict
                for atom in self.atoms:
                    a = atom.index
                    d2P_anivv[a][n, 0] = d2P_avv[a]
            
            for a in self.indices:
    
                H_ii = unpack(dH_asp[a][0])
                P_ni = P_ani[a]
                dP_niv = -1 * dP_aniv[a]
                d2P_nivv = d2P_anivv[a]
                
                # Term with second-order derivative of projector
                HP_ni = np.dot(P_ni, H_ii)
                d2PHP_nvv = (d2P_nivv.conj() *
                             HP_ni[:, :, np.newaxis, np.newaxis]).sum(1)
                d2PHP_nvv *= kpt.weight
                A_vv = d2PHP_nvv.sum(0)
    
                # Term with first-order derivative of the projectors
                HdP_inv = np.dot(H_ii, dP_niv.conj())
                HdP_niv = np.swapaxes(HdP_inv, 0, 1)
                HdP_niv *= kpt.weight
    
                B_vv = (dP_niv[:, :, np.newaxis, :] * 
                        HdP_niv[:, :, :, np.newaxis]).sum(0).sum(0)

                for C_aavv in self.C_qaavv:
                    C_aavv[a][a] += (A_vv + B_vv) + (A_vv + B_vv).conj()

    def core_corrections(self):
        """Contribution from the derivative of the core density."""

        raise NotImplementedError
    
    def density_derivative(self, perturbation, response_calc):
        """Contributions involving the first-order density derivative."""

        # Get attributes from the phononperturbation
        a = perturbation.a
        v = perturbation.v
        #XXX careful here, Gamma calculation has q=-1
        q = perturbation.q

        # Matrix of force constants to be updated; q=-1 for Gamma calculation!
        C_aavv = self.C_qaavv[q]
        
        # Localized functions 
        ghat = perturbation.ghat
        vbar = perturbation.vbar
        # Compensation charge coefficients
        Q_aL = perturbation.Q_aL

        # Density derivative
        nt1_g = response_calc.nt1_g
        
        # Hartree potential derivative including compensation charges
        vH1_g = response_calc.vH1_g.copy()
        vH1_g += perturbation.vghat1_g

        # Integral of Hartree potential derivative times ghat derivative
        dghat_aLv = ghat.dict(derivative=True)
        # Integral of density derivative times vbar derivative
        dvbar_av = vbar.dict(derivative=True)
        
        # Evaluate integrals
        ghat.derivative(vH1_g, dghat_aLv, q=q)
        vbar.derivative(nt1_g, dvbar_av, q=q)

        # Add to force constant matrix attribute
        for a_ in self.indices:
            # Minus sign comes from lfc member function derivative
            C_aavv[a][a_][v] -= np.dot(Q_aL[a_], dghat_aLv[a_])
            C_aavv[a][a_][v] -= dvbar_av[a_][0]

    def wfs_derivative(self, perturbation, response_calc):
        """Contributions from the non-local part of the PAW potential."""

        # Get attributes from the phononperturbation
        a = perturbation.a
        v = perturbation.v
        q = perturbation.q

        # Matrix of force constants to be updated
        C_aavv = self.C_qaavv[q]
           
        # Projector functions
        pt = response_calc.wfs.pt
        # Projector coefficients
        dH_asp = perturbation.dH_asp
        
        # K-point
        kpt_u = response_calc.wfs.kpt_u
        nbands = response_calc.nbands

        # Get k+q indices
        if perturbation.has_q():
            q_c = perturbation.get_q()
            kplusq_k = response_calc.wfs.kd.find_k_plus_q(q_c)
        else:
            kplusq_k = range(len(kpt_u))
            
        for kpt in kpt_u:

            # Indices of k and k+q
            k = kpt.k
            kplusq = kplusq_k[k]

            # Projector coefficients
            P_ani = kpt.P_ani
            dP_aniv = kpt.dP_aniv
            
            psit_nG = kpt.psit_nG
            psit1_nG = kpt.psit1_nG

            # Overlap between wave-function derivative and projectors
            Pdpsi_ani = pt.dict(shape=nbands, zero=True)
            pt.integrate(psit1_nG, Pdpsi_ani, q=kplusq)
            # Overlap between wave-function derivative and derivative of projectors
            dPdpsi_aniv = pt.dict(shape=nbands, derivative=True)
            pt.derivative(psit1_nG, dPdpsi_aniv, q=kplusq)

            for a_ in self.indices:
                # Coefficients from atom a
                Pdpsi_ni = Pdpsi_ani[a]
                dPdpsi_niv = -1 * dPdpsi_aniv[a]
                # Coefficients from atom a_
                H_ii = unpack(dH_asp[a_][0])
                P_ni = P_ani[a_]
                dP_niv = -1 * dP_aniv[a_]
                
                # Term with dPdpsi and P coefficients
                HP_ni = np.dot(P_ni, H_ii)
                dPdpsiHP_nv = (dPdpsi_niv.conj() * HP_ni[:, :, np.newaxis]).sum(1)
                dPdpsiHP_nv *= kpt.weight
                A_v = dPdpsiHP_nv.sum(0)
    
                # Term with dP and Pdpsi coefficients
                HPdpsi_ni = np.dot(Pdpsi_ni.conj(), H_ii)
                dPHPdpsi_nv = (dP_niv * HPdpsi_ni[:, :, np.newaxis]).sum(1)
                dPHPdpsi_nv *= kpt.weight
                B_v = dPHPdpsi_nv.sum(0)

                # Factor of 2 from time-reversal symmetry
                C_aavv[a][a_][v] += 2 * (A_v + B_v)
