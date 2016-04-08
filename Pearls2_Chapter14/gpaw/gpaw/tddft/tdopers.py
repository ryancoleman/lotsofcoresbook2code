# Written by Lauri Lehtovaara, 2007

"""This module implements classes for time-dependent variables and
operators."""

import numpy as np

from gpaw.external_potential import ExternalPotential
from gpaw.utilities import pack2, unpack
from gpaw.mpi import run
from gpaw.fd_operators import Laplace, Gradient
from gpaw.overlap import Overlap
from gpaw.wavefunctions.fd import FDWaveFunctions
from gpaw.tddft.abc import *

# Hamiltonian
class TimeDependentHamiltonian:
    """Time-dependent Hamiltonian, H(t)
    
    This class contains information required to apply time-dependent
    Hamiltonian to a wavefunction.
    """
    
    def __init__(self, wfs, atoms, hamiltonian, td_potential):
        """Create the TimeDependentHamiltonian-object.
        
        The time-dependent potential object must (be None or) have a member
        function strength(self,time), which provides the strength of the
        time-dependent external potential to x-direction at the given time.
        
        Parameters
        ----------
        wfs: FDWaveFunctions
            time-independent grid-based wavefunctions
        hamiltonian: Hamiltonian
            time-independent Hamiltonian
        td_potential: TimeDependentPotential
            time-dependent potential
        """

        self.wfs = wfs
        self.hamiltonian = hamiltonian
        self.td_potential = td_potential
        self.time = self.old_time = 0
        
        # internal smooth potential
        self.vt_sG = hamiltonian.gd.zeros(hamiltonian.nspins)

        # Increase the accuracy of Poisson solver
        if self.hamiltonian.poisson.eps > 1e-12: 
           self.hamiltonian.poisson.eps = 1e-12 

        # external potential
        #if hamiltonian.vext_g is None:
        #    hamiltonian.vext_g = hamiltonian.finegd.zeros()

        #self.ti_vext_g = hamiltonian.vext_g
        #self.td_vext_g = hamiltonian.finegd.zeros(n=hamiltonian.nspins)

        self.P = None

        self.spos_ac = atoms.get_scaled_positions() % 1.0
        self.absorbing_boundary = None
        

    def update(self, density, time):
        """Updates the time-dependent Hamiltonian.
    
        Parameters
        ----------
        density: Density
            the density at the given time  
            (TimeDependentDensity.get_density())
        time: float
            the current time

        """

        self.old_time = self.time = time
        self.hamiltonian.update(density)
        
    def half_update(self, density, time):
        """Updates the time-dependent Hamiltonian, in such a way, that a
        half of the old Hamiltonian is kept and the other half is updated.
        
        Parameters
        ----------
        density: Density
            the density at the given time 
            (TimeDependentDensity.get_density())
        time: float
            the current time

        """
        
        self.old_time = self.time
        self.time = time

        # copy old
        self.vt_sG[:] = self.hamiltonian.vt_sG
        self.dH_asp = {}
        for a, dH_sp in self.hamiltonian.dH_asp.items():
            self.dH_asp[a] = dH_sp.copy()
        # update
        self.hamiltonian.update(density)
        # average and difference
        self.hamiltonian.vt_sG[:], self.vt_sG[:] = \
            0.5*(self.hamiltonian.vt_sG + self.vt_sG), \
            self.hamiltonian.vt_sG - self.vt_sG
        for a, dH_sp in self.hamiltonian.dH_asp.items():
            dH_sp[:], self.dH_asp[a][:] = 0.5*(dH_sp + self.dH_asp[a]), \
                dH_sp - self.dH_asp[a] #pack/unpack is linear for real values

    def half_apply_local_potential(self, psit_nG, Htpsit_nG, s):
        """Apply the half-difference Hamiltonian operator to a set of vectors.
        
        Parameters:

        psit_nG: ndarray
            set of vectors to which the overlap operator is applied.
        psit_nG: ndarray, output
            resulting H applied to psit_nG vectors.
        s: int
            spin index of k-point object defined in kpoint.py.
        
        """
        # Does exactly the same as Hamiltonian.apply_local_potential
        # but uses the difference between vt_sG at time t and t+dt.
        vt_G = self.vt_sG[s]
        if psit_nG.ndim == 3:
            Htpsit_nG += psit_nG * vt_G
        else:
            for psit_G, Htpsit_G in zip(psit_nG, Htpsit_nG):
                Htpsit_G += psit_G * vt_G


    def half_apply(self, kpt, psit, hpsit, calculate_P_ani=True):
        """Applies the half-difference of the time-dependent Hamiltonian
        to the wavefunction psit of the k-point kpt.
        
        Parameters
        ----------
        kpt: Kpoint
            the current k-point (kpt_u[index_of_k-point])
        psit: List of coarse grid
            the wavefuntions (on coarse grid) 
            (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])
        hpsit: List of coarse grid
            the resulting "operated wavefunctions" (H psit)
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | psit> are calculated.
            When False, existing P_uni are used

        """

        hpsit.fill(0.0)
        self.half_apply_local_potential(psit, hpsit, kpt.s)

        # Does exactly the same as last part of Hamiltonian.apply but
        # uses the difference between dH_asp at time t and t+dt.
        shape = psit.shape[:-3]
        P_axi = self.wfs.pt.dict(shape)

        if calculate_P_ani:
            self.wfs.pt.integrate(psit, P_axi, kpt.q)
        else:
            for a, P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        for a, P_xi in P_axi.items():
            dH_ii = unpack(self.dH_asp[a][kpt.s])
            P_axi[a][:] = np.dot(P_xi, dH_ii)
        self.wfs.pt.add(hpsit, P_axi, kpt.q)

        if self.td_potential is not None:
            # FIXME: add half difference here... but maybe it's not important
            # as this will be used only for getting initial guess. So, should
            # not affect to the results, only to the speed of convergence.
            #raise NotImplementedError
            pass

    def apply(self, kpt, psit, hpsit, calculate_P_ani=True):
        """Applies the time-dependent Hamiltonian to the wavefunction psit of
        the k-point kpt.
        
        Parameters
        ----------
        kpt: Kpoint
            the current k-point (kpt_u[index_of_k-point])
        psit: List of coarse grid
            the wavefuntions (on coarse grid) 
            (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])
        hpsit: List of coarse grid
            the resulting "operated wavefunctions" (H psit)
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | psit> are calculated.
            When False, existing P_uni are used

        """

        self.hamiltonian.apply(psit, hpsit, self.wfs, kpt, calculate_P_ani)

        # PAW correction
        if self.P is not None:
            self.P.add(psit, hpsit, self.wfs, kpt)

        # Absorbing boundary conditions

        # Imaginary potential
        if self.absorbing_boundary is not None \
               and self.absorbing_boundary.type == 'IPOT':
            hpsit[:] += self.absorbing_boundary.get_potential_matrix() * psit

        # Perfectly matched layers
        if self.absorbing_boundary is not None \
               and self.absorbing_boundary.type == 'PML':
            # Perfectly matched layer is applied as potential Vpml = Tpml-T
            # Where  T = -0.5*\nabla^{2}\psi  (Use latex for these equations)
            # See abc.py for details
            # This is probably not the most optimal approach and slows
            # the propagation.
            if self.lpsit is None:
                self.lpsit = self.hamiltonian.gd.empty( n=len(psit),
                                                        dtype=complex )
            self.laplace.apply(psit, self.lpsit, kpt.phase_cd)
            hpsit[:] -= (.5 * (self.absorbing_boundary.get_G()**2 - 1.0)
                         * self.lpsit)
            for i in range(3):
                self.gradient[i].apply(psit, self.lpsit, kpt.phase_cd)
                hpsit[:] -= (.5 * self.absorbing_boundary.get_G()
                             * self.absorbing_boundary.get_dG()[i]
                             * self.lpsit)


        # Time-dependent dipole field
        if self.td_potential is not None:
            #TODO on shaky ground here...
            strength = self.td_potential.strength
            ExternalPotential().add_linear_field(self.wfs, self.spos_ac,
                                                 psit, hpsit,
                                                 0.5 * strength(self.time) +
                                                 0.5 * strength(self.old_time),
                                                 kpt)

            
    def set_absorbing_boundary(self, absorbing_boundary):
        """ Sets up the absorbing boundary.            
            Parameters:
            absorbing_boundary: absorbing boundary object of any kind.  
        """
        
        self.absorbing_boundary = absorbing_boundary
        self.absorbing_boundary.set_up(self.hamiltonian.gd)
        if self.absorbing_boundary.type == 'PML':
            gd = self.hamiltonian.gd
            self.laplace = Laplace(gd, n=2, dtype=complex)
            self.gradient = np.array((Gradient(gd,0, n=2, dtype=complex),
                                       Gradient(gd,1, n=2, dtype=complex),
                                       Gradient(gd,2, n=2, dtype=complex)))
            self.lpsit=None

    def calculate_paw_correction(self, psit_nG, hpsit, wfs, kpt, v_atom, calculate_P_ani=True):
        """ Operates on psit_nG with the P-term that is present in PAW based Ehrenfest dynamics

        Parameters
        ----------
        psit_nG: list of coarse grid wavefunctions

        hpsit: array to which P psit_nG is added

        wfs: Wavefunctions object

        kpt: Kpoint object

        v_atom: atomic velocities

        """

        shape = psit_nG.shape[:-3]
        P_axi = wfs.pt.dict(shape)
        wfs.pt.integrate(psit_nG, P_axi, kpt.q)
        
        #G_LLL = gaunt # G_LLL[L1,L2,L3] = \int Y_L1 Y_L2 Y_L3
            
        #Coefficients for calculating P \psi_n
        # P = -i sum_a v_a P^a, P^a = T^{\dagger} \nabla_{R_a} T
        w_ani = wfs.pt.dict(wfs.bd.mynbands, zero=True)
        #projector derivatives < nabla pt_i^a | psit_n >
        dpt_aniv = wfs.pt.dict(wfs.bd.mynbands, derivative=True)     
        wfs.pt.derivative(psit_nG, dpt_aniv, kpt.q)
        #wfs.calculate_forces(paw.hamiltonian, F_av)
        for a in dpt_aniv.keys():
            #ni = wfs.pt.get_function_count(a)
            for c in range(3):

                P_xi = P_axi[a]
                #nabla_iiv contains terms < \phi_i1^a | d / d v phi_i2^a >
                #- < phit_i1^a | d / dv phit_i2^a>, where v is either x,y or z              
                nabla_ii = wfs.setups[a].nabla_iiv[:,:,c]
                dpt_ni = dpt_aniv[a][:,:,c]
                dO_ii = wfs.setups[a].dO_ii
                #dphi_aniv[a] = np.dot(P_xi, nabla_ii.transpose())
                dphi_ni = np.dot(P_xi, nabla_ii.transpose())
                pt_ni = np.dot(dpt_ni, dO_ii)
                #pt_aniv[a] = np.dot(Dpt_ni, dO_ii)
                w_ani[a] += (dphi_ni + pt_ni) * v_atom[a,c]
                        
            w_ani[a] *= complex(0,1)
            #dO_ani[a] *= complex(0,1)

        #wfs.pt.add(ppsit, W_ani, kpt.q)
        wfs.pt.add(hpsit, w_ani, kpt.q)


# AbsorptionKickHamiltonian
class AbsorptionKickHamiltonian:
    """Absorption kick Hamiltonian, p.r
    
    This class contains information required to apply absorption kick
    Hamiltonian to a wavefunction.
    """
    
    def __init__(self, wfs, atoms, strength=[0.0, 0.0, 1e-3]):
        """Create the AbsorptionKickHamiltonian-object.

        Parameters
        ----------
        wfs: FDWaveFunctions
            time-independent grid-based wavefunctions
        atoms: Atoms
            list of atoms
        strength: float[3]
            strength of the delta field to different directions

        """

        self.wfs = wfs
        self.spos_ac = atoms.get_scaled_positions() % 1.0
        
        # magnitude
        magnitude = np.sqrt(strength[0]*strength[0] 
                             + strength[1]*strength[1] 
                             + strength[2]*strength[2])
        # iterations
        self.iterations = int(round(magnitude / 1.0e-4))
        if self.iterations < 1:
            self.iterations = 1
        # delta p
        self.dp = strength / self.iterations

        # hamiltonian
        self.abs_hamiltonian = np.array([self.dp[0], self.dp[1], self.dp[2]])
        

    def update(self, density, time):
        """Dummy function = does nothing. Required to have correct interface.
        
        Parameters
        ----------
        density: Density or None
            the density at the given time or None (ignored)
        time: Float or None
            the current time (ignored)

        """
        pass
        
    def half_update(self, density, time):
        """Dummy function = does nothing. Required to have correct interface.
        
        Parameters
        ----------
        density: Density or None
            the density at the given time or None (ignored)
        time: float or None
            the current time (ignored)

        """
        pass
        
    def apply(self, kpt, psit, hpsit, calculate_P_ani=True):
        """Applies the absorption kick Hamiltonian to the wavefunction psit of
        the k-point kpt.
        
        Parameters
        ----------
        kpt: Kpoint
            the current k-point (kpt_u[index_of_k-point])
        psit: List of coarse grids
            the wavefuntions (on coarse grid) 
            (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])
        hpsit: List of coarse grids
            the resulting "operated wavefunctions" (H psit)
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | psit> are calculated.
            When False, existing P_uni are used

        """
        hpsit[:] = 0.0

        #TODO on shaky ground here...
        ExternalPotential().add_linear_field(self.wfs, self.spos_ac,
                                             psit, hpsit,
                                             self.abs_hamiltonian, kpt)


# Overlap
class TimeDependentOverlap(Overlap):
    """Time-dependent overlap operator S(t)
    
    This class contains information required to apply time-dependent
    overlap operator to a set of wavefunctions.
    """
    
    def __init__(self, ksl, timer):
        """Creates the TimeDependentOverlap-object.
        
        Parameters
        ----------
        XXX TODO

        """
        Overlap.__init__(self, ksl, timer)

    def update_k_point_projections(self, wfs, kpt, psit=None):
        """Updates the projector function overlap integrals
        with the wavefunctions of a given k-point.
        
        Parameters
        ----------
        wfs: TimeDependentWaveFunctions
            time-independent grid-based wavefunctions
        kpt: Kpoint
            the current k-point (kpt_u[index_of_k-point])
        psit: List of coarse grids (optional)
            the wavefuntions (on coarse grid) 
            (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])

        """
        if psit is not None:
            wfs.pt.integrate(psit, kpt.P_ani, kpt.q)
        else:
            wfs.pt.integrate(kpt.psit_nG, kpt.P_ani, kpt.q)

    def update(self, wfs):
        """Updates the time-dependent overlap operator.
        
        Parameters
        ----------
        wfs: TimeDependentWaveFunctions
            time-independent grid-based wavefunctions

        """
        for kpt in wfs.kpt_u:
            self.update_k_point_projections(wfs, kpt)
    
    def half_update(self, wfs):
        """Updates the time-dependent overlap operator, in such a way,
        that a half of the old overlap operator is kept and the other half
        is updated. !Currently does nothing!

        Parameters
        ----------
        wfs: TimeDependentWaveFunctions
            time-independent grid-based wavefunctions

        """
        #for kpt in wfs.kpt_u:
        #    # copy old
        #    P_ani = {}
        #    for a,P_ni in kpt.P_ani.items():
        #        P_ani[a] = P_ni.copy()
        #    # update
        #    self.update_k_point_projections(wfs, kpt)
        #    # average
        #    for a,P_ni in P_ani.items():
        #        kpt.P_ani[a] += P_ni
        #        kpt.P_ani[a] *= .5

        # !!! FIX ME !!! update overlap operator/projectors/...
        pass
    
    #def apply(self, psit, spsit, wfs, kpt, calculate_P_ani=True):
    #    """Apply the time-dependent overlap operator to the wavefunction
    #    psit of the k-point kpt.
    #    
    #    Parameters
    #    ----------
    #    psit: List of coarse grids
    #        the wavefuntions (on coarse grid) 
    #        (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])
    #    spsit: List of coarse grids
    #        the resulting "operated wavefunctions" (S psit)
    #    wfs: TimeDependentWaveFunctions
    #        time-independent grid-based wavefunctions
    #    kpt: Kpoint
    #        the current k-point (kpt_u[index_of_k-point])
    #    calculate_P_ani: bool
    #        When True, the integrals of projector times vectors
    #        P_ni = <p_i | psit> are calculated.
    #        When False, existing P_ani are used
    #
    #    """
    #    self.overlap.apply(psit, spsit, wfs, kpt, calculate_P_ani)
    #
    def apply_inverse(self, a_nG, b_nG, wfs, kpt, calculate_P_ani=True, use_cg=True):
        """Apply the approximative time-dependent inverse overlap operator
        to the wavefunction psit of the k-point kpt.
    
        Parameters
        ----------
        a_nG: List of coarse grids
            the wavefuntions (on coarse grid) 
            (kpt_u[index_of_k-point].psit_nG[indices_of_wavefunc])
        b_nG: List of coarse grids
            the resulting "operated wavefunctions" (S^(-1) psit)
        wfs: TimeDependentWaveFunctions
            time-independent grid-based wavefunctions
        kpt: Kpoint
            the current k-point (kpt_u[index_of_k-point])
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | psit> are calculated.
            When False, existing P_uni are used
        use_cg: bool
            When True, use conjugate gradient method to solve for inverse.
    
        """
        if not use_cg:
            self.timer.start('Apply approximate inverse overlap')
            Overlap.apply_inverse(self, a_nG, b_nG, wfs, kpt, calculate_P_ani)
            self.timer.stop('Apply approximate inverse overlap')
            return            

        self.timer.start('Apply exact inverse overlap')
        from gpaw.utilities.blas import dotu, axpy, dotc
        #from gpaw.tddft.cscg import multi_zdotu, multi_scale, multi_zaxpy
        #initialization
          # Multivector dot product, a^T b, where ^T is transpose
        def multi_zdotu(s, x,y, nvec):
            for i in range(nvec):
                s[i] = dotu(x[i],y[i])
            wfs.gd.comm.sum(s)
            return s
        # Multivector ZAXPY: a x + y => y
        def multi_zaxpy(a,x,y, nvec):
            for i in range(nvec):
                axpy(a[i]*(1+0J), x[i], y[i])
        # Multiscale: a x => x
        def multi_scale(a,x, nvec):
            for i in range(nvec):
                x[i] *= a[i]
        nvec = len(a_nG)
        r = wfs.gd.zeros(nvec, dtype=wfs.dtype)
        z  = wfs.gd.zeros((nvec,), dtype=wfs.dtype)
        sx = wfs.gd.zeros(nvec, dtype=wfs.dtype)
        p = wfs.gd.zeros(nvec, dtype=wfs.dtype)
        q = wfs.gd.zeros(nvec, dtype=wfs.dtype)
        alpha = np.zeros((nvec,), dtype=wfs.dtype)
        beta = np.zeros((nvec,), dtype=wfs.dtype)
        scale = np.zeros((nvec,), dtype=wfs.dtype)
        normr2 = np.zeros((nvec,), dtype=wfs.dtype)
        rho  = np.zeros((nvec,), dtype=wfs.dtype) 
        rho_prev  = np.zeros((nvec,), dtype=wfs.dtype)
        rho_prev[:] = 1.0
        tol_cg = 1e-14
        multi_zdotu(scale, a_nG, a_nG, nvec)
        scale = np.abs(scale)

        x = b_nG #XXX TODO rename this

        #x0 = S^-1_approx a_nG
        self.apply_inverse(a_nG, x, wfs, kpt, calculate_P_ani, use_cg=False)
        #r0 = a_nG - S x_0
        self.apply(-x, r, wfs, kpt, calculate_P_ani)
        r += a_nG
        #print 'r.max() =', abs(r).max()

        max_iter = 50

        for i in range(max_iter):

            #print 'iter =', i

            self.apply_inverse(r, z, wfs, kpt, calculate_P_ani, use_cg=False)

            multi_zdotu(rho, r, z, nvec)

            beta = rho / rho_prev
            multi_scale(beta, p, nvec)
            p += z

            self.apply(p,q,wfs,kpt, calculate_P_ani)

            multi_zdotu(alpha, p, q, nvec)
            alpha = rho/alpha

            multi_zaxpy(alpha, p, x, nvec)
            multi_zaxpy(-alpha, q, r, nvec)

            multi_zdotu(normr2, r,r, nvec)

            #rhoc = rho.copy()
            rho_prev[:] = rho.copy()
            #rho.copy()
            #rho_prev = rho.copy()

            #print '||r|| =', np.sqrt(np.abs(normr2/scale))
            if ( (np.sqrt(np.abs(normr2) / scale)) < tol_cg ).all():
                break

        self.timer.stop('Apply exact inverse overlap')


class TimeDependentWaveFunctions(FDWaveFunctions):
    def __init__(self, stencil, diagksl, orthoksl, initksl, gd, nvalence, setups,
                 bd, world, kd, kptband_comm, timer=None):
        FDWaveFunctions.__init__(self, stencil, diagksl, orthoksl, initksl,
                                 gd, nvalence, setups, bd, complex, world,
                                 kd, kptband_comm, timer)

    def make_overlap(self):
        return TimeDependentOverlap(self.orthoksl, self.timer)

    def calculate_forces(self, hamiltonian, F_av):
        """ Calculate wavefunction forces with optional corrections for
            Ehrenfest dynamics
        """  

            
        #If td_correction is not none, we replace the overlap part of the
        #force, sum_n f_n eps_n < psit_n | dO / dR_a | psit_n>, with
        #sum_n f_n <psit_n | H S^-1 D^a + c.c. | psit_n >, with D^a
        #defined as D^a = sum_{i1,i2} | pt_i1^a > [O_{i1,i2} < d pt_i2^a / dR_a |
        #+ (< phi_i1^a | d phi_i2^a / dR_a > - < phit_i1^a | d phit_i2^a / dR_a >) < pt_i1^a|].
        #This is required in order to conserve the total energy also when electronic
        #excitations start to play a significant role.

        #TODO: move the corrections into the tddft directory

        # Calculate force-contribution from k-points:
        F_av.fill(0.0)
        F_aniv = self.pt.dict(self.bd.mynbands, derivative=True)
        #print 'self.dtype =', self.dtype
        for kpt in self.kpt_u:
            self.pt.derivative(kpt.psit_nG, F_aniv, kpt.q)

            #self.overlap.update_k_point_projections(kpt)
            self.pt.integrate(kpt.psit_nG, kpt.P_ani, kpt.q)
            hpsit = self.gd.zeros(len(kpt.psit_nG), dtype=self.dtype)
            #eps_psit = self.gd.zeros(len(kpt.psit_nG), dtype=self.dtype)
            sinvhpsit = self.gd.zeros(len(kpt.psit_nG), dtype=self.dtype)
            hamiltonian.apply(kpt.psit_nG, hpsit, self, kpt, calculate_P_ani=True)
            self.overlap.apply_inverse(hpsit, sinvhpsit, self, kpt, calculate_P_ani=True)
            #print 'sinvhpsit_0_cg - epspsit_0.max', abs(sinvhpsit[0]-eps_psit[0]).max()
            #print 'sinvhpsit_0 - epspsit_0.max', abs(sinvhpsit2[0]-eps_psit[0]).max()

            G_axi = self.pt.dict(self.bd.mynbands)
            self.pt.integrate(sinvhpsit, G_axi, kpt.q)

            for a, F_niv in F_aniv.items():
                F_niv = F_niv.conj()
                F_niv *= kpt.f_n[:, np.newaxis, np.newaxis]
                FdH1_niv = F_niv.copy()
                dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                P_ni = kpt.P_ani[a]
                dO_ii = hamiltonian.setups[a].dO_ii
                F_vii = np.dot(np.dot(F_niv.transpose(), P_ni), dH_ii)

                fP_ni = P_ni * kpt.f_n[:,np.newaxis]
                G_ni = G_axi[a]
                nabla_iiv = hamiltonian.setups[a].nabla_iiv
                F_vii_sinvh_dpt = -np.dot(np.dot(FdH1_niv.transpose(), G_ni), dO_ii)
                F_vii_sinvh_dphi = -np.dot(nabla_iiv.transpose(2,0,1), np.dot(fP_ni.conj().transpose(), G_ni))
                F_vii += F_vii_sinvh_dpt + F_vii_sinvh_dphi
                                    
                #F_av_dO[a] += 2 * F_vii_dO.real.trace(0,1,2)
                #F_av_dH[a] += 2 * F_vii_dH.real.trace(0,1,2)
                F_av[a] += 2 * F_vii.real.trace(0, 1, 2)

            # Hack used in delta-scf calculations:
            if hasattr(kpt, 'c_on'):
                assert self.bd.comm.size == 1
                self.pt.derivative(kpt.psit_nG, F_aniv, kpt.q)  #XXX again
                d_nn = np.zeros((self.bd.mynbands, self.bd.mynbands),
                                dtype=complex)
                for ne, c_n in zip(kpt.ne_o, kpt.c_on):
                    d_nn += ne * np.outer(c_n.conj(), c_n)
                for a, F_niv in F_aniv.items():
                    F_niv = F_niv.conj()
                    dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                    Q_ni = np.dot(d_nn, kpt.P_ani[a])
                    F_vii = np.dot(np.dot(F_niv.transpose(), Q_ni), dH_ii)
                    F_niv *= kpt.eps_n[:, np.newaxis, np.newaxis]
                    dO_ii = hamiltonian.setups[a].dO_ii
                    F_vii -= np.dot(np.dot(F_niv.transpose(), Q_ni), dO_ii)
                    F_av[a] += 2 * F_vii.real.trace(0, 1, 2)

        self.bd.comm.sum(F_av, 0)

        if self.bd.comm.rank == 0:
            self.kd.comm.sum(F_av, 0)


# DummyDensity
class DummyDensity:
    """Implements dummy (= does nothing) density for AbsorptionKick."""

    def __init__(self, wfs):
        """Placeholder Density object for AbsorptionKick.

        Parameters
        ----------
        wfs: FDWaveFunctions
            time-independent grid-based wavefunctions

        """
        self.wfs = wfs

    def update(self):
        pass

    def get_wavefunctions(self):
        return self.wfs

    def get_density(self):
        return None


# Density
class TimeDependentDensity(DummyDensity):
    """Time-dependent density rho(t)
    
    This class contains information required to get the time-dependent
    density.
    """
    
    def __init__(self, paw):
        """Creates the TimeDependentDensity-object.
        
        Parameters
        ----------
        paw: PAW
            the PAW-object
        """
        DummyDensity.__init__(self, paw.wfs)
        self.density = paw.density

    def update(self):
        """Updates the time-dependent density.
        
        Parameters
        ----------
        None

        """
        #for kpt in self.wfs.kpt_u:
        #    self.wfs.pt.integrate(kpt.psit_nG, kpt.P_ani)
        self.density.update(self.wfs)
       
    def get_density(self):
        """Returns the current density.
        
        Parameters
        ----------
        None

        """
        return self.density
