# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
Ref. to Kresse-paper ... XXX
"""

import numpy as np

from gpaw.utilities.blas import axpy
from gpaw.fd_operators import FDOperator


class BaseMixer:
    """Pulay density mixer."""
    
    def __init__(self, beta=0.1, nmaxold=3, weight=50.0, dtype=float):
        """Construct density-mixer object.

        Parameters
        ----------
        beta: float
            Mixing parameter between zero and one (one is most
            aggressive).
        nmaxold: int
            Maximum number of old densities.
        weight: float
            Weight parameter for special metric (for long wave-length
            changes).

        """

        self.beta = beta
        self.nmaxold = nmaxold
        self.weight = weight
        self.dtype = dtype

        self.dNt = None

        self.mix_rho = False

    def initialize_metric(self, gd):
        self.gd = gd

        if self.weight == 1:
            self.metric = None
        else:
            a = 0.125 * (self.weight + 7)
            b = 0.0625 * (self.weight - 1)
            c = 0.03125 * (self.weight - 1)
            d = 0.015625 * (self.weight - 1)
            self.metric = FDOperator([a,
                                      b, b, b, b, b, b,
                                      c, c, c, c, c, c, c, c, c, c, c, c,
                                      d, d, d, d, d, d, d, d],
                                     [(0, 0, 0),
                                      (-1, 0, 0), (1, 0, 0),                 #b
                                      (0, -1, 0), (0, 1, 0),                 #b
                                      (0, 0, -1), (0, 0, 1),                 #b
                                      (1, 1, 0), (1, 0, 1), (0, 1, 1),       #c
                                      (1, -1, 0), (1, 0, -1), (0, 1, -1),    #c
                                      (-1, 1, 0), (-1, 0, 1), (0, -1, 1),    #c
                                      (-1, -1, 0), (-1, 0, -1), (0, -1, -1), #c
                                      (1, 1, 1), (1, 1, -1), (1, -1, 1),     #d
                                      (-1, 1, 1), (1, -1, -1), (-1, -1, 1),  #d
                                      (-1, 1, -1), (-1, -1, -1)              #d
                                      ],
                                     gd, self.dtype).apply
            self.mR_G = gd.empty(dtype=self.dtype)
        
    def initialize(self, density):
        self.initialize_metric(density.gd)

    def reset(self):
        """Reset Density-history.

        Called at initialization and after each move of the atoms.

        my_nuclei:   All nuclei in local domain.
        """
        
        # History for Pulay mixing of densities:
        self.nt_iG = [] # Pseudo-electron densities
        self.R_iG = []  # Residuals
        self.A_ii = np.zeros((0, 0))
        self.dNt = None
        
        self.D_iap = []
        self.dD_iap = []

    def get_charge_sloshing(self):
        """Return number of electrons moving around.

        Calculated as the integral of the absolute value of the change
        of the density from input to output."""
        
        return self.dNt

    def set_charge_sloshing(self, dNt):
        self.dNt = dNt
        
    def mix(self, nt_G, D_ap, phase_cd=None):
        iold = len(self.nt_iG)
        if iold > 0:
            if iold > self.nmaxold:
                # Throw away too old stuff:
                del self.nt_iG[0]
                del self.R_iG[0]
                del self.D_iap[0]
                del self.dD_iap[0]
                # for D_p, D_ip, dD_ip in self.D_a:
                #     del D_ip[0]
                #     del dD_ip[0]
                iold = self.nmaxold

            # Calculate new residual (difference between input and output)
            R_G = nt_G - self.nt_iG[-1]
            # Use np.absolute instead of np.fabs
            self.dNt = self.gd.integrate(np.absolute(R_G))
            self.R_iG.append(R_G)
            self.dD_iap.append([])
            for D_p, D_ip in zip(D_ap, self.D_iap[-1]):
                self.dD_iap[-1].append(D_p - D_ip)

            # Update matrix:
            A_ii = np.zeros((iold, iold))
            i1 = 0
            i2 = iold - 1
            
            if self.metric is None:
                mR_G = R_G
            else:
                mR_G = self.mR_G
                self.metric(R_G, mR_G, phase_cd=phase_cd)
                
            for R_1G in self.R_iG:
                # Inner product between new and old residues
                # XXX For now, use only real part of residues
                # For complex quantities a .conjugate should be added ??
                a = self.gd.comm.sum(np.vdot(R_1G.real, mR_G.real))
                if self.dtype == complex:
                    a += self.gd.comm.sum(np.vdot(R_1G.imag, mR_G.imag))
                    
                A_ii[i1, i2] = a
                A_ii[i2, i1] = a
                i1 += 1
            A_ii[:i2, :i2] = self.A_ii[-i2:, -i2:]
            self.A_ii = A_ii

            try:
                B_ii = np.linalg.inv(A_ii)
            except np.linalg.LinAlgError:
                alpha_i = np.zeros(iold)
                alpha_i[-1] = 1.0
            else:
                alpha_i = B_ii.sum(1)
                try:
                    # Normalize:
                    alpha_i /= alpha_i.sum()
                except ZeroDivisionError:
                    alpha_i[:] = 0.0
                    alpha_i[-1] = 1.0

            # Calculate new input density:
            nt_G[:] = 0.0

            for D in D_ap:
                D[:] = 0.0
            beta = self.beta

            for i, alpha in enumerate(alpha_i):
                axpy(alpha, self.nt_iG[i], nt_G)
                axpy(alpha * beta, self.R_iG[i], nt_G)
                for D_p, D_ip, dD_ip in zip(D_ap, self.D_iap[i],
                                            self.dD_iap[i]):
                    axpy(alpha, D_ip, D_p)
                    axpy(alpha * beta, dD_ip, D_p)

        # Store new input density (and new atomic density matrices):
        self.nt_iG.append(nt_G.copy())
        self.D_iap.append([])
        for D_p in D_ap:
            self.D_iap[-1].append(D_p.copy())

    def estimate_memory(self, mem, gd):
        gridbytes = gd.bytecount()
        mem.subnode('nt_iG, R_iG', 2 * self.nmaxold * gridbytes)

    def __repr__(self):
        classname = self.__class__.__name__
        template = '%s(beta=%f, nmaxold=%d, weight=%f)'
        string = template % (classname, self.beta, self.nmaxold, self.weight)
        return string


class DummyMixer(BaseMixer):
    """Dummy mixer for TDDFT, i.e., it does not mix."""
    def mix(self, nt_G):
        pass

    def estimate_memory(self, mem, gd):
        pass


class Mixer(BaseMixer):
    """Mix spin up and down densities separately"""

    def initialize(self, density):
        self.mixers = []
        for s in range(density.nspins):
            mixer = BaseMixer(self.beta, self.nmaxold, self.weight)
            mixer.initialize_metric(density.gd)
            self.mixers.append(mixer)
    
    def mix(self, density):
        """Mix pseudo electron densities."""

        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()
        D_sap = []
        for s in range(density.nspins):
            D_sap.append([D_sp[s] for D_sp in D_asp])
        for nt_G, D_ap, mixer in zip(nt_sG, D_sap, self.mixers):
            mixer.mix(nt_G, D_ap)

    def reset(self):
        for mixer in self.mixers:
            mixer.reset()

    def get_charge_sloshing(self):
        """Return number of electrons moving around.

        Calculated as the integral of the absolute value of the change
        of the density from input to output."""

        if self.mixers[0].dNt is None:
            return None
        return sum([mixer.dNt for mixer in self.mixers])

    def set_charge_sloshing(self, dNt):
        for mixer in self.mixers:
            mixer.set_charge_sloshing(dNt / len(self.mixers))


class MixerSum(BaseMixer):
    """For pseudo electron densities, mix the total charge density and for
    density matrices, mix spin up and densities separately.
    Magnetization density is not mixed, i.e new magnetization density is used
    """

    def mix(self, density):
        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        BaseMixer.mix(self, nt_G, D_asp)

        # Only new magnetization for spin density
        dnt_G = nt_sG[0] - nt_sG[1]
        #dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]

        # Construct new spin up/down densities 
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)


class MixerSum2(BaseMixer):
    """Mix the total pseudo electron density and the total density matrices.
    Magnetization density is not mixed, i.e new magnetization density is used.
    """

    def mix(self, density):

        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        D_ap = [D_p[0] + D_p[1] for D_p in D_asp]
        BaseMixer.mix(self, nt_G, D_ap)

        # Only new magnetization for spin density
        dnt_G = nt_sG[0] - nt_sG[1]
        dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]

        # Construct new spin up/down densities 
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)
        for D_sp, D_p, dD_p in zip(D_asp, D_ap, dD_ap):
            D_sp[0] = 0.5 * (D_p + dD_p)
            D_sp[1] = 0.5 * (D_p - dD_p)


class MixerDif(BaseMixer):
    """Mix the charge density and magnetization density separately"""
    
    def __init__(self, beta=0.1, nmaxold=3, weight=50.0,
                 beta_m=0.7, nmaxold_m=2, weight_m=10.0):
        """Construct density-mixer object.

        Parameters
        ----------
        beta: float
            Mixing parameter between zero and one (one is most
            aggressive).
        nmaxold: int
            Maximum number of old densities.
        weight: float
            Weight parameter for special metric (for long wave-length
            changes).
        """

        self.beta = beta
        self.nmaxold = nmaxold
        self.weight = weight

        self.beta_m = beta_m
        self.nmaxold_m = nmaxold_m
        self.weight_m = weight_m
        self.dNt = None

        self.mix_rho = False


    def initialize(self, density):
        assert density.nspins == 2
        self.mixer = BaseMixer(self.beta, self.nmaxold, self.weight)
        self.mixer.initialize_metric(density.gd)
        self.mixer_m = BaseMixer(self.beta_m, self.nmaxold_m, self.weight_m)
        self.mixer_m.initialize_metric(density.gd)

    def reset(self):
        self.mixer.reset()
        self.mixer_m.reset()

    def mix(self, density):

        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        D_ap = [D_sp[0] + D_sp[1] for D_sp in D_asp]
        self.mixer.mix(nt_G, D_ap)

        # Mix magnetization
        dnt_G = nt_sG[0] - nt_sG[1]
        dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]
        self.mixer_m.mix(dnt_G, dD_ap)

        # Construct new spin up/down densities 
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)
        for D_sp, D_p, dD_p in zip(D_asp, D_ap, dD_ap):
            D_sp[0] = 0.5 * (D_p + dD_p)
            D_sp[1] = 0.5 * (D_p - dD_p)
            

    def get_charge_sloshing(self):
        if self.mixer.dNt is None:
            return None
        return self.mixer.dNt


class MixerRho(BaseMixer):
    def initialize(self, density):
    
        self.mix_rho = True
        self.initialize_metric(density.finegd)
    
    def mix(self, density):
        """Mix pseudo electron densities."""

        rhot_g = density.rhot_g
        BaseMixer.mix(self, rhot_g, [])


class MixerRho2(BaseMixer):
    def initialize(self, density):
    
        self.mix_rho = True
        self.initialize_metric(density.finegd)
    
    def mix(self, density):
        """Mix pseudo electron densities."""

        rhot_g = density.rhot_g
        BaseMixer.mix(self, rhot_g, density.D_asp.values())


class BaseMixer_Broydn:
    def __init__(self, beta=0.1, nmaxold=6):
        self.verbose = False
        self.beta = beta
        self.nmaxold = nmaxold
        self.weight = 1
        self.mix_rho = False
        
    def initialize(self, density):
        self.gd = density.gd

    def reset(self):
        self.step = 0
        self.d_nt_G = []
        self.d_D_ap = []
        self.nt_iG = []
        self.D_iap = []
        self.c_G =  []
        self.v_G = []
        self.u_G = []
        self.u_D = []        
        self.dNt = None
        
    def get_charge_sloshing(self):
        return self.dNt 
    
    def mix(self, nt_G, D_ap):
        if self.step > 2:
            del self.d_nt_G[0]
            for d_Dp in self.d_D_ap:
                del d_Dp[0]
        if self.step > 0:
            self.d_nt_G.append(nt_G - self.nt_iG[-1])
            for d_Dp, D_p, D_ip in zip(self.d_D_ap, D_ap, self.D_iap):
                d_Dp.append(D_p - D_ip[-1])
            fmin_G = self.gd.integrate(self.d_nt_G[-1] * self.d_nt_G[-1])
            self.dNt = self.gd.integrate(np.fabs(self.d_nt_G[-1]))
            if self.verbose:
                print('Mixer: broydn: fmin_G = %f fmin_D = %f'% fmin_G)
        if self.step == 0:
            self.eta_G = np.empty(nt_G.shape)
            self.eta_D = []
            for D_p in D_ap:
                self.eta_D.append(0)
                self.u_D.append([])
                self.D_iap.append([])
                self.d_D_ap.append([])
        else:
            if self.step >= 2:
                del self.c_G[:]
                if len(self.v_G) >= self.nmaxold:
                    del self.u_G[0]
                    del self.v_G[0]
                    for u_D in self.u_D:
                        del u_D[0]
                temp_nt_G = self.d_nt_G[1] - self.d_nt_G[0]
                self.v_G.append(temp_nt_G / self.gd.integrate(temp_nt_G
                                                                 * temp_nt_G))
                if len(self.v_G) < self.nmaxold:
                    nstep = self.step - 1
                else:
                    nstep = self.nmaxold 
                for i in range(nstep):
                    self.c_G.append(self.gd.integrate(self.v_G[i] *
                                                             self.d_nt_G[1]))
                self.u_G.append(self.beta  * temp_nt_G + self.nt_iG[1] - self.nt_iG[0])
                for d_Dp, u_D, D_ip in zip(self.d_D_ap, self.u_D, self.D_iap):
                    temp_D_ap = d_Dp[1] - d_Dp[0]
                    u_D.append(self.beta * temp_D_ap + D_ip[1] - D_ip[0])
                usize = len(self.u_G)
                for i in range(usize - 1):
                    a_G = self.gd.integrate(self.v_G[i] * temp_nt_G)
                    axpy(-a_G, self.u_G[i], self.u_G[usize - 1])
                    for u_D in self.u_D:
                        axpy(-a_G, u_D[i], u_D[usize - 1])
            self.eta_G = self.beta * self.d_nt_G[-1]
            for i, d_Dp in enumerate(self.d_D_ap):
                self.eta_D[i] = self.beta * d_Dp[-1]
            usize = len(self.u_G) 
            for i in range(usize):
                axpy(-self.c_G[i], self.u_G[i], self.eta_G)
                for eta_D, u_D in zip(self.eta_D, self.u_D):
                    axpy(-self.c_G[i], u_D[i], eta_D)
            axpy(-1.0, self.d_nt_G[-1], nt_G)
            axpy(1.0, self.eta_G, nt_G)
            for D_p, d_Dp, eta_D in zip(D_ap, self.d_D_ap, self.eta_D):            
                axpy(-1.0, d_Dp[-1], D_p)
                axpy(1.0, eta_D, D_p)
            if self.step >= 2:
                del self.nt_iG[0]
                for D_ip in self.D_iap:
                    del D_ip[0]
        self.nt_iG.append(np.copy(nt_G))
        for D_ip, D_p in zip(self.D_iap, D_ap):
            D_ip.append(np.copy(D_p))
        self.step += 1

        
class Mixer_Broydn(BaseMixer_Broydn):
    """Mix spin up and down densities separately"""

    def initialize(self, density):
        self.mixers = []
        for s in range(density.nspins):
            mixer = BaseMixer_Broydn()
            mixer.initialize(density)
            #mixer.initialize_metric(density.gd)
            self.mixers.append(mixer)
    
    def mix(self, density):
        """Mix pseudo electron densities."""

        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()
        D_sap = []
        for s in range(density.nspins):
            D_sap.append([D_sp[s] for D_sp in D_asp])
        for nt_G, D_ap, mixer in zip(nt_sG, D_sap, self.mixers):
            mixer.mix(nt_G, D_ap)

    def reset(self):
        for mixer in self.mixers:
            mixer.reset()

    def get_charge_sloshing(self):
        """Return number of electrons moving around.

        Calculated as the integral of the absolute value of the change
        of the density from input to output."""

        if self.mixers[0].dNt is None:
            return None
        return sum([mixer.dNt for mixer in self.mixers])

    def set_charge_sloshing(self, dNt):
        for mixer in self.mixers:
            mixer.set_charge_sloshing(dNt / len(self.mixers))


class MixerSum_Broydn(BaseMixer_Broydn):
    def mix(self, density):
        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        BaseMixer_Broydn.mix(self, nt_G, D_asp)

        # Only new magnetization for spin density
        dnt_G = nt_sG[0] - nt_sG[1]
        #dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]

        # Construct new spin up/down densities 
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)
