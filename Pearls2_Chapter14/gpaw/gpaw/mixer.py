# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
See Kresse, Phys. Rev. B 54, 11169 (1996)
"""

import numpy as np
from numpy.fft import fftn, ifftn

from gpaw.utilities.blas import axpy
from gpaw.fd_operators import FDOperator
from gpaw.utilities.tools import construct_reciprocal


class BaseMixer:
    """Pulay density mixer."""
    
    def __init__(self, beta=0.1, nmaxold=3, weight=50.0, dotprod=None):
        """Construct density-mixer object.

        Parameters:

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

        self.dNt = None

        self.mix_rho = False
        
        if dotprod is not None:  # slightly ugly way to override
            self.dotprod = dotprod

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
                                     [(0, 0, 0),  # a
                                      (-1, 0, 0), (1, 0, 0),  # b
                                      (0, -1, 0), (0, 1, 0),
                                      (0, 0, -1), (0, 0, 1),
                                      (1, 1, 0), (1, 0, 1), (0, 1, 1),  # c
                                      (1, -1, 0), (1, 0, -1), (0, 1, -1),
                                      (-1, 1, 0), (-1, 0, 1), (0, -1, 1),
                                      (-1, -1, 0), (-1, 0, -1), (0, -1, -1),
                                      (1, 1, 1), (1, 1, -1), (1, -1, 1),  # d
                                      (-1, 1, 1), (1, -1, -1), (-1, -1, 1),
                                      (-1, 1, -1), (-1, -1, -1)],
                                     gd, float).apply
            self.mR_G = gd.empty()
        
    def initialize(self, density):
        self.initialize_metric(density.gd)

    def reset(self):
        """Reset Density-history.

        Called at initialization and after each move of the atoms.

        my_nuclei:   All nuclei in local domain.
        """
        
        # History for Pulay mixing of densities:
        self.nt_iG = []  # Pseudo-electron densities
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

    def calculate_charge_sloshing(self, R_G):
        return self.gd.integrate(np.fabs(R_G))

    def set_charge_sloshing(self, dNt):
        self.dNt = dNt
        
    def mix(self, nt_G, D_ap):
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

            # Calculate new residual (difference between input and
            # output density):
            R_G = nt_G - self.nt_iG[-1]
            self.dNt = self.calculate_charge_sloshing(R_G)
            self.R_iG.append(R_G)
            self.dD_iap.append([])
            for D_p, D_ip in zip(D_ap, self.D_iap[-1]):
                self.dD_iap[-1].append(D_p - D_ip)

            # Update matrix:
            A_ii = np.zeros((iold, iold))
            i2 = iold - 1
            
            if self.metric is None:
                mR_G = R_G
            else:
                mR_G = self.mR_G
                self.metric(R_G, mR_G)
                
            for i1, R_1G in enumerate(self.R_iG):
                a = self.gd.comm.sum(self.dotprod(R_1G, mR_G, self.dD_iap[i1],
                                                  self.dD_iap[-1]))
                A_ii[i1, i2] = a
                A_ii[i2, i1] = a
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
            #for D_p, D_ip, dD_ip in self.D_a:
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

    # may presently be overridden by passing argument in constructor
    def dotprod(self, R1_G, R2_G, dD1_ap, dD2_ap):
        return np.vdot(R1_G, R2_G).real

    def estimate_memory(self, mem, gd):
        gridbytes = gd.bytecount()
        mem.subnode('nt_iG, R_iG', 2 * self.nmaxold * gridbytes)

    def __repr__(self):
        classname = self.__class__.__name__
        template = '%s(beta=%f, nmaxold=%d, weight=%f)'
        string = template % (classname, self.beta, self.nmaxold, self.weight)
        return string


class ExperimentalDotProd:
    def __init__(self, calc):
        self.calc = calc
    
    def __call__(self, R1_G, R2_G, dD1_ap, dD2_ap):
        prod = np.vdot(R1_G, R2_G).real
        setups = self.calc.wfs.setups
        # okay, this is a bit nasty because it depends on dD1_ap
        # and its friend having come from D_asp.values() and the dictionaries
        # not having been modified.  This is probably true... for now.
        avalues = self.calc.density.D_asp.keys()
        for a, dD1_p, dD2_p in zip(avalues, dD1_ap, dD2_ap):
            I4_pp = setups[a].four_phi_integrals()
            dD4_pp = np.outer(dD1_p, dD2_p)  # not sure if corresponds quite
            prod += (I4_pp * dD4_pp).sum()
        return prod

        
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
            mixer = BaseMixer(self.beta, self.nmaxold, self.weight,
                              self.dotprod)
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

        Parameters:

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


class BroydenBaseMixer:
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
        #self.d_nt_G = []
        #self.d_D_ap = []
      
        self.R_iG = []
        self.dD_iap = []

        self.nt_iG = []
        self.D_iap = []
        self.c_G = []
        self.v_G = []
        self.u_G = []
        self.u_D = []
        self.dNt = None
        
    def get_charge_sloshing(self):
        return self.dNt
    
    def mix(self, nt_G, D_ap):
        if self.step > 2:
            del self.R_iG[0]
            for d_Dp in self.dD_iap:
                del d_Dp[0]
        if self.step > 0:
            self.R_iG.append(nt_G - self.nt_iG[-1])
            for d_Dp, D_p, D_ip in zip(self.dD_iap, D_ap, self.D_iap):
                d_Dp.append(D_p - D_ip[-1])
            fmin_G = self.gd.integrate(self.R_iG[-1] * self.R_iG[-1])
            self.dNt = self.gd.integrate(np.fabs(self.R_iG[-1]))
            if self.verbose:
                print('Mixer: broyden: fmin_G = %f fmin_D = %f' % fmin_G)
        if self.step == 0:
            self.eta_G = np.empty(nt_G.shape)
            self.eta_D = []
            for D_p in D_ap:
                self.eta_D.append(0)
                self.u_D.append([])
                self.D_iap.append([])
                self.dD_iap.append([])
        else:
            if self.step >= 2:
                del self.c_G[:]
                if len(self.v_G) >= self.nmaxold:
                    del self.u_G[0]
                    del self.v_G[0]
                    for u_D in self.u_D:
                        del u_D[0]
                temp_nt_G = self.R_iG[1] - self.R_iG[0]
                self.v_G.append(temp_nt_G / self.gd.integrate(temp_nt_G *
                                                              temp_nt_G))
                if len(self.v_G) < self.nmaxold:
                    nstep = self.step - 1
                else:
                    nstep = self.nmaxold
                for i in range(nstep):
                    self.c_G.append(self.gd.integrate(self.v_G[i] *
                                                      self.R_iG[1]))
                self.u_G.append(self.beta * temp_nt_G + self.nt_iG[1] -
                                self.nt_iG[0])
                for d_Dp, u_D, D_ip in zip(self.dD_iap, self.u_D, self.D_iap):
                    temp_D_ap = d_Dp[1] - d_Dp[0]
                    u_D.append(self.beta * temp_D_ap + D_ip[1] - D_ip[0])
                usize = len(self.u_G)
                for i in range(usize - 1):
                    a_G = self.gd.integrate(self.v_G[i] * temp_nt_G)
                    axpy(-a_G, self.u_G[i], self.u_G[usize - 1])
                    for u_D in self.u_D:
                        axpy(-a_G, u_D[i], u_D[usize - 1])
            self.eta_G = self.beta * self.R_iG[-1]
            for i, d_Dp in enumerate(self.dD_iap):
                self.eta_D[i] = self.beta * d_Dp[-1]
            usize = len(self.u_G)
            for i in range(usize):
                axpy(-self.c_G[i], self.u_G[i], self.eta_G)
                for eta_D, u_D in zip(self.eta_D, self.u_D):
                    axpy(-self.c_G[i], u_D[i], eta_D)
            axpy(-1.0, self.R_iG[-1], nt_G)
            axpy(1.0, self.eta_G, nt_G)
            for D_p, d_Dp, eta_D in zip(D_ap, self.dD_iap, self.eta_D):
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

        
class BroydenMixer(BroydenBaseMixer):
    """Mix spin up and down densities separately"""

    def initialize(self, density):
        self.mixers = []
        for s in range(density.nspins):
            mixer = BroydenBaseMixer()
            mixer.initialize(density)
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


class BroydenMixerSum(BroydenBaseMixer):
    def mix(self, density):
        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        BroydenBaseMixer.mix(self, nt_G, D_asp)

        # Only new magnetization for spin density
        dnt_G = nt_sG[0] - nt_sG[1]
        #dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]

        # Construct new spin up/down densities
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)

        
class FFTBaseMixer(BaseMixer):
    """Mix the density in Fourier space"""
    def __init__(self, beta=0.4, nmaxold=3, weight=20.0):
        BaseMixer.__init__(self, beta, nmaxold, weight)

    def initialize(self, density):
        if density.gd.comm.size > 1:
            err = 'FFT Mixer cannot be used with domain decomposition'
            raise NotImplementedError(err)
        self.initialize_metric(density.gd)

    def initialize_metric(self, gd):
        self.gd = gd
        k2_Q, N3 = construct_reciprocal(self.gd)

        self.metric = ReciprocalMetric(self.weight, k2_Q)
        self.mR_G = gd.empty(dtype=complex)

    def calculate_charge_sloshing(self, R_Q):
        return self.gd.integrate(np.fabs(ifftn(R_Q).real))

    def mix(self, nt_G, D_ap):
        # Transform real-space density to Fourier space
        nt_Q = np.ascontiguousarray(fftn(nt_G))

        BaseMixer.mix(self, nt_Q, D_ap)

        # Return density in real space
        nt_G[:] = np.ascontiguousarray(ifftn(nt_Q).real)

        
class FFTMixer(FFTBaseMixer):
    """Mix spin up and down densities separately"""

    def initialize(self, density):
        self.mixers = []
        for s in range(density.nspins):
            mixer = FFTBaseMixer(self.beta, self.nmaxold, self.weight)
            mixer.initialize(density)
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

            
class FFTMixerSum(FFTBaseMixer):
    """For pseudo electron densities, mix the total charge density and for
    density matrices, mix spin up and densities separately.
    Magnetization density is not mixed, i.e new magnetization density is used
    """

    def mix(self, density):
        nt_sG = density.nt_sG
        D_asp = density.D_asp.values()

        # Mix density
        nt_G = density.nt_sG.sum(0)
        FFTBaseMixer.mix(self, nt_G, D_asp)

        # Only new magnetization for spin density
        dnt_G = nt_sG[0] - nt_sG[1]
        #dD_ap = [D_sp[0] - D_sp[1] for D_sp in D_asp]

        # Construct new spin up/down densities
        nt_sG[0] = 0.5 * (nt_G + dnt_G)
        nt_sG[1] = 0.5 * (nt_G - dnt_G)

        
class FFTMixerDif(FFTBaseMixer):
    """Mix the charge density and magnetization density separately"""
    
    def __init__(self, beta=0.1, nmaxold=3, weight=20.0,
                 beta_m=0.7, nmaxold_m=2, weight_m=1.0):
        """Construct density-mixer object.

        Parameters:

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
        self.beta1 = beta
        self.nmaxold = nmaxold
        self.weight = weight

        self.beta_m = beta_m
        self.nmaxold_m = nmaxold_m
        self.weight_m = weight_m
        self.dNt = None

        self.mix_rho = False

    def initialize(self, density):
        assert density.nspins == 2
        self.mixer = FFTBaseMixer(self.beta, self.nmaxold, self.weight)
        self.mixer.initialize(density)
        self.mixer_m = FFTBaseMixer(self.beta_m, self.nmaxold_m, self.weight_m)
        self.mixer_m.initialize(density)

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

        
class ReciprocalMetric:
    def __init__(self, weight, k2_Q):
        self.k2_Q = k2_Q
        k2_min = np.min(self.k2_Q)
        self.q1 = (weight - 1) * k2_min

    def __call__(self, R_Q, mR_Q):
            mR_Q[:] = R_Q * (1.0 + self.q1 / self.k2_Q)
