
import numpy as np

from gpaw import debug
from gpaw.io.tar import Reader, Writer
from gpaw.utilities import is_contiguous
from gpaw.analyse.observers import Observer
from gpaw.transformers import Transformer
from gpaw.tddft import attosec_to_autime, eV_to_aufrequency

# -------------------------------------------------------------------

class DensityFourierTransform(Observer):
    def __init__(self, timestep, frequencies, width=None, interval=1):
        """
        Parameters
        ----------
        timestep: float
            Time step in attoseconds (10^-18 s), e.g., 4.0 or 8.0
        frequencies: NumPy array or list of floats
            Frequencies in eV for Fourier transforms
        width: float or None
            Width of Gaussian envelope in eV, otherwise no envelope
        interval: int
            Number of timesteps between calls (used when attaching)
        """

        Observer.__init__(self, interval)
        self.timestep = interval * timestep * attosec_to_autime # autime
        self.omega_w = np.asarray(frequencies) * eV_to_aufrequency # autime^(-1)

        if width is None:
            self.sigma = None
        else:
            self.sigma = width * eV_to_aufrequency # autime^(-1)

        self.nw = len(self.omega_w)
        self.dtype = complex # np.complex128 really, but hey...
        self.Fnt_wsG = None
        self.Fnt_wsg = None

        self.Ant_sG = None
        self.Ant_sg = None

    def initialize(self, paw, allocate=True):
        self.allocated = False

        assert hasattr(paw, 'time') and hasattr(paw, 'niter'), 'Use TDDFT!'
        self.time = paw.time
        self.niter = paw.niter

        self.world = paw.wfs.world
        self.gd = paw.density.gd
        self.finegd = paw.density.finegd
        self.nspins = paw.density.nspins
        self.stencil = paw.input_parameters.stencils[1] # i.e. tar['InterpolationStencil']
        self.interpolator = paw.density.interpolator
        self.cinterpolator = Transformer(self.gd, self.finegd, self.stencil, \
                                        dtype=self.dtype)
        self.phase_cd = np.ones((3, 2), dtype=complex)

        self.Ant_sG = paw.density.nt_sG.copy() # TODO in allocate instead?

        # Attach to PAW-type object
        paw.attach(self, self.interval, density=paw.density)

        if allocate:
            self.allocate()

    def allocate(self):
        if not self.allocated:
            self.Fnt_wsG = self.gd.zeros((self.nw, self.nspins), \
                                        dtype=self.dtype)
            self.Fnt_wsg = None
            #self.Ant_sG = ...
            self.Ant_sg = None
            self.gamma_w = np.ones(self.nw, dtype=complex) * self.timestep
            self.allocated = True

        if debug:
            assert is_contiguous(self.Fnt_wsG, self.dtype)

    def interpolate_fourier_transform(self):
        if self.Fnt_wsg is None:
            self.Fnt_wsg = self.finegd.empty((self.nw, self.nspins), \
                                            dtype=self.dtype)

        if self.dtype == float:
            intapply = self.interpolator.apply
        else:
            intapply = lambda Fnt_G, Fnt_g: self.cinterpolator.apply(Fnt_G, \
                Fnt_g, self.phase_cd)

        for w in range(self.nw):
            for s in range(self.nspins):
                intapply(self.Fnt_wsG[w,s], self.Fnt_wsg[w,s])

    def interpolate_average(self):
        if self.Ant_sg is None:
            self.Ant_sg = self.finegd.empty(self.nspins, dtype=float)

        for s in range(self.nspins):
            self.interpolator.apply(self.Ant_sG[s], self.Ant_sg[s])
            
    def update(self, density):

        # Update time
        # t[N] = t[N-1] + dt[N-1] #TODO better time-convention?
        self.time += self.timestep

        # Complex exponential with/without finite-width envelope
        f_w = np.exp(1.0j*self.omega_w*self.time)
        if self.sigma is not None:
            f_w *= np.exp(-self.time**2*self.sigma**2/2.0)

        # Update Fourier transformed density components
        # Fnt_wG[N] = Fnt_wG[N-1] + 1/sqrt(pi) * (nt_G[N]-avg_nt_G[N-1]) \
        #     * (f[N]*t[N] - gamma[N-1]) * dt[N]/(t[N]+dt[N])
        for w in range(self.nw):
            self.Fnt_wsG[w] += 1/np.pi**0.5 * (density.nt_sG - self.Ant_sG) \
                * (f_w[w]*self.time - self.gamma_w[w]) * self.timestep \
                / (self.time + self.timestep)

        # Update the cumulative phase factors
        # gamma[N] = gamma[N-1] + f[N]*dt[N]
        self.gamma_w += f_w * self.timestep

        # If dt[N] = dt for all N and sigma = 0, then this simplifies to:
        # gamma[N] = Sum_{n=0}^N exp(i*omega*n*dt) * dt
        # = (1 - exp(i*omega*(N+1)*dt)) / (1 - exp(i*omega*dt)) * dt

        # Update average density
        # Ant_G[N] = (t[N]*Ant_G[N-1] + nt_G[N]*dt[N])/(t[N]+dt[N])
        self.Ant_sG = (self.time*self.Ant_sG + density.nt_sG*self.timestep) \
            / (self.time + self.timestep)

    def get_fourier_transform(self, frequency=0, spin=0, gridrefinement=1):
        if gridrefinement == 1:
            return self.Fnt_wsG[frequency, spin]
        elif gridrefinement == 2:
            if self.Fnt_wsg is None:
                self.interpolate_fourier_transform()
            return self.Fnt_wsg[frequency, spin]
        else:
            raise NotImplementedError('Arbitrary refinement not implemented')

    def get_average(self, spin=0, gridrefinement=1):
        if gridrefinement == 1:
            return self.Ant_sG[spin]
        elif gridrefinement == 2:
            if self.Ant_sg is None:
                self.interpolate_average()
            return self.Ant_sg[spin]
        else:
            raise NotImplementedError('Arbitrary refinement not implemented')

    def read(self, filename, idiotproof=True):
        if idiotproof and not filename.endswith('.ftd'):
            raise IOError('Filename must end with `.ftd`.')

        tar = Reader(filename)

        # Test data type
        dtype = {'Float':float, 'Complex':complex}[tar['DataType']]
        if dtype != self.dtype:
            raise IOError('Data is an incompatible type.')

        # Test time
        time = tar['Time']
        if idiotproof and abs(time-self.time) >= 1e-9:
            raise IOError('Timestamp is incompatible with calculator.')

        # Test timestep (non-critical)
        timestep = tar['TimeStep']
        if abs(timestep - self.timestep) > 1e-12:
            print('Warning: Time-step has been altered. (%lf -> %lf)' \
                % (self.timestep, timestep))
        self.timestep = timestep

        # Test dimensions
        nw = tar.dimension('nw')
        nspins = tar.dimension('nspins')
        ng = (tar.dimension('ngptsx'), tar.dimension('ngptsy'), \
              tar.dimension('ngptsz'),)

        if (nw != self.nw or nspins != self.nspins or
            (ng != self.gd.get_size_of_global_array()).any()):
            raise IOError('Data has incompatible shapes.')

        # Test width (non-critical)
        sigma = tar['Width']
        if ((sigma is None)!=(self.sigma is None) or # float <-> None
            (sigma is not None and self.sigma is not None and \
             abs(sigma - self.sigma) > 1e-12)): # float -> float
            print('Warning: Width has been altered. (%s -> %s)' \
                % (self.sigma, sigma))
        self.sigma = sigma

        # Read frequencies
        self.omega_w[:] = tar.get('Frequency')

        # Read cumulative phase factors
        self.gamma_w[:] = tar.get('PhaseFactor')

        # Read average densities on master and distribute
        for s in range(self.nspins):
            all_Ant_G = tar.get('Average', s)
            self.gd.distribute(all_Ant_G, self.Ant_sG[s])

        # Read fourier transforms on master and distribute
        for w in range(self.nw):
            for s in range(self.nspins):
                all_Fnt_G = tar.get('FourierTransform', w, s)
                self.gd.distribute(all_Fnt_G, self.Fnt_wsG[w,s])

        # Close for good measure
        tar.close()

    def write(self, filename, idiotproof=True):
        if idiotproof and not filename.endswith('.ftd'):
            raise IOError('Filename must end with `.ftd`.')

        master = self.world.rank == 0

        # Open writer on master and set parameters/dimensions
        if master:
            tar = Writer(filename)
            tar['DataType'] = {float:'Float', complex:'Complex'}[self.dtype]
            tar['Time'] = self.time
            tar['TimeStep'] = self.timestep #non-essential
            tar['Width'] = self.sigma

            tar.dimension('nw', self.nw)
            tar.dimension('nspins', self.nspins)

            # Create dimensions for varioius netCDF variables:
            ng = self.gd.get_size_of_global_array()
            tar.dimension('ngptsx', ng[0])
            tar.dimension('ngptsy', ng[1])
            tar.dimension('ngptsz', ng[2])

            # Write frequencies
            tar.add('Frequency', ('nw',), self.omega_w, dtype=float)

            # Write cumulative phase factors
            tar.add('PhaseFactor', ('nw',), self.gamma_w, dtype=self.dtype)

        # Collect average densities on master and write
        if master:
            tar.add('Average', ('nspins', 'ngptsx', 'ngptsy', 
                'ngptsz', ), dtype=float)
        for s in range(self.nspins):
            big_Ant_G = self.gd.collect(self.Ant_sG[s])
            if master:
                tar.fill(big_Ant_G)

        # Collect fourier transforms on master and write
        if master:
            tar.add('FourierTransform', ('nw', 'nspins', 'ngptsx', 'ngptsy', \
                'ngptsz', ), dtype=self.dtype)
        for w in range(self.nw):
            for s in range(self.nspins):
                big_Fnt_G = self.gd.collect(self.Fnt_wsG[w,s])
                if master:
                    tar.fill(big_Fnt_G)

        # Close to flush changes
        if master:
            tar.close()

        # Make sure slaves don't return before master is done
        self.world.barrier()

    def dump(self, filename):
        if debug:
            assert is_contiguous(self.Fnt_wsG, self.dtype)
            assert is_contiguous(self.Ant_sG, float)

        all_Fnt_wsG = self.gd.collect(self.Fnt_wsG)
        all_Ant_sG = self.gd.collect(self.Ant_sG)

        if self.world.rank == 0:
            all_Fnt_wsG.dump(filename)
            all_Ant_sG.dump(filename+'_avg') # crude but easy
            self.omega_w.dump(filename+'_omega') # crude but easy
            self.gamma_w.dump(filename+'_gamma') # crude but easy

    def load(self, filename):
        if self.world.rank == 0:
            all_Fnt_wsG = np.load(filename)
            all_Ant_sG = np.load(filename+'_avg') # crude but easy
        else:
            all_Fnt_wsG = None
            all_Ant_sG = None

        if debug:
            assert all_Fnt_wsG is None or is_contiguous(all_Fnt_wsG, self.dtype)
            assert all_Ant_sG is None or is_contiguous(all_Ant_sG, float)

        if not self.allocated:
            self.allocate()

        self.gd.distribute(all_Fnt_wsG, self.Fnt_wsG)
        self.gd.distribute(all_Ant_sG, self.Ant_sG)

        self.omega_w = np.load(filename+'_omega') # crude but easy
        self.gamma_w = np.load(filename+'_gamma') # crude but easy

