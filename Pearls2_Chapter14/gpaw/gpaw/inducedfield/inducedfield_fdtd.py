import numpy as np

from ase.parallel import parprint
from ase.units import Bohr

from gpaw import debug
from gpaw.analyse.observers import Observer
from gpaw.transformers import Transformer
from gpaw.lfc import BasisFunctions
from gpaw.utilities import unpack2, is_contiguous

from gpaw.inducedfield.inducedfield_base import BaseInducedField, sendreceive_dict


class FDTDInducedField(BaseInducedField, Observer):
    """Induced field class for FDTD.
    
    Attributes (see also ``BaseInducedField``):
    -------------------------------------------
    time: float
        Current time
    Fn_wsG: ndarray (complex)
        Fourier transform of induced polarization charge density
    n0t_sG: ndarray (float)
        Ground state charge density
    FD_awsp: dict of ndarray (complex)
        Fourier transform of induced D_asp
    D0_asp: dict of ndarray (float)
        Ground state D_asp
    """
    
    def __init__(self, filename=None, paw=None, ws='all',
                  frequencies=None, folding='Gauss', width=0.08,
                  interval=1, restart_file=None
                  ):
        """
        Parameters (see also ``BaseInducedField``):
        -------------------------------------------
        paw: TDDFT object
            TDDFT object for time propagation
        width: float
            Width in eV for the Gaussian (sigma) or Lorentzian (eta) folding
            Gaussian   = exp(- (1/2) * sigma^2 * t^2)
            Lorentzian = exp(- eta * t)
        interval: int
            Number of timesteps between calls (used when attaching)
        restart_file: string
            Name of the restart file
        """

        Observer.__init__(self, interval)
        # From observer:
        # self.niter
        # self.interval
        
        # Restart file
        self.restart_file = restart_file
        
        # These are allocated in allocate()
        self.Fn_wsG = None
        self.n0_G = None


        self.readwritemode_str_to_list = \
        {'': ['Fn', 'n0', 'FD', 'atoms'],
         'all': ['Fn', 'n0', 'FD',
                 'Frho', 'Fphi', 'Fef', 'Ffe', 'eps0', 'atoms'],
         'field': ['Frho', 'Fphi', 'Fef', 'Ffe', 'eps0', 'atoms']}

        BaseInducedField.__init__(self, filename, paw, ws,
                                  frequencies, folding, width)
        
    def initialize(self, paw, allocate=True):
        BaseInducedField.initialize(self, paw, allocate)
        
        # FDTD replacements and overwrites
        self.fdtd   = paw.hamiltonian.poisson
        #self.gd     = self.fdtd.cl.gd.refine()
        self.gd     = self.fdtd.cl.gd
               
        assert hasattr(paw, 'time') and hasattr(paw, 'niter'), 'Use TDDFT!'
        self.time = paw.time
        self.niter = paw.niter
        
        # TODO: remove this requirement
        assert np.count_nonzero(paw.kick_strength) > 0, \
        'Apply absorption kick before %s' % self.__class__.__name__
        
        # Background electric field
        self.Fbgef_v = paw.kick_strength

        # Attach to PAW-type object
        paw.attach(self, self.interval)
        # TODO: write more details (folding, freqs, etc)
        parprint('%s: Attached ' % self.__class__.__name__)

    def set_folding(self, folding, width):
        BaseInducedField.set_folding(self, folding, width)
        
        if self.folding is None:
            self.envelope = lambda t: 1.0
        else:
            if self.folding == 'Gauss':
                self.envelope = lambda t: np.exp(- 0.5 * self.width**2 * t**2)
            elif self.folding == 'Lorentz':
                self.envelope = lambda t: np.exp(- self.width * t)
            else:
                raise RuntimeError('unknown folding "' + self.folding + '"')
        
    def allocate(self):
        if not self.allocated:
            
            # Ground state charge density
            self.n0_G = (-1.0) * self.paw.hamiltonian.poisson.classical_material.sign * \
                                 self.paw.hamiltonian.poisson.classical_material.charge_density.copy()
            
            # Fourier transformed charge density
            self.Fn_wsG = self.paw.hamiltonian.poisson.cl.gd.zeros((self.nw, self.nspins),
                                                                   dtype=self.dtype)
    
            self.allocated = True
            
        if debug:
            assert is_contiguous(self.Ft_wG, self.dtype)

    def deallocate(self):
        BaseInducedField.deallocate(self)
        self.n0t_G = None
        self.Fnt_wG = None
        
    def update(self):
        # Update time
        self.time = self.paw.time
        time_step = self.paw.time_step

        # Complex exponential with envelope
        f_w = np.exp(1.0j * self.omega_w * self.time) * \
              self.envelope(self.time) * time_step

        # Time-dependent quantities
        n_G = (-1.0) * self.fdtd.classical_material.charge_density * self.fdtd.classical_material.sign

        # Update Fourier transforms
        for w in range(self.nw):
            self.Fn_wsG[w] += (n_G - self.n0_G) * f_w[w]

        # Restart file
        if self.restart_file is not None and \
           self.niter % self.paw.dump_interval == 0:
            self.write(self.restart_file)
            parprint('%s: Wrote restart file %s' % (self.__class__.__name__, self.restart_file))
    
        
    def get_induced_density(self, from_density, gridrefinement):
        #Frho_wg = -self.Fn_wG.sum(axis=1).copy()
        #Fn_wG_global = -self.gd.collect(self.Fn_wG).sum(axis=1)
        #if self.gd.comm.rank==0:
        
        Frho_wg = []
        for w in range(self.nw):
            if self.gd == self.fdtd.cl.gd:
                Frho_wg.append(self.Fn_wsG[w].sum(axis=0))
            else:
                Frho_wg.append(Transformer(self.fdtd.cl.gd,
                                           self.gd,
                                           self.stencil,
                                           dtype=self.dtype).apply(self.Fn_wsG[w].sum(axis=0)))
        
        return Frho_wg, self.gd
    
    def _read(self, tar, reads, ws):
        BaseInducedField._read(self, tar, reads, ws)
        
        # Test time
        time = tar['time']
        if abs(time - self.time) >= 1e-9:
            raise IOError('Timestamp is incompatible with calculator.')

        # Allocate
        self.allocate()

        # Read arrays
        if 'n0' in reads:
            big_g = tar.get('n0_G')
            self.gd.distribute(big_g, self.n0_G)
        
        if 'Fn' in reads:
            for w, wread in enumerate(ws):
                big_g = tar.get('Fn_wsG', wread)
                self.gd.distribute(big_g, self.Fn_wsG[w])
        
        if 'eps0' in reads:
            self.eps0_G = self.gd.empty(dtype=float)
            big_g = tar.get('eps0_G')
            self.gd.distribute(big_g, self.eps0_G)
        else:
            self.eps0_G = None


    def _write(self, tar, writes, ws, masters):
        # Swap classical and quantum cells, and shift atom positions for the time of writing
        qmcell = self.atoms.get_cell()
        self.atoms.set_cell(self.fdtd.cl.cell*Bohr) # Set classical cell
        self.atoms.positions = self.atoms.get_positions() + self.fdtd.qm.corner1*Bohr
        BaseInducedField._write(self, tar, writes, ws, masters)
        self.atoms.set_cell(qmcell) # Restore quantum cell to the atoms object
        self.atoms.positions = self.atoms.get_positions() - self.fdtd.qm.corner1*Bohr

        master, domainmaster, bandmaster, kptmaster = masters
        
        # Write time propagation status
        if master:
            tar['time'] = self.time
        
        # Mask, interpolation approach:
        #self.eps0_G = self.fdtd.classical_material.permittivityValue(omega=0.0) - self.fdtd.classical_material.epsInfty
        #self.eps0_G= -interpolator.apply(self.eps0_G)
        
        # Mask, direct approach:
        self.eps0_G = self.fdtd.cl.gd.zeros()
        for component in self.fdtd.classical_material.components:
            self.eps0_G += 1.0 * component.get_mask(gd = self.fdtd.cl.gd, verbose = False)
        
        # Write time propagation arrays
        if 'n0' in writes:
            
            if master:
                tar.add('n0_G',
                        ('ng0', 'ng1', 'ng2'),
                        dtype=float)
            
            if self.fdtd.cl.gd == self.gd:
                big_g = self.gd.collect(self.n0_G)
            else:
                big_g = self.gd.collect(Transformer(self.fdtd.cl.gd, self.gd, self.stencil, float).apply(self.n0_G))
            
            if master:
                tar.fill(big_g)
                
        if 'Fn' in writes:
            if master:
                tar.add('Fn_wsG',
                        ('nw', 'nspins', 'ng0', 'ng1', 'ng2'),
                        dtype=self.dtype)
            for w in ws:
                #big_g = self.gd.collect(self.Fn_wG[w])
                if self.fdtd.cl.gd == self.gd:
                    big_g = self.gd.collect(self.Fn_wsG[w])
                else:
                    big_g = self.gd.collect(Transformer(self.fdtd.cl.gd, self.gd, self.stencil, self.dtype).apply(self.Fn_wsG[w]))
                #big_g = self.gd.collect(interpolator_float.apply(self.n0_G.copy()))
                if master:
                    tar.fill(big_g)
         
        if 'eps0' in writes:
            if master:
                tar.add('eps0_G',
                        ('ng0', 'ng1', 'ng2'),
                        dtype=float)
            #big_g = self.fieldgd.collect(self.eps0_G)
            if self.fdtd.cl.gd == self.gd:
                big_g = self.gd.collect(self.eps0_G)
            else:
                big_g = self.gd.collect(Transformer(self.fdtd.cl.gd, self.gd, self.stencil, float).apply(self.eps0_G))
            if master:
                tar.fill(big_g)
        
        if master and hasattr(self.fdtd, 'qm') and hasattr(self.fdtd.qm, 'corner1'):
            tar.add('corner1_v', ('nv',), self.fdtd.qm.corner1, dtype=float)
            tar.add('corner2_v', ('nv',), self.fdtd.qm.corner2, dtype=float)
        
        self.fdtd.cl.gd.comm.barrier()
            
            