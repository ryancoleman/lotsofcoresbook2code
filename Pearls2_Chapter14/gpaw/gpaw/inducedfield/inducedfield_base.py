import numpy as np

from ase.parallel import parprint

from ase.units import Bohr, Hartree
from gpaw.tddft import eV_to_aufrequency, aufrequency_to_eV  # TODO: remove

from gpaw.io.tar import Reader, Writer
from gpaw.poisson import PoissonSolver
from gpaw.fd_operators import Gradient

from gpaw.inducedfield.extend_grid import extend_grid, extend_array, deextend_array, move_atoms


def sendreceive_dict(comm, a_i, dest, b_i, src_i, iitems):
    """
    Send iitems of dictionary b_i distributed according to src_i
    to dictionary a_i in dest.
    """
    requests = []
    for i in iitems:
        # Send and receive
        if comm.rank == dest:
            # Check if dest has it already
            if dest == src_i[i]:
                a_i[i] = b_i[i].copy()
            else:
                requests.append(comm.receive(a_i[i], src_i[i], tag=i,
                                             block=False))
        elif comm.rank == src_i[i]:
            requests.append(comm.send(b_i[i], dest, tag=i, block=False))
    comm.waitall(requests)


class BaseInducedField(object):
    """Virtual base class for induced field calculations.
    
    Attributes:
    -----------
    omega_w: ndarray
        Frequencies for Fourier transform in atomic units
    folding: string
        Folding type ('Gauss' or 'Lorentz' or None)
    width: float
        Width parameter for folding
    Frho_wg: ndarray (complex)
        Fourier transform of induced charge density
    Fphi_wg: ndarray (complex)
        Fourier transform of induced electric potential
    Fef_wvg: ndarray (complex)
        Fourier transform of induced electric field
    Ffe_wg: ndarray (float)
        Fourier transform of field enhancement
    Fbgef_v: ndarray (float)
        Fourier transform of background electric field
    """
    
    def __init__(self, filename=None, paw=None, ws='all',
                  frequencies=None, folding='Gauss', width=0.08
                  ):
        """
        Parameters:
        -----------
        filename: string
            Filename of a previous ``InducedField`` to be loaded.
            Setting filename disables parameters ``frequencies``,
            ``folding`` and ``width``.
        paw: PAW object
            PAW object for InducedField
        ws: ndarray or list of ints
            Indices for frequencies that are read from given file.
            This parameter is neglected if ``filename`` is not given.
        frequencies: ndarray or list of floats
            Frequencies in eV for Fourier transforms.
            This parameter is neglected if ``filename`` is given.
        folding: string
            Folding type: 'Gauss' or 'Lorentz' or None:
            Gaussian or Lorentzian folding or no folding.
            This parameter is neglected if ``filename`` is given.
        width: float
            Width in eV for the Gaussian (sigma) or Lorentzian (eta) folding
            This parameter is neglected if ``filename`` is given.
        """
        self.dtype = complex
        
        if filename is not None:
            assert paw is not None, 'Provide PAW object also'
            self.initialize(paw, allocate=False)
            self.read(filename, ws=ws)
            return
        
        self.folding = folding
        self.width = width * eV_to_aufrequency
        self.set_folding(folding, width * eV_to_aufrequency)
        
        self.omega_w = np.asarray(frequencies) * eV_to_aufrequency
        self.nw = len(self.omega_w)
        
        self.nv = 3  # dimensionality of the space
        
        # These are allocated when calculated
        self.fieldgd = None
        self.Frho_wg = None
        self.Fphi_wg = None
        self.Fef_wvg = None
        self.Ffe_wg = None
        self.Fbgef_v = None
        
        if paw is not None:
            self.initialize(paw, allocate=True)
    
    def initialize(self, paw, allocate=True):
        self.allocated = False
        self.paw = paw
        self.world = paw.wfs.world
        self.domain_comm = paw.wfs.gd.comm
        self.band_comm = paw.wfs.band_comm
        self.kpt_comm = paw.wfs.kd.comm
        self.rank_a = paw.wfs.rank_a
        self.nspins = paw.density.nspins
        self.setups = paw.wfs.setups
        self.density = paw.density
        self.atoms = paw.atoms
        self.na = len(self.atoms.get_atomic_numbers())
        self.gd = self.density.gd
        self.stencil = self.paw.input_parameters.stencils[1]
        
        if allocate:
            self.allocate()
    
    def allocate(self):
        self.allocated = True
    
    def deallocate(self):
        self.fieldgd = None
        self.Frho_wg = None
        self.Fphi_wg = None
        self.Fef_wvg = None
        self.Ffe_wg = None
        self.Fbgef_v = None
        
    def set_folding(self, folding, width):
        """
        width: float
            Width in atomic units
        """
        if width is None:
            folding = None

        self.folding = folding
        if self.folding is None:
            self.width = None
        else:
            self.width = width
        
    def get_induced_density(self, from_density, gridrefinement):
        raise RuntimeError('Virtual member function called')

    def calculate_induced_field(self, from_density='comp',
                                   gridrefinement=2,
                                   extend_N_cd=None,
                                   deextend=False,
                                   poisson_nn=3, poisson_relax='J',
                                   gradient_n=3):

        Frho_wg, gd = self.get_induced_density(from_density,
                                               gridrefinement)
        
        # Always extend a bit to get field without jumps
        if extend_N_cd is None:
            extend_N_cd = gradient_n * np.ones(shape=(3, 2), dtype=np.int)
            deextend = True
        
        # Extend grid
        oldgd = gd
        egd, cell_cv, move_c = extend_grid(gd, extend_N_cd)
        Frho_we = egd.zeros((self.nw,), dtype=self.dtype)
        for w in range(self.nw):
            extend_array(Frho_wg[w], gd, Frho_we[w], egd)
        Frho_wg = Frho_we
        gd = egd
        if not deextend:
            # TODO: this will make atoms unusable with original grid
            self.atoms.set_cell(cell_cv, scale_atoms=False)
            move_atoms(self.atoms, move_c)
        
        # Allocate arrays
        Fphi_wg = gd.zeros((self.nw,), dtype=self.dtype)
        Fef_wvg = gd.zeros((self.nw, self.nv,), dtype=self.dtype)
        Ffe_wg = gd.zeros((self.nw,), dtype=float)
        
        for w in range(self.nw):
            # TODO: better output of progress
            #parprint('%d' % w)
            calculate_field(gd, Frho_wg[w], self.Fbgef_v,
                            Fphi_wg[w], Fef_wvg[w], Ffe_wg[w],
                            nv=self.nv,
                            poisson_nn=poisson_nn,
                            poisson_relax=poisson_relax,
                            gradient_n=gradient_n)

        # De-extend grid
        if deextend:
            Frho_wo = oldgd.zeros((self.nw,), dtype=self.dtype)
            Fphi_wo = oldgd.zeros((self.nw,), dtype=self.dtype)
            Fef_wvo = oldgd.zeros((self.nw, self.nv,), dtype=self.dtype)
            Ffe_wo = oldgd.zeros((self.nw,), dtype=float)
            for w in range(self.nw):
                deextend_array(Frho_wo[w], oldgd, Frho_wg[w], gd)
                deextend_array(Fphi_wo[w], oldgd, Fphi_wg[w], gd)
                deextend_array(Ffe_wo[w], oldgd, Ffe_wg[w], gd)
                for v in range(self.nv):
                    deextend_array(Fef_wvo[w][v], oldgd, Fef_wvg[w][v], gd)
            Frho_wg = Frho_wo
            Fphi_wg = Fphi_wo
            Fef_wvg = Fef_wvo
            Ffe_wg = Ffe_wo
            gd = oldgd
    
        # Store results
        self.field_from_density = from_density
        self.fieldgd = gd
        self.Frho_wg = Frho_wg
        self.Fphi_wg = Fphi_wg
        self.Fef_wvg = Fef_wvg
        self.Ffe_wg = Ffe_wg
        
    def _parse_readwritemode(self, mode):
        if type(mode) == str:
            try:
                readwrites = self.readwritemode_str_to_list[mode]
            except KeyError:
                raise IOError('unknown readwrite mode string')
        elif type(mode) == list:
            readwrites = mode
        else:
            raise IOError('unknown readwrite mode type')
        
        if any(k in readwrites for k in ['Frho', 'Fphi', 'Fef', 'Ffe']):
            readwrites.append('field')

        return readwrites

    def read(self, filename, mode='', ws='all', idiotproof=True):
        if idiotproof and not filename.endswith('.ind'):
            raise IOError('Filename must end with `.ind`.')
        
        reads = self._parse_readwritemode(mode)
        
        # Open reader (handles masters)
        tar = Reader(filename)

        # Actual read
        self.nw = tar.dimension('nw')
        if ws == 'all':
            ws = range(self.nw)
        self.nw = len(ws)
        self._read(tar, reads, ws)

        # Close
        tar.close()
        self.world.barrier()

    def _read(self, tar, reads, ws):
        # Test data type
        dtype = tar['datatype']
        if dtype != self.dtype:
            raise IOError('Data is an incompatible type.')

        # Read dimensions
        na = tar['na']
        self.nv = tar.dimension('nv')
        nspins = tar.dimension('nspins')
        ng = (tar.dimension('ng0'), tar.dimension('ng1'), tar.dimension('ng2'))
        
        # Test dimensions
        if na != self.na:
            raise IOError('natoms is incompatible with calculator')
        if nspins != self.nspins:
            raise IOError('nspins is incompatible with calculator')
        if (ng != self.gd.get_size_of_global_array()).any():
            raise IOError('grid is incompatible with calculator')

        # Folding
        folding = tar['folding']
        width = tar['width']
        self.set_folding(folding, width)
        
        # Frequencies
        self.omega_w = tar.get('omega_w')[ws]

        # Read field
        if 'field' in reads:
            nfieldg = (tar.dimension('nfieldg0'),
                       tar.dimension('nfieldg1'),
                       tar.dimension('nfieldg2'))
            self.field_from_density = tar['field_from_density']
            self.fieldgd = self.gd.new_descriptor(N_c=np.array(nfieldg) +\
                                                      np.array([1, 1, 1]))
    
        if 'Frho' in reads:
            self.Frho_wg = self.fieldgd.empty((self.nw), dtype=self.dtype)
            for w, wread in enumerate(ws):
                big_g = tar.get('Frho_wg', wread)
                self.fieldgd.distribute(big_g, self.Frho_wg[w])
        
        if 'Fphi' in reads:
            self.Fphi_wg = self.fieldgd.empty((self.nw), dtype=self.dtype)
            for w, wread in enumerate(ws):
                big_g = tar.get('Fphi_wg', wread)
                self.fieldgd.distribute(big_g, self.Fphi_wg[w])

        if 'Fef' in reads:
            self.Fef_wvg = self.fieldgd.empty((self.nw, self.nv),
                                              dtype=self.dtype)
            for w, wread in enumerate(ws):
                for v in range(self.nv):
                    big_g = tar.get('Fef_wvg', wread, v)
                    self.fieldgd.distribute(big_g, self.Fef_wvg[w][v])

        if 'Ffe' in reads:
            self.Ffe_wg = self.fieldgd.empty((self.nw), dtype=float)
            for w, wread in enumerate(ws):
                big_g = tar.get('Ffe_wg', wread)
                self.fieldgd.distribute(big_g, self.Ffe_wg[w])

    def write(self, filename, mode='', ws='all', idiotproof=True):
        """
        Parameters
        ----------
        mode: string or list of strings
            
        
        """
        if idiotproof and not filename.endswith('.ind'):
            raise IOError('Filename must end with `.ind`.')

        # Masters
        domainmaster = self.domain_comm.rank == 0
        bandmaster = self.band_comm.rank == 0
        kptmaster = self.kpt_comm.rank == 0
        master = domainmaster and bandmaster and kptmaster
#        master = self.world.rank == 0

        if master:
            if self.world.rank != 0:
                raise IOError('master not world master')

        writes = self._parse_readwritemode(mode)

        if ws == 'all':
            ws = range(self.nw)

        if 'field' in writes and self.fieldgd is None:
            raise IOError('field variables cannot be written ' +
                          'before they are calculated')

        # Open writer on master
        if master:
            tar = Writer(filename)
        else:
            tar = None
       
        # Actual write
        self._write(tar, writes, ws,
                     (master, domainmaster, bandmaster, kptmaster))
       
        # Close to flush changes
        if master:
            tar.close()

        # Make sure slaves don't return before master is done
        self.world.barrier()

    def _write(self, tar, writes, ws, masters):
        master, domainmaster, bandmaster, kptmaster = masters

        nw = len(ws)

        if master:
            # Write parameters/dimensions
            tar['datatype'] = {float: 'float', complex: 'complex'}[self.dtype]
            tar['folding'] = self.folding
            tar['width'] = self.width
            
            tar['na'] = self.na
            tar.dimension('nv', self.nv)
            tar.dimension('nw', nw)
            tar.dimension('nspins', self.nspins)
            
            # Write grid
            ng = self.gd.get_size_of_global_array()
            tar.dimension('ng0', ng[0])
            tar.dimension('ng1', ng[1])
            tar.dimension('ng2', ng[2])
            
            # Write field grid
            if 'field' in writes:
                nfieldg = self.fieldgd.get_size_of_global_array()
                tar.dimension('nfieldg0', nfieldg[0])
                tar.dimension('nfieldg1', nfieldg[1])
                tar.dimension('nfieldg2', nfieldg[2])
                tar['field_from_density'] = self.field_from_density

            # Write frequencies
            tar.add('omega_w', ('nw',), self.omega_w[ws], dtype=float)

            # Write background electric field
            tar.add('Fbgef_v', ('nv',), self.Fbgef_v, dtype=float)
       
            # Write system description
            # TODO: remove this and use ASE's atoms object instead
            if 'atoms' in writes:
                atomnum_a = self.atoms.get_atomic_numbers()
                atompos_a = self.atoms.get_positions()
                atomcell_cv = self.atoms.get_cell()
                
                tar.dimension('na', self.na)
                tar.add('atomnum_a', ('na',), atomnum_a, dtype=int)
                tar.add('atompos_a', ('na', 'nv',), atompos_a, dtype=float)
                tar.add('atomcell_cv', ('nv', 'nv',), atomcell_cv, dtype=float)
                
        if 'Frho' in writes:
            if master:
                tar.add('Frho_wg',
                        ('nw', 'nfieldg0', 'nfieldg1', 'nfieldg2', ),
                        dtype=self.dtype)
            for w in ws:
                big_g = self.fieldgd.collect(self.Frho_wg[w])
                if master:
                    tar.fill(big_g)

        if 'Fphi' in writes:
            if master:
                tar.add('Fphi_wg',
                        ('nw', 'nfieldg0', 'nfieldg1', 'nfieldg2', ),
                        dtype=self.dtype)
            for w in ws:
                big_g = self.fieldgd.collect(self.Fphi_wg[w])
                if master:
                    tar.fill(big_g)

        if 'Fef' in writes:
            if master:
                tar.add('Fef_wvg',
                        ('nw', 'nv', 'nfieldg0', 'nfieldg1', 'nfieldg2', ),
                        dtype=self.dtype)
            for w in ws:
                for v in range(self.nv):
                    big_g = self.fieldgd.collect(self.Fef_wvg[w][v])
                    if master:
                        tar.fill(big_g)
                        
        if 'Ffe' in writes:
            if master:
                tar.add('Ffe_wg',
                        ('nw', 'nfieldg0', 'nfieldg1', 'nfieldg2', ),
                        dtype=float)
            for w in ws:
                big_g = self.fieldgd.collect(self.Ffe_wg[w])
                if master:
                    tar.fill(big_g)


def calculate_field(gd, rho_g, bgef_v,
                      phi_g, ef_vg, fe_g,  # preallocated numpy arrays
                      nv=3, poisson_nn=3, poisson_relax='J', gradient_n=3):
    
    dtype = rho_g.dtype
    yes_complex = dtype == complex
    
    phi_g[:] = 0.0
    ef_vg[:] = 0.0
    fe_g[:] = 0.0
    tmp_g = gd.zeros(dtype=float)
    
    # Poissonsolver
    poissonsolver = PoissonSolver(nn=poisson_nn, relax=poisson_relax)
    poissonsolver.set_grid_descriptor(gd)
    poissonsolver.initialize()
    
    # Potential, real part
    poissonsolver.solve(tmp_g, rho_g.real.copy())
    phi_g += tmp_g
    # Potential, imag part
    if yes_complex:
        tmp_g[:] = 0.0
        poissonsolver.solve(tmp_g, rho_g.imag.copy())
        phi_g += 1.0j * tmp_g
    
    # Gradient
    gradient = [Gradient(gd, v, scale=1.0, n=gradient_n) for v in range(nv)]
    for v in range(nv):
        # Electric field, real part
        gradient[v].apply(-phi_g.real, tmp_g)
        ef_vg[v] += tmp_g
        # Electric field, imag part
        if yes_complex:
            gradient[v].apply(-phi_g.imag, tmp_g)
            ef_vg[v] += 1.0j * tmp_g
    
    # Electric field enhancement
    tmp_g[:] = 0.0  # total electric field norm
    bgefnorm = 0.0  # background electric field norm
    for v in range(nv):
        tmp_g += np.absolute(bgef_v[v] + ef_vg[v])**2
        bgefnorm += np.absolute(bgef_v[v])**2
        
    tmp_g = np.sqrt(tmp_g)
    bgefnorm = np.sqrt(bgefnorm)

    fe_g[:] = tmp_g / bgefnorm


def zero_pad(a):
    """Add zeros to both sides of all axis on array a.
    
    Not parallel safe.
    """
    
    z_shape = np.zeros(shape=len(a.shape))
    z_shape[-3:] = 2
    b_shape = np.array(a.shape) + z_shape
    
    b = np.zeros(shape=b_shape, dtype=a.dtype)
    
    b[..., 1:-1, 1:-1, 1:-1] = a
    return b


def read_data(filename, keys=None, ws='all'):
    """
    Read data arrays for post processing.
    
    Not parallel safe. No GridDescriptor, only numpy arrays.
    
    Parameters
    ----------
    filename: string
        File to be read.
    keys: list of strings
        Keys to be read.
    ws: list of ints
        Indices of frequencies to be read.
    """
    
    key_to_tarname = {'n0t_sG': 'n0t_sG',
                      'Fnt_wsG': 'Fnt_wsG',
                      'Frho_wg': 'Frho_wg',
                      'Fphi_wg': 'Fphi_wg',
                      'Fef_wvg': 'Fef_wvg',
                      'Ffe_wg': 'Ffe_wg',
                      'eps0_G': 'eps0_G'
                      }
    
    print('Reading %s' % (filename))
    
    if keys is None:
        keys = key_to_tarname.keys()  # all keys
    
    tar = Reader(filename)

    omega_w = tar.get('omega_w')

    if ws == 'all':
        ws = range(len(omega_w))
    else:
        omega_w = omega_w[ws]
    
    freq_w = omega_w * aufrequency_to_eV
    
    try:
        nspins = tar.dimension('nspins')
    except  KeyError:
        nspins = None
    
    na = tar['na']
    try:
        atomnum_a = tar.get('atomnum_a')
        atompos_av = tar.get('atompos_a')
        atomcell_cv = tar.get('atomcell_cv')
        Fbgef_v = tar.get('Fbgef_v')
        
        atom_a = []
        for a in range(na):
            atom_a.append({'atom': atomnum_a[a], 'pos': atompos_av[a]})
    except KeyError:
        atom_a = None
        atomcell_cv = None
        Fbgef_v = None
        print('no atoms')
    
    data = dict()
    data['freq_w'] = freq_w
    data['nspins'] = nspins
    data['na'] = na
    data['atom_a'] = atom_a
    data['cell_cv'] = atomcell_cv
    data['Fbgef_v'] = Fbgef_v
   
    try:
        data['corner1_v'] = tar.get('corner1_v')
        data['corner2_v'] = tar.get('corner2_v')
    except:
        print('no corners')

    try:
        FD_awsp = {}
        D0_asp = {}
        for a in range(na):
            FD_awsp[a] = tar.get('FD_%dwsp' % a)
            D0_asp[a] = tar.get('D0_%dsp' % a)
        data['FD_awsp'] = FD_awsp
        data['D0_asp'] = D0_asp
    except KeyError:
        print('no FD_awsp')

    for key in keys:
        try:
            if '_w' in key:
                tmp = zero_pad(tar.get(key_to_tarname[key], ws[0]))
                data[key] = np.empty((len(ws),) + tmp.shape, tmp.dtype)
                data[key][0] = tmp
                for w, wread in enumerate(ws[1:], 1):
                    data[key][w] = zero_pad(tar.get(key_to_tarname[key],
                                                    wread))
            else:
                data[key] = zero_pad(tar.get(key_to_tarname[key]))
        except KeyError:
            print('no %s' % key)
            pass

    tar.close()
    
    return data


# TOOD: remove/edit this function
def calculate_oscstr(Fn_wg, omega_w, box, kick):
    omega_w = np.array(omega_w)
    omega_w *= eV_to_aufrequency  # to a.u.
    box = box.copy() / Bohr  # to a.u
    volume = box.prod()
    ng = np.array(Fn_wg[0].shape)
    dV = volume / (ng - 1).prod()  # Note: data is zeropadded to both sides
    
    Fn_wg_im = Fn_wg.imag
    if kick[0] != 0.0:
        osc_w = -2 * omega_w / np.pi / kick[0] * ((Fn_wg_im.sum(axis=2)).sum(axis=2) * np.linspace(0, box[0], ng[0])).sum(axis=1) * dV
    elif kick[1] != 0.0:
        osc_w = -2 * omega_w / np.pi / kick[1] * ((Fn_wg_im.sum(axis=1)).sum(axis=2) * np.linspace(0, box[1], ng[1])).sum(axis=1) * dV
    elif kick[2] != 0.0:
        osc_w = -2 * omega_w / np.pi / kick[2] * ((Fn_wg_im.sum(axis=1)).sum(axis=1) * np.linspace(0, box[2], ng[2])).sum(axis=1) * dV
    return osc_w / Hartree  # to 1/eV


# TOOD: remove/edit this function
def calculate_polarizability(Fn_wg, box, kick):
    box = box.copy() / Bohr  # to a.u
    volume = box.prod()
    ng = np.array(Fn_wg[0].shape)
    dV = volume / (ng - 1).prod()  # Note: data is zeropadded to both sides
    
    if kick[0] != 0.0:
        pol_w = -1.0 / kick[0] * ((Fn_wg.sum(axis=2)).sum(axis=2) * np.linspace(0, box[0], ng[0])).sum(axis=1) * dV
    elif kick[1] != 0.0:
        pol_w = -1.0 / kick[1] * ((Fn_wg.sum(axis=1)).sum(axis=2) * np.linspace(0, box[1], ng[1])).sum(axis=1) * dV
    elif kick[2] != 0.0:
        pol_w = -1.0 / kick[2] * ((Fn_wg.sum(axis=1)).sum(axis=1) * np.linspace(0, box[2], ng[2])).sum(axis=1) * dV
    return pol_w / Hartree**2  # to 1/eV**2
