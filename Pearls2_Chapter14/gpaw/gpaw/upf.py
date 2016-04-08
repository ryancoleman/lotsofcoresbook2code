"""Unified pseudopotential format (UPF).

UPF refers to a collection of different formats to store
pseudopotentials or PAW setups.  This file attempts to parse them and
provide pseudopotential objects for use with GPAW.
"""

from optparse import OptionParser
from xml.etree.ElementTree import parse as xmlparse, ParseError, fromstring

import numpy as np
from ase.data import atomic_numbers

from gpaw.basis_data import Basis, BasisFunction
from gpaw.pseudopotential import PseudoPotential, screen_potential, \
    figure_out_valence_states
from gpaw.spline import Spline
from gpaw.utilities import pack2, divrl


class UPFStateSpec:
    def __init__(self, index, label, l, values, occupation=None, n=None):
        self.index = index
        self.label = label
        self.l = l
        self.values = values
        self.occupation = occupation
        self.n = n


class UPFCompensationChargeSpec: # Not used right now.......
    def __init__(self, i, j, l, qint, values, coefs):
        self.i = i
        self.j = j
        self.l = l
        self.qint = qint
        self.values = values
        self.coefs = coefs

    def __repr__(self):
        return ('CompensationCharge(i=%d, j=%d, l=%d, qint=%s, '
                'ccvalues=[...], coefs=[...(%d)...]') % (self.i, self.j,
                                                         self.l, self.qint,
                                                         len(self.coefs))


def trim_outer_zeros(values, threshold=0.0):
    values = list(values)
    if abs(values[-1]) <= threshold:
        while abs(values[-2]) <= threshold:
            values.pop()
        assert abs(values[-1]) <= threshold
    # This will give an error if the array length becomes smaller than 2.
    # Which would be fine I guess.
    return np.array(values)


def parse_upf(fname):
    """Parse UPF pseudopotential from file descriptor.

    Return a dictionary with the parsed data.

    UPF is a custom format with some XML tags in it.  The files are
    not well-formed XML.  We will just do some custom parsing rather
    than involving an xml library.
    """

    pp = {}

    try:
        root = xmlparse(fname).getroot()
    except ParseError:
        # Typical!  Non-well-formed file full of probable FORTRAN output.
        # We'll try to insert our own header and see if things go well.
        root = fromstring('\n'.join(['<xml>', open(fname).read(), '</xml>']))

    pp['fname'] = fname

    header_element = root.find('PP_HEADER')
    header = header_element.attrib

    # v201 is the sensible version.
    # There are other files out there and we will have to deal with their
    # inconsistencies in the most horrendous way conceivable: by if-statements.
    v201 = (root.attrib.get('version') == '2.0.1')

    def toarray(element):
        attr = element.attrib
        numbers = map(float, element.text.split())
        if attr:
            assert attr['type'] == 'real'
            assert int(attr['size']) == len(numbers)
            icut = attr.get('cutoff_radius_index')
            # XXX There's also something called ultrasoft_cutoff_radius
            if icut is not None:
                numbers = numbers[:int(icut)]

        return np.array(numbers)

    if v201:
        for key, val in header.items():
            header[key] = val.strip() # some values have whitespace...
        for key in ['is_paw', 'is_coulomb', 'has_so', 'has_wfc', 'has_gipaw',
                    'paw_as_gipaw', 'core_correction']:
            header[key] = {'F': False, 'T': True}[header[key]]
        for key in ['z_valence', 'total_psenergy', 'wfc_cutoff', 'rho_cutoff']:
            header[key] = float(header[key])
        for key in ['l_max', 'l_max_rho', 'mesh_size', 'number_of_wfc',
                    'number_of_proj']:
            header[key] = int(header[key])
    else:
        assert len(header) == 0
        # This is the crappy format.
        headerlines = [line for line
                       in header_element.text.strip().splitlines()]
        #header['version'] = headerlines[0].split()[0]
        header['element'] = headerlines[1].split()[0]
        #header['has_nlcc'] = {'T': True, 'F': False} .........
        #assert not header['has_nlcc']
        # XC functional?
        header['z_valence'] = float(headerlines[5].split()[0])
        header['total_psenergy'] = float(headerlines[6].split()[0])
        # cutoffs?
        header['l_max'] = int(headerlines[8].split()[0])
        header['mesh_size'] = int(headerlines[9].split()[0])
        nwfs, nprojectors = headerlines[10].split()[:2]
        header['number_of_wfc'] = int(nwfs)
        header['number_of_proj'] = int(nprojectors)

    info_element = root.find('PP_INFO')
    pp['info'] = info_element.text

    pp['header'] = header
    mesh_element = root.find('PP_MESH')
    pp['r'] = toarray(mesh_element.find('PP_R'))
    pp['rab'] = toarray(mesh_element.find('PP_RAB'))
    
    # Convert to Hartree from Rydberg.
    pp['vlocal'] = 0.5 * toarray(root.find('PP_LOCAL'))
    
    nonlocal = root.find('PP_NONLOCAL')
    
    pp['projectors'] = []
    for element in nonlocal:
        if element.tag.startswith('PP_BETA'):
            if '.' in element.tag:
                name, num = element.tag.split('.')
                assert name == 'PP_BETA'
                attr = element.attrib
                proj = UPFStateSpec(int(attr['index']),
                                    attr['label'],
                                    int(attr['angular_momentum']),
                                    toarray(element))
                assert num == attr['index']
            else:
                tokens = element.text.split()
                metadata = tokens[:5]
                values = map(float, tokens[5:])
                
                npts = int(metadata[4])
                assert npts == len(values), (npts, len(values))
                while values[-1] == 0.0:
                    values.pop()

                proj = UPFStateSpec(int(tokens[0]), '??', int(tokens[1]),
                                    np.array(values))

            pp['projectors'].append(proj)
        else:
            if v201:
                assert element.tag == 'PP_DIJ', element.tag
                # XXX probably measured in Rydberg.
                pp['DIJ'] = 0.5 * toarray(element)
            else:
                lines = element.text.strip().splitlines()
                nnonzero = int(lines[0].split()[0])
                nproj = pp['header']['number_of_proj']
                D_ij = np.zeros((nproj, nproj))
                for n in range(1, nnonzero + 1):
                    tokens = lines[n].split()
                    i = int(tokens[0]) - 1
                    j = int(tokens[1]) - 1
                    D = float(tokens[2])
                    D_ij[i, j] = D
                    D_ij[j, i] = D
                assert len(lines) == 1 + nnonzero
                # XXX probably measured in Rydberg.
                pp['DIJ'] = 0.5 * D_ij

    pswfc_element = root.find('PP_PSWFC')
    pp['states'] = []
    if v201:
        for element in pswfc_element:
            attr = element.attrib
            name = element.tag
            state = UPFStateSpec(int(attr['index']),
                                 attr['label'],
                                 int(attr['l']),
                                 toarray(element),
                                 float(attr['occupation']),
                                 int(attr['n']))
            pp['states'].append(state)
    else:
        state_data = []
        for line in pswfc_element.text.splitlines()[1:]:
            if line.endswith('Wavefunction'):
                values = []
                state_data.append((line, values))
            else:
                values.extend(map(float, line.split()))
        for header, values in state_data:
            label, l, occupation, wfstring = header.split()
            assert wfstring == 'Wavefunction'
            state = UPFStateSpec(None, label, int(l), np.array(values),
                                 occupation=float(occupation))
            pp['states'].append(state)
        #print repr(lines[0])
    #assert len(pp['states']) > 0

    pp['rhoatom'] = toarray(root.find('PP_RHOATOM'))
    return pp


class UPFSetupData:
    def __init__(self, data):
        # data can be string (filename)
        # or dict (that's what we are looking for).
        # Maybe just a symbol would also be fine if we know the
        # filename to look for.
        if isinstance(data, basestring):
            data = parse_upf(data)
        
        assert isinstance(data, dict)
        self.data = data # more or less "raw" data from the file

        self.name = 'upf'

        header = data['header']
        
        keys = header.keys()
        keys.sort()

        self.symbol = header['element']
        self.Z = atomic_numbers[self.symbol]
        self.Nv = header['z_valence']
        self.Nc = self.Z - self.Nv
        self.Delta0 = -self.Nv / np.sqrt(4.0 * np.pi) # like hgh
        self.rcut_j = [data['r'][len(proj.values) - 1]
                       for proj in data['projectors']]

        #beta = 0.1 # XXX nice parameters?
        #N = 4 * 450 # XXX
        # This is "stolen" from hgh.  Figure out something reasonable
        #rgd = AERadialGridDescriptor(beta / N, 1.0 / N, N,
        #                             default_spline_points=100)
        
        from gpaw.atom.radialgd import EquidistantRadialGridDescriptor
        rgd = EquidistantRadialGridDescriptor(0.02)
        self.rgd = rgd
        
        # Whyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy???
        # What abominable part of the code requires the states
        # to be ordered like this?
        self.jargs = self.get_jargs()

        projectors = [data['projectors'][jarg] for jarg in self.jargs]
        if projectors:
            self.l_j = [proj.l for proj in projectors]
            self.pt_jg = []
            for proj in projectors:
                val = proj.values.copy()
                val /= 2.
                val[1:] /= data['r'][1:len(val)]
                val[0] = val[1]
                pt_g = self._interp(val) #* np.sqrt(4.0 * np.pi)
                #sqrnorm = (pt_g**2 * self.rgd.dr_g).sum()
                self.pt_jg.append(pt_g)
        
        else:
            self.l_j = [0]
            gcut = self.rgd.r2g(1.0)
            self.pt_jg = [np.zeros(gcut)] # XXX yet another "null" function

        self.rcgauss = 0.0 # XXX .... what is this used for?
        self.ni = sum([2 * l + 1 for l in self.l_j])

        self.filename = None # remember filename?
        self.fingerprint = None # XXX hexdigest the file?
        self.HubU = None # XXX
        self.lq = None # XXX

        #if data['states']:
        #    states_lmax = max([state.l for state in data['states']])
        #else:
        #    states_lmax = 1 # XXXX

        if data['states']:
            states_lmax = max([state.l for state in data['states']])
            f_ln = [[] for _ in range(1 + states_lmax)]
            electroncount = 0.0
            for state in data['states']:
                # Where should the electrons be in the inner list??
                # This is probably wrong and will lead to bad initialization
                f_ln[state.l].append(state.occupation)
                electroncount += state.occupation
                # The Cl.pz-hgh.UPF from quantum espresso has only 6
                # but should have 7 electrons.  Oh well....
            #err = abs(electroncount - self.Nv)
            self.f_j = [state.occupation for state in data['states']]
            self.n_j = [state.n for state in data['states']]
            self.l_orb_j = [state.l for state in data['states']]
            self.f_ln = f_ln
        else:
            self.n_j, self.l_orb_j, self.f_j, self.f_ln = \
                figure_out_valence_states(self)

        vlocal_unscreened = data['vlocal']

        # The UPF representation of HGH setups should be equal to that
        # used with setups='hgh'.  But the UPF files do not contain
        # info on the localization of the compensation charges!  That
        # means the screen_potential code will choose the shape of the
        # compensation charges, resulting in some numerical differences.
        #
        # One could hack that it looks up some radii to compare e.g. H2O.
        # The results of the below comments still have some numerical error.
        # But it's not huge.
        #a = None
        #if self.symbol == 'H':
        #    a = 0.2
        #elif self.symbol == 'O':
        #    a = .247621

        vbar_g, ghat_g = screen_potential(data['r'], vlocal_unscreened,
                                          self.Nv) #, a=a)
        
        self.vbar_g = self._interp(vbar_g) * np.sqrt(4.0 * np.pi)
        self.ghat_lg = [4.0 * np.pi / self.Nv * self._interp(ghat_g)]

        # XXX Subtract Hartree energy of compensation charge as reference

    def get_jargs(self):
        projectors = list(self.data['projectors'])
        jargs = []
        
        for n in range(4):
            for l in range(4):
                for i, proj in enumerate(projectors):
                    if proj.l == l:
                        jargs.append(proj.index - 1)
                        projectors.pop(i)
                        break
        return jargs

    def tostring(self):
        lines = []
        indent = 0
        
        def add(line):
            lines.append(indent * ' ' + line)

        add('Norm-conserving UPF setup:')
        indent = 2
        add('Element: %4s' % self.symbol)
        add('Z:       %4s' % self.Z)
        add('Valence: %4s' % self.Nv)
        indent -= 2
        add('Projectors:')
        indent += 2
        for j, proj in enumerate(self.data['projectors']):
            add('l=%d rcut=%s' % (proj.l, self.rcut_j[j]))
        indent -= 2
        if len(self.data['states']) == 0:
            add('No states stored on this setup')
        else:
            add('States:')
            indent += 2
            for j, state in enumerate(self.data['states']):
                add('l=%d f=%s' % (state.l, state.occupation))
            indent -= 2
        add('Local potential cutoff: %s'
            % self.get_local_potential().get_cutoff())
        add('Comp charge cutoff:     %s'
            % self.rgd.r_g[len(self.ghat_lg[0]) - 1])
        add('')
        return '\n'.join(lines)

    def print_info(self, txt, setup):
        txt(self.tostring())

    # XXX something can perhaps be stolen from HGH
    def expand_hamiltonian_matrix(self):
        """Construct K_p from individual h_nn for each l."""
        ni = sum([2 * l + 1 for l in self.l_j])
        nj = len(self.l_j)

        H_ii = np.zeros((ni, ni))
        if len(self.data['DIJ']) == 0:
            return pack2(H_ii)
        
        # Multiply by 4.
        # I think the factor of 4 compensates for the fact that the projectors
        # all had square norms of 4, but we brought them back down to 1
        # because that's more sensible.
        H_jj = 4.0 * self.data['DIJ'].reshape((nj, nj))
        m1start = 0
        for j1, l1 in enumerate(self.l_j):
            j1 = self.jargs[j1]
            m1stop = m1start + 2 * l1 + 1
            m2start = 0
            for j2, l2 in enumerate(self.l_j):
                j2 = self.jargs[j2]
                m2stop = m2start + 2 * l2 + 1
                if l1 == l2:
                    dm = m1stop - m1start
                    H_ii[m1start:m1stop, m2start:m2stop] = \
                        np.eye(dm) * H_jj[j1, j2]
                else:
                    assert H_jj[j1, j2] == 0.0
                m2start = m2stop
            m1start = m1stop
        return pack2(H_ii)
    
    def get_local_potential(self):
        vbar = Spline(0, self.rgd.r_g[len(self.vbar_g) - 1], self.vbar_g)
        return vbar

    # XXXXXXXXXXXXXXXXX stolen from hghsetupdata
    def get_projectors(self):
        # XXX equal-range projectors still required for some reason
        maxlen = max([len(pt_g) for pt_g in self.pt_jg])
        pt_j = []
        for l, pt1_g in zip(self.l_j, self.pt_jg):
            pt2_g = self.rgd.zeros()[:maxlen]
            pt2_g[:len(pt1_g)] = divrl(pt1_g, l, self.rgd.r_g[:len(pt1_g)])
            spline = Spline(l, self.rgd.r_g[maxlen - 1], pt2_g)
            pt_j.append(spline)
        return pt_j

    def _interp(self, func):
        r_g = self.rgd.r_g
        gcut = len(func)
        rcut = self.data['r'][gcut - 1]
        gcutnew = np.searchsorted(self.rgd.r_g, rcut)
        rnew = r_g[:gcutnew]
        ynew = np.interp(rnew, self.data['r'][:gcut], func, right=0.0)
        return ynew

    def get_compensation_charge_functions(self):
        assert len(self.ghat_lg) == 1
        ghat_g = self.ghat_lg[0]
        ng = len(ghat_g)
        rcutcc = self.rgd.r_g[ng - 1] # correct or not?
        r = np.linspace(0.0, rcutcc, 50)
        ghat_g[-1] = 0.0
        ghatnew_g = Spline(0, rcutcc, ghat_g).map(r)
        return r, [0], [ghatnew_g]

    def create_basis_functions(self):
        if len(self.data['states']) > 0:
            return self.get_stored_basis_functions()
        else:
            from gpaw.pseudopotential import generate_basis_functions
            return generate_basis_functions(self)
    
    def get_stored_basis_functions(self, ):
        b = Basis(self.symbol, 'upf', readxml=False)
        b.generatordata = 'upf-pregenerated'
        
        states = self.data['states']
        maxlen = max([len(state.values) for state in states])
        orig_r = self.data['r']
        rcut = min(orig_r[maxlen - 1], 12.0) # XXX hardcoded 12 max radius
        
        b.d = 0.02
        b.ng = int(1 + rcut / b.d)
        rgd = b.get_grid_descriptor()
        
        for j, state in enumerate(states):
            val = state.values
            phit_g = np.interp(rgd.r_g, orig_r, val)
            phit_g = divrl(phit_g, 1, rgd.r_g)
            icut = len(phit_g) - 1 # XXX correct or off-by-one?
            rcut = rgd.r_g[icut]
            bf = BasisFunction(state.l, rcut, phit_g, 'pregenerated')
            b.bf_j.append(bf)
        return b

    def build(self, xcfunc, lmax, basis, filter=None):
        if basis is None:
            basis = self.create_basis_functions()
        return PseudoPotential(self, basis)


def main_plot():
    p = OptionParser(usage='%prog [OPTION] [FILE...]',
                     description='plot upf pseudopotential from file.')
    opts, args = p.parse_args()

    import pylab as pl
    for fname in args:
        pp = parse_upf(fname)
        print('--- %s ---' % fname)
        print(UPFSetupData(pp).tostring())
        print(pp['info'])
        upfplot(pp, show=False)
    pl.show()


def upfplot(setup, show=True):
    # A version of this, perhaps nicer, is in pseudopotential.py.
    # Maybe it is not worth keeping this version
    if isinstance(setup, dict):
        setup = UPFSetupData(setup)
    pp = setup.data
    r0 = pp['r'].copy()
    r0[0] = 1e-8
    
    def rtrunc(array, rdividepower=0):
        r = r0[:len(array)]
        arr = divrl(array, rdividepower, r)
        return r, arr
    
    import pylab as pl
    fig = pl.figure()
    fig.canvas.set_window_title('%s - UPF setup for %s' % (pp['fname'],
                                                           setup.symbol))

    vax = fig.add_subplot(221)
    pax = fig.add_subplot(222)
    rhoax = fig.add_subplot(223)
    wfsax = fig.add_subplot(224)

    r, v = rtrunc(pp['vlocal'])
    
    vax.plot(r, v, label='vloc')

    vscreened, rhocomp = screen_potential(r, v, setup.Nv)
    icut = len(rhocomp)
    vcomp = v.copy()
    vcomp[:icut] -= vscreened
    vax.axvline(r[icut], ls=':', color='k')

    vax.plot(r, vcomp, label='vcomp')
    vax.plot(r[:icut], vscreened, label='vscr')
    vax.axis(xmin=0.0, xmax=6.0)
    rhoax.plot(r[:icut], rhocomp[:icut], label='rhocomp')

    for j, proj in enumerate(pp['projectors']):
        r, p = rtrunc(proj.values, 0)
        pax.plot(r, p,
                 label='p%d [l=%d]' % (j + 1, proj.l))

    for j, st in enumerate(pp['states']):
        r, psi = rtrunc(st.values, 1)
        wfsax.plot(r, psi, label='wf%d %s' % (j, st.label))

    r, rho = rtrunc(pp['rhoatom'], 2)
    wfsax.plot(r, rho, label='rho')

    vax.legend(loc='best')
    rhoax.legend(loc='best')
    pax.legend(loc='best')
    wfsax.legend(loc='best')

    for ax in [vax, rhoax, pax, wfsax]:
        ax.set_xlabel('r [Bohr]')
        ax.axhline(0.0, ls=':', color='black')

    vax.set_ylabel('potential')
    pax.set_ylabel('projectors')
    wfsax.set_ylabel('WF / density')
    rhoax.set_ylabel('Comp charges')

    fig.subplots_adjust(wspace=0.3)

    if show:
        pl.show()
