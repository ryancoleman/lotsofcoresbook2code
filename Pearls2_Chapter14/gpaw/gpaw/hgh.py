import numpy as np
from ase.data import atomic_numbers

from gpaw.utilities import pack2
from gpaw.utilities.tools import md5_new
from gpaw.atom.radialgd import AERadialGridDescriptor
from gpaw.atom.configurations import configurations
from gpaw.pseudopotential import PseudoPotential

setups = {}  # Filled out during parsing below
sc_setups = {}  # Semicore


# Tabulated values of Gamma(m + 1/2)
half_integer_gamma = [np.sqrt(np.pi)]
for m in range(20):
    half_integer_gamma.append(half_integer_gamma[m] * (m + 0.5))


class HGHSetupData:
    """Setup-compatible class implementing HGH pseudopotential.

    To the PAW code this will appear as a legit PAW setup, but is
    in fact considerably simpler.  In particular, all-electron and
    pseudo partial waves are all zero, and compensation charges do not
    depend on environment.

    A HGH setup has the following form::

                  ----
                   \
      V = Vlocal +  )  | p  > h   < p |
                   /      i    ij    j
                  ----
                   ij

    Vlocal contains a short-range term which is Gaussian-shaped and
    implemented as vbar of a PAW setup, along with a long-range term
    which goes like 1/r and is implemented in terms of a compensation
    charge.

    The non-local part contains KB projector functions which are
    essentially similar to those in PAW, while h_ij are constants.
    h_ij are provided by setting the K_p variable of the normal
    setup.

    Most other properties of a PAW setup do not exist for HGH setups, for
    which reason they are generally set to zero:

    * All-electron partial waves: always zero
    * Pseudo partial waves: always zero
    * Projectors: HGH projectors
    * Zero potential (vbar): Gaussian times polynomial
    * Compensation charges: One Gaussian-shaped spherically symmetric charge
    * All-electron core density: Delta function corresponding to core electron
      charge
    * Pseudo core density: always zero
    * Pseudo valence density: always zero
    * PS/AE Kinetic energy density: always zero
    * The mysterious constants K_p of a setup correspond to h_ij.

    Note that since the pseudo partial waves are set to zero,
    initialization of atomic orbitals requires loading a custom basis
    set.

    Absolute energies become numerically large since no atomic
    reference is subtracted.
    """
    def __init__(self, hghdata):
        if isinstance(hghdata, str):
            symbol = hghdata
            if symbol.endswith('.sc'):
                hghdata = sc_setups[symbol[:-3]]
            else:
                hghdata = setups[symbol]
        self.hghdata = hghdata
        
        chemsymbol = hghdata.symbol
        if '.' in chemsymbol:
            chemsymbol, sc = chemsymbol.split('.')
            assert sc == 'sc'
        self.symbol = chemsymbol
        self.type = hghdata.symbol
        self.name = 'LDA'
        self.initialize_setup_data()

    def initialize_setup_data(self):
        hghdata = self.hghdata
        beta = 0.1
        N = 450
        rgd = AERadialGridDescriptor(beta / N, 1.0 / N, N,
                                     default_spline_points=100)
        #rgd = EquidistantRadialGridDescriptor(0.001, 10000)
        self.rgd = rgd

        self.Z = hghdata.Z
        self.Nc = hghdata.Z - hghdata.Nv
        self.Nv = hghdata.Nv
        self.rcgauss = np.sqrt(2.0) * hghdata.rloc

        threshold = 1e-8
        if len(hghdata.c_n) > 0:
            vloc_g = create_local_shortrange_potential(rgd.r_g, hghdata.rloc,
                                                       hghdata.c_n)
            gcutvbar, rcutvbar = self.find_cutoff(rgd.r_g, rgd.dr_g, vloc_g,
                                                  threshold)
            self.vbar_g = np.sqrt(4.0 * np.pi) * vloc_g[:gcutvbar]
        else:
            rcutvbar = 0.5
            gcutvbar = rgd.ceil(rcutvbar)
            self.vbar_g = np.zeros(gcutvbar)

        nj = sum([v.nn for v in hghdata.v_l])
        if nj == 0:
            nj = 1  # Code assumes nj > 0 elsewhere, we fill out with zeroes

        if not hghdata.v_l:
            # No projectors.  But the remaining code assumes that everything
            # has projectors!  We'll just add the zero function then
            hghdata.v_l = [VNonLocal(0, 0.01, [[0.]])]

        n_j = []
        l_j = []

        # j ordering is significant, must be nl rather than ln
        for n, l in self.hghdata.nl_iter():
            n_j.append(n + 1)  # Note: actual n must be positive!
            l_j.append(l)
        assert nj == len(n_j)
        self.nj = nj
        self.l_j = l_j
        self.l_orb_j = l_j
        self.n_j = n_j

        self.rcut_j = []
        self.pt_jg = []

        for n, l in zip(n_j, l_j):
            # Note: even pseudopotentials without projectors will get one
            # projector, but the coefficients h_ij should be zero so it
            # doesn't matter
            pt_g = create_hgh_projector(rgd.r_g, l, n, hghdata.v_l[l].r0)
            norm = np.sqrt(np.dot(rgd.dr_g, pt_g**2 * rgd.r_g**2))
            assert np.abs(1 - norm) < 1e-5, str(1 - norm)
            gcut, rcut = self.find_cutoff(rgd.r_g, rgd.dr_g, pt_g, threshold)
            if rcut < 0.5:
                rcut = 0.5
                gcut = rgd.ceil(rcut)
            pt_g = pt_g[:gcut].copy()
            rcut = max(rcut, 0.5)
            self.rcut_j.append(rcut)
            self.pt_jg.append(pt_g)

        # This is the correct magnitude of the otherwise normalized
        # compensation charge
        self.Delta0 = -self.Nv / np.sqrt(4.0 * np.pi)

        f_ln = self.hghdata.get_occupation_numbers()
        f_j = [0] * nj
        for j, (n, l) in enumerate(self.hghdata.nl_iter()):
            try:
                f_j[j] = f_ln[l][n]
            except IndexError:
                pass
        self.f_ln = f_ln
        self.f_j = f_j

    def find_cutoff(self, r_g, dr_g, f_g, sqrtailnorm=1e-5):
        g = len(r_g)
        acc_sqrnorm = 0.0
        while acc_sqrnorm <= sqrtailnorm:
            g -= 1
            acc_sqrnorm += (r_g[g] * f_g[g])**2.0 * dr_g[g]
            if r_g[g] < 0.5:  # XXX
                return g, r_g[g]
        return g, r_g[g]

    def expand_hamiltonian_matrix(self):
        """Construct K_p from individual h_nn for each l."""
        ni = sum([2 * l + 1 for l in self.l_j])

        H_ii = np.zeros((ni, ni))

        # The H_ii used in gpaw is much larger and more general than the one
        # required for HGH pseudopotentials.  This means a lot of the elements
        # must be assigned the same value.  Not a performance issue though,
        # since these are small matrices
        M1start = 0
        for n1, l1 in self.hghdata.nl_iter():
            M1end = M1start + 2 * l1 + 1
            M2start = 0
            v = self.hghdata.v_l[l1]
            for n2, l2 in self.hghdata.nl_iter():
                M2end = M2start + 2 * l2 + 1
                if l1 == l2:
                    h_nn = v.expand_hamiltonian_diagonal()
                    H_mm = np.identity(M2end - M2start) * h_nn[n1, n2]
                    H_ii[M1start:M1end, M2start:M2end] += H_mm
                M2start = M2end
            M1start = M1end
        K_p = pack2(H_ii)
        return K_p

    def __str__(self):
        return "HGHSetupData('%s')" % self.type

    def __repr__(self):
        return self.__str__()

    def print_info(self, text, _setup):
        self.hghdata.print_info(text)

    def plot(self):
        """Plot localized functions of HGH setup."""
        import pylab as pl
        rgd = self.rgd

        pl.subplot(211)  # vbar, compensation charge
        gcutvbar = len(self.vbar_g)
        pl.plot(rgd.r_g[:gcutvbar], self.vbar_g, 'r', label='vloc',
                linewidth=3)
        rcc, gcc = self.get_compensation_charge_functions()
        gcc = gcc[0]

        pl.plot(rcc, gcc * self.Delta0, 'b--', label='Comp charge [arb. unit]',
                linewidth=3)
        pl.legend(loc='best')

        pl.subplot(212)  # projectors
        for j, (n, l, pt_g) in enumerate(zip(self.n_j, self.l_j, self.pt_jg)):
            label = 'n=%d, l=%d' % (n, l)
            pl.ylabel('$p_n^l(r)$')
            ng = len(pt_g)
            r_g = rgd.r_g[:ng]
            pl.plot(r_g, pt_g, label=label)
        pl.legend()

    def get_projectors(self):
        # XXX equal-range projectors still required for some reason
        maxlen = max([len(pt_g) for pt_g in self.pt_jg])
        pt_j = []
        for l, pt1_g in zip(self.l_j, self.pt_jg):
            pt2_g = self.rgd.zeros()[:maxlen]
            pt2_g[:len(pt1_g)] = pt1_g
            pt_j.append(self.rgd.spline(pt2_g, self.rgd.r_g[maxlen - 1], l))
        return pt_j

    def create_basis_functions(self):
        from gpaw.pseudopotential import generate_basis_functions
        return generate_basis_functions(self)

    def get_compensation_charge_functions(self):
        alpha = self.rcgauss**-2
        
        rcutgauss = self.rcgauss * 5.0
        # smaller values break charge conservation
        
        r = np.linspace(0.0, rcutgauss, 100)
        g = alpha**1.5 * np.exp(-alpha * r**2) * 4.0 / np.sqrt(np.pi)
        g[-1] = 0.0
        return r, [0], [g]

    def get_local_potential(self):
        n = len(self.vbar_g)
        return self.rgd.spline(self.vbar_g, self.rgd.r_g[n - 1])

    def build(self, xcfunc, lmax, basis, filter=None):
        if basis is None:
            basis = self.create_basis_functions()
        setup = PseudoPotential(self, basis)
        setup.fingerprint = md5_new(str(self.hghdata)).hexdigest()
        return setup


def create_local_shortrange_potential(r_g, rloc, c_n):
    rr_g = r_g / rloc  # "Relative r"
    rr2_g = rr_g**2
    rr4_g = rr2_g**2
    rr6_g = rr4_g * rr2_g

    gaussianpart = np.exp(-.5 * rr2_g)
    polypart = np.zeros(r_g.shape)
    for c, rrn_g in zip(c_n, [1, rr2_g, rr4_g, rr6_g]):
        polypart += c * rrn_g

    vloc_g = gaussianpart * polypart
    return vloc_g


def create_hgh_projector(r_g, l, n, r0):
    poly_g = r_g**(l + 2 * (n - 1))
    gauss_g = np.exp(-.5 * r_g**2 / r0**2)
    A = r0**(l + (4 * n - 1) / 2.0)
    assert (4 * n - 1) % 2 == 1
    B = half_integer_gamma[l + (4 * n - 1) // 2]**.5
    pt_g = 2.**.5 / A / B * poly_g * gauss_g
    return pt_g


# Coefficients determining off-diagonal elements of h_nn for l = 0...2
# given the diagonal elements
hcoefs_l = [
    [-.5 * (3. / 5.)**.5, .5 * (5. / 21.)**.5, -.5 * (100. / 63.)**.5],
    [-.5 * (5. / 7.)**.5, 1. / 6. * (35. / 11.)**.5, -1. / 6. * 14. / 11.**.5],
    [-.5 * (7. / 9.)**.5, .5 * (63. / 143)**.5, -.5 * 18. / 143.**.5]
    ]


class VNonLocal:
    """Wrapper class for one nonlocal term of an HGH potential."""
    def __init__(self, l, r0, h_n):
        self.l = l
        self.r0 = r0
        h_n = np.array(h_n)
        nn = len(h_n)
        self.nn = nn
        self.h_n = h_n

    def expand_hamiltonian_diagonal(self):
        """Construct full atomic Hamiltonian from diagonal elements."""
        nn = self.nn
        h_n = self.h_n
        h_nn = np.zeros((nn, nn))
        for n, h in enumerate(h_n):
            h_nn[n, n] = h
        if self.l > 2:
            #print 'Warning: no diagonal elements for l=%d' % l
            # Some elements have projectors corresponding to l=3, but
            # the HGH article only specifies how to calculate the
            # diagonal elements of the atomic hamiltonian for l = 0, 1, 2 !
            return
        coefs = hcoefs_l[self.l]
        if nn > 2:
            h_nn[0, 2] = h_nn[2, 0] = coefs[1] * h_n[2]
            h_nn[1, 2] = h_nn[2, 1] = coefs[2] * h_n[2]
        if nn > 1:
            h_nn[0, 1] = h_nn[1, 0] = coefs[0] * h_n[1]
        return h_nn

    def copy(self):
        return VNonLocal(self.l, self.r0, self.h_n.copy())

    def serialize(self):  # no spin-orbit part
        return ' '.join(['    ', '%-10s' % self.r0] +
                        ['%10f' % h for h in self.h_n])


class HGHParameterSet:
    """Wrapper class for HGH-specific data corresponding to one element."""
    def __init__(self, symbol, Z, Nv, rloc, c_n, v_l):
        self.symbol = symbol  # Identifier, e.g. 'Na', 'Na.sc', ...
        self.Z = Z  # Actual atomic number
        self.Nv = Nv  # Valence electron count
        self.rloc = rloc  # Characteristic radius of local part
        self.c_n = np.array(c_n)  # Polynomial coefficients for local part
        self.v_l = list(v_l)  # Non-local parts

        Z, nlfe_j = configurations[self.symbol.split('.')[0]]
        self.configuration = nlfe_j

    def __str__(self):
        strings = ['HGH setup for %s\n' % self.symbol,
                   '    Valence Z=%d, rloc=%.05f\n' % (self.Nv, self.rloc)]

        if len(self.c_n) > 0:
            coef_string = ', '.join(['%.05f' % c for c in self.c_n])
        else:
            coef_string = 'zeros'
        strings.append('    Local part coeffs: %s\n' % coef_string)
        strings.append('    Projectors:\n')
        if not self.v_l:
            strings.append('        None\n')
        for v in self.v_l:
            strings.append('        l=%d, rc=%.05f\n' % (v.l, v.r0))
        strings.append('    Diagonal coefficients of nonlocal parts:')
        if not self.v_l:
            strings.append('\n        None\n')
        for v in self.v_l:
            strings.append('\n')
            strings.append('        l=%d: ' % v.l +
                           ', '.join(['%8.05f' % h for h in v.h_n]))
        return ''.join(strings)

    def copy(self):
        other = HGHParameterSet(self.symbol, self.Z, self.Nv, self.rloc,
                                self.c_n, self.v_l)
        return other

    def print_info(self, txt):
        txt(str(self))
        txt()

    def nl_iter(self):
        for n in range(4):
            for l, v in enumerate(self.v_l):
                if n < v.nn:
                    yield n, l

    def get_occupation_numbers(self):
        nlfe_j = list(self.configuration)
        nlfe_j.reverse()
        f_ln = [[], [], []]  # [[s], [p], [d]]
        # f states will be ignored as the atomic Hamiltonians
        # of those are, carelessly, not defined in the article.
        lmax = len(self.v_l) - 1
        Nv = 0
        # Right.  We need to find the occupation numbers of each state and
        # put them into a nice list of lists f_ln.
        #
        # We loop over states starting with the least bound one
        # (i.e. reversed nlfe_j), adding the occupation numbers of each state
        # as appropriate.  Once we have the right number of electrons, we
        # end the loop.
        #
        # Some states in the standard configuration might
        # be f-type; these should be skipped (unless the HGH setup actually
        # has a valence f-state; however as noted above, some of the
        # parameters are undefined in that case so are ignored anyway).  More
        # generally if for some state l > lmax,
        # we can skip that state.
        for n, l, f, e in nlfe_j:
            if l > lmax:
                continue
            Nv += f
            f_n = f_ln[l]
            assert f_n == [] or self.symbol.endswith('.sc')
            f_n.append(f)
            if Nv >= self.Nv:
                break
        assert Nv == self.Nv
        return f_ln

    def zeropad(self):
        """Return a new HGHParameterSet with all arrays zero padded so they
        have the same (max) length for all such HGH setups.  Makes
        plotting multiple HGH setups easier because they have compatible
        arrays."""
        c_n = np.zeros(4)
        for n, c in enumerate(self.c_n):
            c_n[n] = c
        v_l = []
        for l, v in enumerate(self.v_l):
            h_n = np.zeros(3)
            h_n[:len(v.h_n)] = list(v.h_n)
            v2 = VNonLocal(l, v.r0, h_n)
            v_l.append(v2)
        for l in range(len(self.v_l), 3):
            v_l.append(VNonLocal(l, 0.5, np.zeros(3)))
        copy = HGHParameterSet(self.symbol, self.Z, self.Nv, self.rloc, c_n,
                               v_l)
        return copy

    def serialize(self):
        string1 = '%-5s %-12s %10s ' % (self.symbol, self.Z, self.rloc)
        string2 = ' '.join(['%.10s' % c for c in self.c_n])
        nonlocal_strings = [v.serialize() for v in self.v_l]
        return '\n'.join([string1 + string2] + nonlocal_strings)


def parse_local_part(string):
    """Create HGHParameterSet object with local part initialized."""
    tokens = iter(string.split())
    symbol = tokens.next()
    actual_chemical_symbol = symbol.split('.')[0]
    Z = atomic_numbers[actual_chemical_symbol]
    Nv = int(tokens.next())
    rloc = float(tokens.next())
    c_n = [float(token) for token in tokens]
    return symbol, Z, Nv, rloc, c_n


class HGHBogusNumbersError(ValueError):
    """Error which is raised when the HGH parameters contain f-type
    or higher projectors.  The HGH article only defines atomic Hamiltonian
    matrices up to l=2, so these are meaningless."""
    pass


def parse_hgh_setup(lines):
    """Initialize HGHParameterSet object from text representation."""
    lines = iter(lines)
    symbol, Z, Nv, rloc, c_n = parse_local_part(lines.next())

    def pair_up_nonlocal_lines(lines):
        yield lines.next(), ''
        while True:
            yield lines.next(), lines.next()

    v_l = []
    for l, (nonlocal, spinorbit) in enumerate(pair_up_nonlocal_lines(lines)):
        # we discard the spinorbit 'k_n' data so far
        nltokens = nonlocal.split()
        r0 = float(nltokens[0])
        h_n = [float(token) for token in nltokens[1:]]

        #if h_n[-1] == 0.0: # Only spin-orbit contributes.  Discard.
        #    h_n.pop()
        # Actually the above causes trouble.  Probably it messes up state
        # ordering or something else that shouldn't have any effect.
        
        vnl = VNonLocal(l, r0, h_n)
        v_l.append(vnl)
        if l > 2:
            raise HGHBogusNumbersError

    hgh = HGHParameterSet(symbol, Z, Nv, rloc, c_n, v_l)
    return hgh


def str2hgh(string):
    return parse_hgh_setup(string.splitlines())


def hgh2str(hgh):
    return hgh.serialize()


def parse_setups(lines):
    """Read HGH data from file."""
    setups = {}
    entry_lines = [i for i in xrange(len(lines))
                   if lines[i][0].isalpha()]
    lines_by_element = [lines[entry_lines[i]:entry_lines[i + 1]]
                        for i in xrange(len(entry_lines) - 1)]
    lines_by_element.append(lines[entry_lines[-1]:])

    for elines in lines_by_element:
        try:
            hgh = parse_hgh_setup(elines)
        except HGHBogusNumbersError:
            continue
        assert hgh.symbol not in setups
        setups[hgh.symbol] = hgh
    return setups


def plot(symbol, extension=None):
    import pylab as pl
    try:
        s = HGHSetupData(symbol)
    except IndexError:
        print('Nooooo')
        return
    s.plot()
    if extension is not None:
        pl.savefig('hgh.%s.%s' % (symbol, extension))


def plot_many(*symbols):
    import pylab as pl
    if not symbols:
        symbols = setups.keys() + [key + '.sc' for key in sc_setups.keys()]
    for symbol in symbols:
        pl.figure(1)
        plot(symbol, extension='png')
        pl.clf()


def parse_default_setups():
    from hgh_parameters import parameters
    lines = parameters.splitlines()
    setups0 = parse_setups(lines)
    for key, value in setups0.items():
        if key.endswith('.sc'):
            sym, sc = key.split('.')
            sc_setups[sym] = value
        else:
            setups[key] = value

parse_default_setups()
