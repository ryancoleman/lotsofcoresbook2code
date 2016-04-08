"""This module is used to generate atomic orbital basis sets."""
from __future__ import print_function
import sys
from StringIO import StringIO

import numpy as np
from ase.units import Hartree
from ase.utils import devnull

from gpaw.version import version
from gpaw.utilities import divrl
from gpaw.atom.generator import Generator
from gpaw.atom.all_electron import AllElectron
from gpaw.atom.configurations import parameters
from gpaw.basis_data import Basis, BasisFunction
from gpaw.atom.radialgd import AERadialGridDescriptor


def get_basis_l(f_j, l_j, Nv):
    N = 0
    lvalues = []
    j_l = []
    Nl_j = len(l_j)
    while N < Nv:
        Nl_j -= 1
        if f_j[Nl_j] > 0:
            lvalues.append(l_j[Nl_j])
            N += f_j[Nl_j]
            j_l.append(Nl_j)
    assert N == Nv

    if min(lvalues) != 0:  # always include s-orbital !
        lvalues = lvalues + [0]
        reversed_l_j = list(l_j)
        reversed_l_j.reverse()
        j = len(reversed_l_j) - reversed_l_j.index(0) - 1
        j_l.append(j)

    return [j_l, lvalues]


def get_gaussianlike_basis_function(rgd, l, rchar, gcut):
    rcut = rgd.r_g[gcut]
    gaussian = QuasiGaussian(1.0 / rchar**2, rcut)
    r_g = rgd.r_g[:gcut + 1]
    phit_g = gaussian(r_g) * r_g**l
    phit_g[-1] = 0.0
    norm = np.dot(rgd.dr_g[:gcut + 1], (r_g * phit_g)**2)**0.5
    phit_g /= norm
    norm2 = np.dot(rgd.dr_g[:gcut + 1], (r_g * phit_g)**2)**0.5
    assert abs(1.0 - norm2) < 1e-10
    return phit_g


def make_split_valence_basis_function(r_g, psi_g, l, gcut):
    """Get polynomial which joins psi smoothly at rcut.

    Returns an array of function values f(r) * r, where::

              l           2
      f(r) = r  * (a - b r ),  r < rcut
      f(r) = psi(r),           r >= rcut

    where a and b are determined such that f(r) is continuous and
    differentiable at rcut.  The parameter psi should be an atomic
    orbital.
    """
    r1 = r_g[gcut]  # ensure that rcut is moved to a grid point
    r2 = r_g[gcut + 1]
    y1 = psi_g[gcut] / r_g[gcut]
    y2 = psi_g[gcut + 1] / r_g[gcut + 1]
    b = - (y2 / r2**l - y1 / r1**l) / (r2**2 - r1**2)
    a = (y1 / r1**l + b * r1**2)
    psi_g2 = r_g**(l + 1) * (a - b * r_g**2)
    psi_g2[gcut:] = psi_g[gcut:]
    return psi_g2


def rsplit_by_norm(rgd, l, u, tailnorm_squared, txt):
    """Find radius outside which remaining tail has a particular norm."""
    norm_squared = np.dot(rgd.dr_g, u * u)
    partial_norm_squared = 0.
    i = len(u) - 1
    absolute_tailnorm_squared = tailnorm_squared * norm_squared
    while partial_norm_squared < absolute_tailnorm_squared:
        # Integrate backwards.  This is important since the pseudo
        # wave functions have strange behaviour near the core.
        partial_norm_squared += rgd.dr_g[i] * u[i]**2
        i -= 1
    rsplit = rgd.r_g[i + 1]
    msg = ('Tail norm %.03f :: rsplit=%.02f Bohr' %
           ((partial_norm_squared / norm_squared)**0.5, rsplit))
    print(msg, file=txt)
    gsplit = rgd.floor(rsplit)
    splitwave = make_split_valence_basis_function(rgd.r_g, u, l, gsplit)
    return rsplit, partial_norm_squared, splitwave


class BasisMaker:
    """Class for creating atomic basis functions."""
    def __init__(self, generator, name=None, run=True, gtxt='-',
                 non_relativistic_guess=False, xc='PBE'):
        if isinstance(generator, str):  # treat 'generator' as symbol
            generator = Generator(generator, scalarrel=True,
                                  xcname=xc, txt=gtxt,
                                  nofiles=True)
            generator.N *= 4
        self.generator = generator
        self.rgd = AERadialGridDescriptor(generator.beta / generator.N,
                                          1.0 / generator.N, generator.N,
                                          default_spline_points=100)
        self.name = name
        if run:
            if non_relativistic_guess:
                ae0 = AllElectron(generator.symbol, scalarrel=False,
                                  nofiles=False, txt=gtxt, xcname=xc)
                ae0.N = generator.N
                ae0.beta = generator.beta
                ae0.run()
                # Now files will be stored such that they can
                # automagically be used by the next run()
            generator.run(write_xml=False, use_restart_file=False,
                          **parameters[generator.symbol])

    def smoothify(self, psi_mg, l):
        """Generate pseudo wave functions from all-electron ones.

        The pseudo wave function is::
        
                                   ___
               ~                   \    /   ~             \    ~     ~
            | psi  > = | psi  > +   )  | | phi > - | phi > ) < p  | psi > ,
                 m          m      /__  \     i         i /     i      m
                                    i

        where the scalar products are found by solving::
        
                            ___
              ~             \     ~             ~     ~
            < p | psi  > =   )  < p  | phi  > < p  | psi  > .
               i     m      /__    i      j      j      m
                             j

        In order to ensure smoothness close to the core, the
        all-electron wave function and partial wave are then
        multiplied by a radial function which approaches 0 near the
        core, such that the pseudo wave function approaches::

                        ___
               ~        \      ~       ~     ~
            | psi  > =   )  | phi >  < p  | psi >    (for r << rcut),
                 m      /__      i      i      m
                         i

        which is exact if the projectors/pseudo partial waves are complete.
        """
        if psi_mg.ndim == 1:
            return self.smoothify(psi_mg[None], l)[0]

        g = self.generator
        u_ng = g.u_ln[l]
        q_ng = g.q_ln[l]
        s_ng = g.s_ln[l]

        Pi_nn = np.dot(g.r * q_ng, u_ng.T)
        Q_nm = np.dot(g.r * q_ng, psi_mg.T)
        Qt_nm = np.linalg.solve(Pi_nn, Q_nm)

        # Weight-function for truncating all-electron parts smoothly near core
        gmerge = g.r2g(g.rcut_l[l])
        w_g = np.ones(g.r.shape)
        w_g[0:gmerge] = (g.r[0:gmerge] / g.r[gmerge])**2.
        w_g = w_g[None]
        
        psit_mg = psi_mg * w_g + np.dot(Qt_nm.T, s_ng - u_ng * w_g)
        return psit_mg

    def make_orbital_vector(self, j, rcut, vconf=None):
        """Returns a smooth basis vector given an all-electron one."""
        l = self.generator.l_j[j]
        psi_g, e = self.generator.solve_confined(j, rcut, vconf)
        psit_g = self.smoothify(psi_g, l)
        return psit_g

    def rcut_by_energy(self, j, esplit=.1, tolerance=.1, rguess=6.,
                       vconf_args=None):
        """Find confinement cutoff corresponding to given orbital energy shift.

        Creates a confinement potential for the orbital given by j,
        such that the confined-orbital energy is (emin to emax) eV larger
        than the free-orbital energy."""
        g = self.generator
        e_base = g.e_j[j]
        rc = rguess

        if vconf_args is None:
            vconf = None
        else:
            amplitude, ri_rel = vconf_args
            vconf = g.get_confinement_potential(amplitude, ri_rel * rc, rc)

        psi_g, e = g.solve_confined(j, rc, vconf)
        de_min, de_max = esplit / Hartree, (esplit + tolerance) / Hartree

        rmin = 0.
        rmax = g.r[-1]

        de = e - e_base
        #print '--------'
        #print 'Start bisection'
        #print'e_base =',e_base
        #print 'e =',e
        #print '--------'
        while de < de_min or de > de_max:
            if de < de_min:  # Move rc left -> smaller cutoff, higher energy
                rmax = rc
                rc = (rc + rmin) / 2.
            else:  # Move rc right
                rmin = rc
                rc = (rc + rmax) / 2.
            if vconf is not None:
                vconf = g.get_confinement_potential(amplitude, ri_rel * rc, rc)
            psi_g, e = g.solve_confined(j, rc, vconf)
            de = e - e_base
            #print 'rc = %.03f :: e = %.03f :: de = %.03f' % (rc, e*Hartree,
            #                                                 de*Hartree)
            #if rmin - rmax < 1e-
            if g.r2g(rmax) - g.r2g(rmin) <= 1:  # adjacent points
                break  # cannot meet tolerance due to grid resolution
        #print 'Done!'
        return psi_g, e, de, vconf, rc

    def generate(self, zetacount=2, polarizationcount=1,
                 tailnorm=(0.16, 0.3, 0.6), energysplit=0.1, tolerance=1.0e-3,
                 referencefile=None, referenceindex=None, rcutpol_rel=1.0,
                 rcutmax=20.0,
                 rcharpol_rel=None,
                 vconf_args=(12.0, 0.6), txt='-',
                 include_energy_derivatives=False,
                 #lvalues=None, # XXX clean up some of these!
                 jvalues=None,
                 l_pol=None):
        """Generate an entire basis set.

        This is a high-level method which will return a basis set
        consisting of several different basis vector types.

        Parameters:

        ===================== =================================================
        ``zetacount``         Number of basis functions per occupied orbital
        ``polarizationcount`` Number of polarization functions
        ``tailnorm``          List of tail norms for split-valence scheme
        ``energysplit``       Energy increase defining confinement radius (eV)
        ``tolerance``         Tolerance of energy split (eV)
        ``referencefile``     gpw-file used to generate polarization function
        ``referenceindex``    Index in reference system of relevant atom
        ``rcutpol_rel``       Polarization rcut relative to largest other rcut
        ``rcutmax``           No cutoff will be greater than this value
        ``vconf_args``        Parameters (alpha, ri/rc) for conf. potential
        ``txt``               Log filename or '-' for stdout
        ===================== =================================================

        Returns a fully initialized Basis object.
        """

        if txt == '-':
            txt = sys.stdout
        elif txt is None:
            txt = devnull

        if isinstance(tailnorm, float):
            tailnorm = (tailnorm,)
        if 1 + len(tailnorm) < max(polarizationcount, zetacount):
            raise ValueError(
                'Needs %d tail norm values, but only %d are specified' %
                (max(polarizationcount, zetacount) - 1, len(tailnorm)))

        textbuffer = StringIO()
        
        class TeeStream:  # quick hack to both write and save output
            def __init__(self, out1, out2):
                self.out1 = out1
                self.out2 = out2
                
            def write(self, string):
                self.out1.write(string)
                self.out2.write(string)
                
        txt = TeeStream(txt, textbuffer)

        if vconf_args is not None:
            amplitude, ri_rel = vconf_args

        g = self.generator
        rgd = self.rgd

        njcore = g.njcore
        n_j = g.n_j[njcore:]
        l_j = g.l_j[njcore:]
        f_j = g.f_j[njcore:]

        if jvalues is None:
            jvalues = []
            sortkeys = []
            for j in range(len(n_j)):
                if f_j[j] == 0 and l_j[j] != 0:
                    continue
                jvalues.append(j)
                sortkeys.append(l_j[j])
            
            # Now order jvalues by l
            #
            # Use a stable sort so the energy ordering within each
            # angular momentum is guaranteed to be preserved
            args = np.argsort(sortkeys, kind='mergesort')
            jvalues = np.array(jvalues)[args]

        fulljvalues = [njcore + j for j in jvalues]
        
        if isinstance(energysplit, float):
            energysplit = [energysplit] * len(jvalues)
        
        title = '%s Basis functions for %s' % (g.xcname, g.symbol)
        print(title, file=txt)
        print('=' * len(title), file=txt)
        
        singlezetas = []
        energy_derivative_functions = []
        multizetas = [[] for i in range(zetacount - 1)]
        polarization_functions = []

        splitvalencedescr = 'split-valence wave, fixed tail norm'
        derivativedescr = 'derivative of sz wrt. (ri/rc) of potential'

        for vj, fullj, esplit in zip(jvalues, fulljvalues, energysplit):
            l = l_j[vj]
            n = n_j[vj]
            assert n > 0
            orbitaltype = str(n) + 'spdf'[l]
            msg = 'Basis functions for l=%d, n=%d' % (l, n)
            print(file=txt)
            print(msg + '\n', '-' * len(msg), file=txt)
            print(file=txt)
            if vconf_args is None:
                adverb = 'sharply'
            else:
                adverb = 'softly'
            print('Zeta 1: %s confined pseudo wave,' % adverb, end=' ', file=txt)

            u, e, de, vconf, rc = self.rcut_by_energy(fullj, esplit,
                                                      tolerance,
                                                      vconf_args=vconf_args)
            if rc > rcutmax:
                rc = rcutmax  # scale things down
                if vconf is not None:
                    vconf = g.get_confinement_potential(amplitude, ri_rel * rc,
                                                        rc)
                u, e = g.solve_confined(fullj, rc, vconf)
                print('using maximum cutoff', file=txt)
                print('rc=%.02f Bohr' % rc, file=txt)
            else:
                print('fixed energy shift', file=txt)
                print('DE=%.03f eV :: rc=%.02f Bohr' % (de * Hartree,
                                                                rc), file=txt)
            if vconf is not None:
                print(('Potential amp=%.02f :: ri/rc=%.02f' %
                               (amplitude, ri_rel)), file=txt)
            phit_g = self.smoothify(u, l)
            bf = BasisFunction(l, rc, phit_g,
                               '%s-sz confined orbital' % orbitaltype)
            norm = np.dot(g.dr, phit_g * phit_g)**.5
            print('Norm=%.03f' % norm, file=txt)
            singlezetas.append(bf)

            zetacounter = iter(xrange(2, zetacount + 1))

            if include_energy_derivatives:
                assert zetacount > 1
                zeta = zetacounter.next()
                print('\nZeta %d: %s' % (zeta, derivativedescr), file=txt)
                vconf2 = g.get_confinement_potential(amplitude,
                                                     ri_rel * rc * .99, rc)
                u2, e2 = g.solve_confined(fullj, rc, vconf2)
                
                phit2_g = self.smoothify(u2, l)
                dphit_g = phit2_g - phit_g
                dphit_norm = np.dot(rgd.dr_g, dphit_g * dphit_g) ** .5
                dphit_g /= dphit_norm
                descr = '%s-dz E-derivative of sz' % orbitaltype
                bf = BasisFunction(l, rc, dphit_g, descr)
                energy_derivative_functions.append(bf)

            for i, zeta in enumerate(zetacounter):
                print('\nZeta %d: %s' % (zeta, splitvalencedescr), file=txt)
                # Unresolved issue:  how does the lack of normalization
                # of the first function impact the tail norm scheme?
                # Presumably not much, since most interesting stuff happens
                # close to the core.
                rsplit, norm, splitwave = rsplit_by_norm(rgd, l, phit_g,
                                                         tailnorm[i]**2.0,
                                                         txt)
                descr = '%s-%sz split-valence wave' % (orbitaltype,
                                                       '0sdtq56789'[zeta])
                bf = BasisFunction(l, rsplit, phit_g - splitwave, descr)
                multizetas[i].append(bf)
            
        if polarizationcount > 0 or l_pol is not None:
            if l_pol is None:
                # Now make up some properties for the polarization orbital
                # We just use the cutoffs from the previous one times a factor
                # Find 'missing' values in lvalues
                lvalues = [l_j[vj] for vj in jvalues]
                for i in range(max(lvalues) + 1):
                    if list(lvalues).count(i) == 0:
                        l_pol = i
                        break
                else:
                        l_pol = max(lvalues) + 1

            # Find the last state with l=l_pol - 1, which will be the state we
            # base the polarization function on
            for vj, fullj, bf in zip(jvalues[::-1], fulljvalues[::-1],
                                     singlezetas[::-1]):
                if bf.l == l_pol - 1:
                    fullj_pol = fullj
                    rcut = bf.rc * rcutpol_rel
                    break
            else:
                raise ValueError('The requested value l_pol=%d requires l=%d '
                                 'among valence states' % (l_pol, l_pol - 1))
            rcut = min(rcut, rcutmax)
            msg = 'Polarization function: l=%d, rc=%.02f' % (l_pol, rcut)
            print('\n' + msg, file=txt)
            print('-' * len(msg), file=txt)
            # Make a single Gaussian for polarization function.
            #
            # It is known that for given l, the sz cutoff defined
            # by some fixed energy is strongly correlated to the
            # value of the characteristic radius which best reproduces
            # the wave function found by interpolation.
            #
            # We know that for e.g. d orbitals:
            #   rchar ~= .37 rcut[sz](.3eV)
            # Since we don't want to spend a lot of time finding
            # these value for other energies, we just find the energy
            # shift at .3 eV now

            u, e, de, vconf, rc_fixed = self.rcut_by_energy(fullj_pol,
                                                            .3, 1e-2,
                                                            6., (12., .6))

            default_rchar_rel = .25
            # Defaults for each l.  Actually we don't care right now
            rchar_rels = {}

            if rcharpol_rel is None:
                rcharpol_rel = rchar_rels.get(l_pol, default_rchar_rel)
            rchar = rcharpol_rel * rc_fixed
            gaussian = QuasiGaussian(1.0 / rchar**2, rcut)
            psi_pol = gaussian(rgd.r_g) * rgd.r_g**(l_pol + 1)
            norm = np.dot(rgd.dr_g, psi_pol * psi_pol) ** .5
            psi_pol /= norm
            print('Single quasi Gaussian', file=txt)
            msg = 'Rchar = %.03f*rcut = %.03f Bohr' % (rcharpol_rel, rchar)
            adjective = 'Gaussian'
            print(msg, file=txt)
            type = '%s-type %s polarization' % ('spdfg'[l_pol], adjective)
            bf_pol = BasisFunction(l_pol, rcut, psi_pol, type)
                                   
            polarization_functions.append(bf_pol)
            for i in range(polarizationcount - 1):
                npol = i + 2
                msg = '\n%s: %s' % (['Secondary', 'Tertiary', 'Quaternary',
                                     'Quintary', 'Sextary', 'Septenary'][i],
                                    splitvalencedescr)
                print(msg, file=txt)
                rsplit, norm, splitwave = rsplit_by_norm(rgd, l_pol, psi_pol,
                                                         tailnorm[i],
                                                         txt)
                descr = ('%s-type split-valence polarization %d'
                         % ('spdfg'[l_pol], npol))
                bf_pol = BasisFunction(l_pol, rsplit, psi_pol - splitwave,
                                       descr)
                polarization_functions.append(bf_pol)
        
        bf_j = []
        bf_j.extend(singlezetas)
        bf_j.extend(energy_derivative_functions)
        for multizeta_list in multizetas:
            bf_j.extend(multizeta_list)
        bf_j.extend(polarization_functions)
        
        rcmax = max([bf.rc for bf in bf_j])

        # The non-equidistant grids are really only suited for AE WFs
        d = 1.0 / 64
        equidistant_grid = np.arange(0.0, rcmax + d, d)
        ng = len(equidistant_grid)

        for bf in bf_j:
            # We have been storing phit_g * r, but we just want phit_g
            bf.phit_g = divrl(bf.phit_g, 1, rgd.r_g)
            
            gcut = min(int(1 + bf.rc / d), ng - 1)
            
            assert equidistant_grid[gcut] >= bf.rc
            assert equidistant_grid[gcut - 1] <= bf.rc
            
            bf.rc = equidistant_grid[gcut]
            # Note: bf.rc *must* correspond to a grid point (spline issues)
            bf.ng = gcut + 1
            # XXX all this should be done while building the basis vectors,
            # not here
            
            # Quick hack to change to equidistant coordinates
            spline = rgd.spline(bf.phit_g, rgd.r_g[rgd.floor(bf.rc)], bf.l,
                                points=100)
            bf.phit_g = np.array([spline(r) * r**bf.l
                                  for r in equidistant_grid[:bf.ng]])
            bf.phit_g[-1] = 0.

        basis = Basis(g.symbol, self.name, False)
        basis.ng = ng
        basis.d = d
        basis.bf_j = bf_j
        basis.generatordata = textbuffer.getvalue().strip()
        basis.generatorattrs = {'version': version}
        textbuffer.close()

        return basis

    def grplot(self, bf_j):
        """Plot basis functions on generator's radial grid."""
        import pylab
        rc = max([bf.rc for bf in bf_j])
        g = self.generator
        r = g.r
        for bf in bf_j:
            label = bf.type
            # XXX times g.r or not times g.r ?
            pylab.plot(r, bf.phit_g / r, label=label[:12])
        axis = pylab.axis()
        newaxis = [0., rc, axis[2], axis[3]]
        pylab.axis(newaxis)
        pylab.legend()
        pylab.show()

    def plot(self, basis, figure=None, title=None, filename=None):
        """Plot basis functions using pylab."""
        # XXX method should no longer belong to a basis maker
        import pylab
        rc = max([bf.rc for bf in basis.bf_j])
        r = np.linspace(0., basis.d * (basis.ng - 1), basis.ng)
        g = self.generator
        if figure is not None:
            pylab.figure(figure)
        else:
            pylab.figure()  # not very elegant
        if title is None:
            title = g.symbol
        pylab.title(title)
        for bf in basis.bf_j:
            label = bf.type
            # XXX times g.r or not times g.r ?
            phit_g = np.zeros_like(r)
            phit_g[:len(bf.phit_g)] = bf.phit_g
            pylab.plot(r, phit_g * r, label=label[:12])
        axis = pylab.axis()
        newaxis = [0., rc, axis[2], axis[3]]
        pylab.axis(newaxis)
        pylab.legend()
        if filename is not None:
            pylab.savefig(filename)


class QuasiGaussian:
    """Gaussian-like functions for expansion of orbitals.

    Implements f(r) = A [G(r) - P(r)] where::

      G(r) = exp{- alpha r^2}
      P(r) = a - b r^2

    with (a, b) such that f(rcut) == f'(rcut) == 0.
    """
    def __init__(self, alpha, rcut, A=1.):
        self.alpha = alpha
        self.rcut = rcut
        expmar2 = np.exp(-alpha * rcut**2)
        a = (1 + alpha * rcut**2) * expmar2
        b = alpha * expmar2
        self.a = a
        self.b = b
        self.A = A
        
    def __call__(self, r):
        """Evaluate function values at r, which is a numpy array."""
        condition = (r < self.rcut) & (self.alpha * r**2 < 700.)
        r2 = np.where(condition, r**2., 0.)  # prevent overflow
        g = np.exp(-self.alpha * r2)
        p = (self.a - self.b * r2)
        y = np.where(condition, g - p, 0.)
        return self.A * y

    def renormalize(self, norm):
        """Divide function by norm."""
        self.A /= norm
