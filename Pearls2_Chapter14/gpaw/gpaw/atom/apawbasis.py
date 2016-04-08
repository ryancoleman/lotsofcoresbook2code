import sys
from optparse import OptionParser

import numpy as np

from gpaw.atom.atompaw import AtomPAW
from gpaw.atom.basis import rsplit_by_norm, QuasiGaussian,\
     get_gaussianlike_basis_function
from gpaw.basis_data import BasisFunction, Basis
from gpaw.hgh import setups as hgh_setups, sc_setups as hgh_sc_setups,\
     HGHSetupData

# XXX
from scipy.optimize import bisect

def atompaw(setup, f_ln, rcut, **kwargs):
    return AtomPAW(setup.symbol,
                   [f_ln],
                   rcut=rcut,
                   setups={setup.symbol : setup},
                   **kwargs)


def get_orbitals_by_energy_shift(opts, setup, **kwargs):
    h = opts.grid
    try:
        f_ln = setup.f_ln
    except AttributeError:
        f_ln = []
        for n, l, f in zip(setup.n_j, setup.l_j, setup.f_j):
            if n < 0:
                continue
            if l == len(f_ln):
                f_ln.append([])
            if f > 0 and n > 0:
                f_ln[l].append(f)
    def calculate(rcut, h=h, txt='-'):
        return atompaw(setup, f_ln, rcut, h=h, txt=txt, **kwargs)

    def get_orbital_energy(l0, n0, rcut):
        calc = calculate(rcut, 0.15, txt=None)
        for l, n, f, eps, psit_G in calc.state_iter():
            if l == l0 and n + 1== n0: # XXX
                return eps
        raise ValueError('No such valence state: l=%d, n=%d' % (l0, n0))

    calc0 = calculate(2.0, 0.2, txt=None) # XXX
    def valence_states():
        for l, n, f, eps, psit_G in calc0.state_iter():
            yield l, n + 1 # XXX

    bf_j = []
    cutoffs = []
    
    for i, (l, n) in enumerate(valence_states()):
        e0 = get_orbital_energy(l, n, 15.0) * 27.211 # 15Ang == infinity
        print('e0', e0)
        def obj(rcut):
            eps = get_orbital_energy(l, n, rcut) * 27.211
            de = eps - opts.energy_shift - e0
            #print rcut, eps, de
            return de
        
        # Not exactly efficient...
        rcut = bisect(obj, 1.0, 15.0, xtol=0.1)
        calc = calculate(h=h, rcut=rcut, txt=None)
        bfs = calc.extract_basis_functions()
        bf_j.append(bfs.bf_j[i])

    basis = Basis(setup.symbol, 'strange', readxml=False)
    basis.d = bfs.d
    basis.ng = max([bf.ng for bf in bf_j])
    basis.bf_j = bf_j
    return basis
    #for (l, n), cutoff in zip(valence_states(), cutoffs):
    #    calculate(
    
    #return
    #calc = calculate(rcut)
    bfs = calc.extract_basis_functions(basis_name=opts.name)
    ldict = dict([(bf.l, bf) for bf in bfs.bf_j])

    rgd = bfs.get_grid_descriptor()

    def get_rsplit(bf, splitnorm):
        if opts.s_approaches_zero and bf.l == 0:
            l = 1 # create a function phi(r) = A * r + O(r^2)
        else:
            l = bf.l
        return rsplit_by_norm(rgd,
                              l,
                              bf.phit_g * rgd.r_g,
                              splitnorm**2,
                              sys.stdout)

    splitvalence_bfs = []
    for splitnorm in opts.splitnorm:
        splitnorm = float(splitnorm)
        for orbital_bf in bfs.bf_j:
            rsplit, normsqr, phit_g = get_rsplit(orbital_bf, splitnorm)
            phit_g[1:] /= rgd.r_g[1:]
            gcut = rgd.ceil(rsplit)
            #tailnorm = np.dot(rgd.dr_g[gcut:],
            #                  (rgd.r_g[gcut:] * orbital_bf.phit_g[gcut:])**2)**0.5
            #print 'tailnorm', tailnorm
            dphit_g = orbital_bf.phit_g[:gcut+1] - phit_g[:gcut+1]
            bf = BasisFunction(l=orbital_bf.l,
                               rc=rgd.r_g[gcut],
                               phit_g=dphit_g,
                               type='%s split-valence' % 'spd'[orbital_bf.l])
            splitvalence_bfs.append(bf)
    bfs.bf_j.extend(splitvalence_bfs)

    #rpol = None
    for l in range(3):
        if not l in ldict:
            lpol = l
            source_bf = ldict[lpol - 1]
            break
    else:
        raise NotImplementedError('f-type polarization not implemented')
    
    for splitnorm in opts.polarization:
        splitnorm = float(splitnorm)
        rchar, normsqr, phit_g = get_rsplit(source_bf, splitnorm)
        gcut = rgd.ceil(3.5 * rchar)
        rcut = rgd.r_g[gcut]
        phit_g = get_gaussianlike_basis_function(rgd, lpol, rchar, gcut)
        N = len(phit_g)
        x = np.dot(rgd.dr_g[:N], (phit_g * rgd.r_g[:N])**2)**0.5
        print('x', x)
        bf = BasisFunction(lpol,
                           rc=rcut,
                           phit_g=phit_g,
                           type='%s polarization' % 'spd'[lpol])
        bf.phit_g
        bfs.bf_j.append(bf)
    bfs.write_xml()


def build_parser():
    usage = '%prog [OPTION] [SYMBOL...]'
    
    description = 'generate basis sets from existing setups.'
    
    p = OptionParser(usage=usage, description=description)
    p.add_option('--grid', default=0.05, type=float, metavar='DR',
                 help='grid spacing in atomic calculation')
    #p.add_option('--rcut', default=10.0, type=float,
    #             help='radial cutoff for atomic calculation')
    p.add_option('-f', '--file', action='store_true',
                 help='load setups directly from files')
    p.add_option('-E', '--energy-shift', type=float, default=0.1,
                 help='use given energy shift to determine cutoff')
    p.add_option('-n', '--name', default='apaw',
                 help='name of basis set')
    #p.add_option('--splitnorm', action='append', default=[], metavar='NORM',
    #             help='add split-valence basis functions at this'
    #             ' tail norm.  Multiple options accumulate')
    #p.add_option('--polarization', action='append', default=[],
    #             metavar='NORM',
    #             help='add polarization function with characteristic radius'
    #             ' determined by tail norm.  Multiple options accumulate')
    #p.add_option('--s-approaches-zero', action='store_true',
    #             help='force s-type split-valence functions to 0 at origin')
    #p.add_option('-t', '--type',
    #             help='string describing extra basis functions')
                 
    return p


def main():
    p = build_parser()
    opts, args = p.parse_args()

    if opts.file:
        for fname in args:
            tokens = fname.split('.')
            symbol = tokens[0]
            xc = tokens[1]
            name = tokens[2]
            if tokens[-1] == 'gz':
                import gzip
                fopen = gzip.open
            else:
                fopen = open
            source = fopen(fname).read()

            from gpaw.setup_data import SetupData

            s = SetupData(symbol, xc, name, readxml=False)
            s.read_xml(source=source)
            basis = get_orbitals_by_energy_shift(opts, s, xc=xc)
            basis.write_xml()
    else:
        for arg in args:
            setup = hgh_setups.get(arg)
            if setup is None:
                setup = hgh_sc_setups.get(arg.split('.')[0])
            if setup is None:
                raise ValueError('Unknown setup %s' % arg)
            print(setup)
            basis = get_orbitals_by_energy_shift(opts, HGHSetupData(setup))
            basis.write_xml()
            #generate_basis(opts, HGHSetupData(setup))
