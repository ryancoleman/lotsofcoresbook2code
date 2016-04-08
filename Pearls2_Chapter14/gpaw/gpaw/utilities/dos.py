from __future__ import print_function
from math import pi, sqrt
import numpy as np
from ase.units import Hartree
from ase.parallel import paropen
from gpaw.utilities import pack
from gpaw.analyse.wignerseitz import wignerseitz
from gpaw.setup_data import SetupData
from gpaw.gauss import Gauss
from gpaw.io.fmf import FMF
from gpaw.utilities.blas import gemmdot
from itertools import izip


def print_projectors(setup):
    """Print information on the projectors of input nucleus object.

    If nucleus is a string, treat this as an element name.
    """
    if type(setup) is str:
        setup = SetupData(setup, 'LDA', 'paw')
        n_j = setup.n_j
        l_j = setup.l_j
    else:
        n_j = setup.n_j
        l_j = setup.l_j
    
    angular = [['1'],
               ['y', 'z', 'x'],
               ['xy', 'yz', '3z^2-r^2', 'xz', 'x^2-y^2'],
               ['3x^2y-y^3', 'xyz', '5yz^2-yr^2', '5z^3-3zr^2',
                '5xz^2-xr^2', 'x^2z-y^2z', 'x^3-3xy^2'],
               ]
    print(' i n l m')
    print('--------')
    i = 0
    for n, l in zip(n_j, l_j):
        for m in range(2 * l + 1):
            if n == -1:
                n = '*'
            print('%2s %s %s_%s' % (i, n, 'spdf'[l], angular[l][m]))
            i += 1

            
def number_of_projectors(setup):
    """Returns the number of the bound state projectors.

    If setup is a string, treat this as an element name.
    """
    if type(setup) is str:
        setup = SetupData(setup, 'LDA', 'paw')
        n_j = setup.n_j
        l_j = setup.l_j
    else:
        n_j = setup.n_j
        l_j = setup.l_j
    
    i = 0
    for n, l in zip(n_j, l_j):
        for m in range(2 * l + 1):
            if n != -1:
                i += 1
    return i

    
def get_angular_projectors(setup, angular, type='bound'):
    """Determine the projector indices which have specified angula
    quantum number.

    angular can be s, p, d, f, or a list of these.
    If type is 'bound', only bound state projectors are considered, otherwise
    all projectors are included.
    """
    # Get the number of relevant j values
    if type == 'bound':
        nj = len([n for n in setup.n_j if n >= 0])
    else:
        nj = len(setup.n_j)
            

    # Choose the relevant projectors
    projectors = []
    i = j = 0
    for j in range(nj):
        m = 2 * setup.l_j[j] + 1
        if 'spdf'[setup.l_j[j]] in angular:
            projectors.extend(range(i, i + m))
        j += 1
        i += m

    return projectors

    
def delta(x, x0, width, mode='Gauss'):
    """Return a gaussian of given width centered at x0."""
    if mode == 'Gauss':
        return np.exp(np.clip(-((x - x0) / width)**2,
                              -100.0, 100.0)) / (sqrt(pi) * width)
    if mode == 'Lorentz':
        return (2 / pi / width) / ((np.clip(((x - x0) / (width / 2))**2,
                                            -100.0, 100.0)) + 1)

        
def fold(energies, weights, npts, width, mode='Gauss'):
    """Take a list of energies and weights, and sum a delta function
    for each."""
    emin = min(energies) - 5 * width
    emax = max(energies) + 5 * width
    e = np.linspace(emin, emax, npts)
    dos_e = np.zeros(npts)
    for e0, w in zip(energies, weights):
        dos_e += w * delta(e, e0, width, mode=mode)
    return e, dos_e

    
def raw_orbital_LDOS(paw, a, spin, angular='spdf'):
    """Return a list of eigenvalues, and their weight on the specified atom.

    angular can be s, p, d, f, or a list of these.
    If angular is None, the raw weight for each projector is returned.

    An integer value for ``angular`` can also be used to specify a specific
    projector function.
    """
    wfs = paw.wfs
    w_k = wfs.kd.weight_k
    nk = len(w_k)
    nb = wfs.bd.nbands

    if a < 0:
        # Allow list-style negative indices; we'll need the positive a for the
        # dictionary lookup later
        a = len(wfs.setups) + a

    setup = wfs.setups[a]
    energies = np.empty(nb * nk)
    weights_xi = np.empty((nb * nk, setup.ni))
    x = 0
    for k, w in enumerate(w_k):
        eps = wfs.collect_eigenvalues(k=k, s=spin)
        print(wfs.world.rank, type(eps))
        if eps is not None:
            energies[x:x + nb] = eps
        u = spin * nk + k
        P_ani = wfs.kpt_u[u].P_ani
        if a in P_ani:
            weights_xi[x:x + nb, :] = w * np.absolute(P_ani[a])**2
        x += nb

    wfs.world.broadcast(energies, 0)
    wfs.world.broadcast(weights_xi, wfs.rank_a[a])

    if angular is None:
        return energies, weights_xi
    elif type(angular) is int:
        return energies, weights_xi[:, angular]
    else:
        projectors = get_angular_projectors(setup, angular, type='bound')
        weights = np.sum(np.take(weights_xi,
                                 indices=projectors, axis=1), axis=1)
        return energies, weights

        
def all_electron_LDOS(paw, mol, spin, lc=None, wf_k=None, P_aui=None):
    """Returns a list of eigenvalues, and their weights on a given molecule
    
       wf_k should be a list of pseudo_wavefunctons of a Kohn-Sham state,
       corresponding to the different kpoints. It should be accompanied by a
       list of projector overlaps: P_aui=[[kpt.P_ani[a][n] for kpt in
       paw.wfs.kpt_u] for a in range(len(molecule))] for the band n. The
       weights are then calculated as the overlap of all-electron
       KS wavefunctions with wf_k

       If wf_k is None, the weights are calculated as linear combinations of
       atomic orbitals using P_uni. lc should then be a list of weights
       for each atom. For example, the pure 2pi_x orbital of a CO or N2
       molecule can be obtained with lc=[[0,0,0,1.0],[0,0,0,-1.0]]. mol
       should be a list of atom numbers contributing to the molecule."""

    w_k = paw.wfs.kd.weight_k
    nk = len(w_k)
    nb = paw.wfs.bd.nbands
    
    P_kn = np.zeros((nk, nb), np.complex)
    if wf_k is None:
        if lc is None:
            lc = [[1,0,0,0] for a in mol]
        for k, kpt in enumerate(paw.wfs.kpt_u[spin * nk:(spin + 1) * nk]):
            N = 0
            for atom, w_a in zip(mol, lc):
                i=0
                for w_o in w_a:
                    P_kn[k] += w_o * kpt.P_ani[atom][:, i]
                    N += abs(w_o)**2
                    i +=1
        P_kn /= sqrt(N)

    else:
        P_aui = [np.array(P_ui).conj() for P_ui in P_aui]
        for k, kpt in enumerate(paw.wfs.kpt_u[spin * nk:(spin + 1) * nk]):
            for n in range(nb):
                P_kn[k][n] = paw.wfs.integrate(wf_k[k], kpt.psit_nG[n])
                for a, b in zip(mol, range(len(mol))):
                    atom = paw.wfs.setups[a]
                    p_i = kpt.P_ani[a][n]
                    for i in range(len(p_i)):
                        for j in range(len(p_i)):
                            P_kn[k][n] += (P_aui[b][spin*nk + k][i] *
                                           atom.dO_ii[i][j] * p_i[j])
            print('# k', k, ' Sum_m |<m|n>|^2 =',  sum(abs(P_kn[k])**2))
          
    energies = np.empty(nb * nk)
    weights = np.empty(nb * nk)
    x = 0
    for k, w in enumerate(w_k):
        energies[x:x + nb] = paw.wfs.collect_eigenvalues(k=k, s=spin)
        weights[x:x + nb] = w * abs(P_kn[k])**2
        x += nb
        
    return energies, weights

    
def get_all_electron_IPR(paw):
    density = paw.density
    wfs = paw.wfs
    n_G = wfs.gd.empty()
    n_g = density.finegd.empty()
    print()
    print("inverse participation function")
    print("-"*35)
    print("%5s %5s %10s %10s" % ("k","band","eps","ipr"))
    print("-"*35)
    for k, kpt in enumerate(paw.wfs.kpt_u):
        for n, (eps, psit_G)  in enumerate(zip(kpt.eps_n, kpt.psit_nG)):
            n_G[:] = 0.0
            wfs.add_orbital_density(n_G, kpt, n)
            density.interpolator.apply(n_G, n_g)
            norm = density.finegd.integrate(n_g)
            n_g = n_g ** 2
            ipr = density.finegd.integrate(n_g)
            for a in kpt.P_ani:
                # Get xccorr for atom a
                setup = paw.density.setups[a]
                xccorr = setup.xc_correction

                # Get D_sp for atom a
                D_sp = np.array(wfs.get_orbital_density_matrix(a, kpt, n))

                # density a function of L and partial wave radial pair density coefficient
                D_sLq = gemmdot(D_sp, xccorr.B_Lqp, trans='t')
                
                # Create pseudo/ae density iterators for integration
                n_iter = xccorr.expand_density(D_sLq, xccorr.n_qg, None)
                nt_iter = xccorr.expand_density(D_sLq, xccorr.nt_qg, None)

                # Take the spherical average of smooth and ae radial xc potentials
                for n_sg, nt_sg, integrator in izip(n_iter,
                                                    nt_iter,
                                                    xccorr.get_integrator(None)):
                    ipr += integrator.weight * np.sum((n_sg[0]**2-nt_sg[0]**2) * xccorr.rgd.dv_g)
                    norm += integrator.weight * np.sum((n_sg[0]-nt_sg[0]) * xccorr.rgd.dv_g)

            print("%5i %5i %10.5f %10.5f" % (k, n, eps, ipr/norm**2))
    print("-"*35)
            
                            
def raw_wignerseitz_LDOS(paw, a, spin):
    """Return a list of eigenvalues, and their weight on the specified atom"""
    wfs = paw.wfs
    assert wfs.dtype == float
    gd = wfs.gd
    atom_index = wignerseitz(gd, paw.atoms)

    w_k = wfs.kd.weight_k
    nk = len(w_k)
    nb = wfs.bd.nbands

    energies = np.empty(nb * nk)
    weights = np.empty(nb * nk)
    x = 0
    for k, w in enumerate(w_k):
        u = spin * nk + k
        energies[x:x + nb] = wfs.collect_eigenvalues(k=k, s=spin)
        for n, psit_G in enumerate(wfs.kpt_u[u].psit_nG):
            P_i = wfs.kpt_u[u].P_ani[a][n]
            P_p = pack(np.outer(P_i, P_i))
            Delta_p = sqrt(4 * pi) * wfs.setups[a].Delta_pL[:, 0]
            weights[x + n] = w * (gd.integrate(abs(
                np.where(atom_index == a, psit_G, 0.0))**2)
                                  + np.dot(Delta_p, P_p))
        x += nb
    return energies, weights


class RawLDOS:
    """Class to get the unfolded LDOS"""
    def __init__(self, calc):
        self.paw = calc
        for setup in calc.wfs.setups.setups.values():
            if not hasattr(setup, 'l_i'):
                # get the mapping
                l_i = []
                for l in setup.l_j:
                    for m in range(2 * l + 1):
                        l_i.append(l)
                setup.l_i = l_i

    def get(self, atom):
        """Return the s,p,d weights for each state"""
        wfs = self.paw.wfs
        nibzkpts = len(wfs.kd.ibzk_kc)
        spd = np.zeros((wfs.nspins, nibzkpts, wfs.bd.nbands, 3))

        if hasattr(atom, '__iter__'):
            # atom is a list of atom indicies 
            for a in atom:
                spd += self.get(a)
            return spd

        l_i = wfs.setups[atom].l_i

        for kpt in self.paw.wfs.kpt_u:
            if atom in kpt.P_ani:
                for i, P_n in enumerate(kpt.P_ani[atom].T):
                    spd[kpt.s, kpt.k, :, l_i[i]] += np.abs(P_n)**2

        wfs.gd.comm.sum(spd)
        wfs.kd.comm.sum(spd)
        return spd

    def by_element(self):
        """Return a dict with elements as keys and LDOS as values."""
        elemi = {}
        for i,a in enumerate(self.paw.atoms):
            symbol = a.symbol
            if symbol in elemi:
                elemi[symbol].append(i)
            else:
                elemi[symbol] = [i]
        for key in elemi:
            elemi[key] = self.get(elemi[key])
        return elemi

    def to_file(self, 
                ldbe,
                filename=None,
                width=None,
                shift=True,
                bound=False):
        """Write the LDOS to a file.

        If a width is given, the LDOS will be Gaussian folded and shifted to set 
        Fermi energy to 0 eV. The latter can be avoided by setting shift=False. 

        If you use fixmagmom=true, you will get two fermi-levels, one for each 
        spin-setting. Normaly these will shifted individually to 0 eV. If you
        want to shift them as pair to the higher energy use bound=True.
        """

        f = paropen(filename, 'w')

        def append_weight_strings(ldbe, data):
            s = ''
            for key in ldbe:
                for l in 'spd':
                    data.append(key + '(' + l + ')-weight')
                if len(key) == 1: 
                    key = ' ' + key
                s +=  ' ' + key + ':s     p        d       '
            return s

        wfs = self.paw.wfs
        
        if width is None:
            # unfolded ldos
            fmf = FMF(['Raw LDOS obtained from projector weights'])
            print(fmf.header(), end=' ', file=f)
            data = ['energy: energy [eV]',
                    'occupation number: occ',
                    'spin index: s',
                    'k-point index: k',
                    'band index: n',
                    'k-point weight: weight']
            string = '# e_i[eV]  occ     s     k      n   kptwght '
            string += append_weight_strings(ldbe, data)
            data.append('summed weights: sum')
            string += ' sum'
            print(fmf.data(data), end=' ', file=f)
            print(string, file=f)
            for k in range(wfs.kd.nibzkpts):
                for s in range(wfs.nspins):
                    e_n = self.paw.get_eigenvalues(kpt=k, spin=s)
                    f_n = self.paw.get_occupation_numbers(kpt=k, spin=s)
                    if e_n is None:
                        continue
                    w = wfs.kd.weight_k[k]
                    for n in range(wfs.bd.nbands):
                        sum = 0.0
                        print('%10.5f %6.4f %2d %5d' % (e_n[n], f_n[n], 
                                                              s, k), end=' ', file=f) 
                        print('%6d %8.4f' % (n, w), end=' ', file=f)
                        for key in ldbe:
                            spd = ldbe[key][s, k, n]
                            for l in range(3):
                                sum += spd[l]
                                print('%8.4f' % spd[l], end=' ', file=f)
                        print('%8.4f' % sum, file=f)
        else:
            # folded ldos
            fmf = FMF(['Folded raw LDOS obtained from projector weights'])
            print(fmf.header(), end=' ', file=f)

            gauss = Gauss(width)
            print(fmf.field('folding',
                                  ['function: Gauss',
                                   'width: ' + str(width) + ' [eV]']), end=' ', file=f)

            data = ['energy: energy [eV]',
                    'spin index: s',
                    'k-point index: k',
                    'band index: n',
                    'k-point weight: weight']

            # minimal and maximal energies
            emin = 1.e32
            emax = -1.e32
            for k in range(wfs.kd.nibzkpts):
                for s in range(wfs.nspins):
                    e_n = self.paw.get_eigenvalues(kpt=k, spin=s,
                                                   broadcast=True)
                    emin = min(e_n.min(), emin)
                    emax = max(e_n.max(), emax)
            emin -= 4 * width
            emax += 4 * width

            # Fermi energy
            try:
                if self.paw.occupations.fixmagmom:
                    efermi = self.paw.get_fermi_levels()
                else:
                    efermi = self.paw.get_fermi_level()
            except:
                # set Fermi level half way between HOMO and LUMO
                hl = self.paw.occupations.get_homo_lumo(wfs)
                efermi = (hl[0] + hl[1]) * Hartree / 2

            eshift = 0.0

            if shift and not self.paw.occupations.fixmagmom:
                eshift = -efermi

            # set de to sample 4 points in the width
            de = width / 4
            
            string = '# e[eV]     s  '
            string += append_weight_strings(ldbe, data)

            for s in range(wfs.nspins):
                if self.paw.occupations.fixmagmom:
                    if not bound:
                        eshift = - efermi[s]
                    else:
                        eshift = - efermi.max()

                print(fmf.data(data), end=' ', file=f)

                print('# Gauss folded, width=%g [eV]' % width, file=f)
                if shift:
                    print('# shifted to Fermi energy = 0', file=f)
                    print('# Fermi energy was', end=' ', file=f) 
                else:
                    print('# Fermi energy', end=' ', file=f)
                print(efermi, 'eV', file=f)
                print(string, file=f)

                # loop over energies
                emax=emax+.5*de
                e=emin
                while e<emax:
                    val = {}
                    for key in ldbe:
                        val[key] = np.zeros((3))
                    for k in range(wfs.kd.nibzkpts):
                        w = wfs.kpt_u[k].weight
                        e_n = self.paw.get_eigenvalues(kpt=k, spin=s,
                                                       broadcast=True)
                        for n in range(wfs.bd.nbands):
                            w_i = w * gauss.get(e_n[n] - e)
                            for key in ldbe:
                                val[key] += w_i * ldbe[key][s, k, n]

                    print('%10.5f %2d' % (e + eshift, s), end=' ', file=f) 
                    for key in val:
                        spd = val[key]
                        for l in range(3):
                            print('%8.4f' % spd[l], end=' ', file=f)
                    print(file=f)
                    e += de
                            

        f.close()

    def by_element_to_file(self, 
                           filename='ldos_by_element.dat',
                           width=None,
                           shift=True,
                           bound=False):
        """Write the LDOS by element to a file.
        """
        ldbe = self.by_element()
        self.to_file(ldbe, filename, width, shift, bound)


class LCAODOS:
    """Class for calculating the projected subspace DOS.

    The projected subspace density of states is defined only in LCAO
    mode.  The advantages to the PDOS based on projectors is that the
    LCAODOS will properly take non-orthogonality and completeness into
    account."""
    def __init__(self, calc):
        self.calc = calc

    def get_orbital_pdos(self, M, ravel=True):
        """Get projected DOS from LCAO basis function M."""
        return self.get_subspace_pdos([M], ravel=ravel)
    
    def get_atomic_subspace_pdos(self, a, ravel=True):
        """Get projected subspace DOS from LCAO basis on atom a."""
        M = self.calc.wfs.basis_functions.M_a[a]
        Mvalues = range(M, M + self.calc.wfs.setups[a].nao)
        return self.get_subspace_pdos(Mvalues, ravel=ravel)

    def get_subspace_pdos(self, Mvalues, ravel=True):
        """Get projected subspace DOS from LCAO basis."""
        wfs = self.calc.wfs
        bd = wfs.bd
        kd = wfs.kd
        gd = wfs.gd
        
        for kpt in wfs.kpt_u:
            assert not np.isnan(kpt.eps_n).any()

        w_skn = np.zeros((kd.nspins, kd.nks, bd.nbands))
        eps_skn = np.zeros((kd.nspins, kd.nks, bd.nbands))
        for u, kpt in enumerate(wfs.kpt_u):
            C_nM = kpt.C_nM
            from gpaw.kohnsham_layouts import BlacsOrbitalLayouts
            if isinstance(wfs.ksl, BlacsOrbitalLayouts):
                S_MM = wfs.ksl.mmdescriptor.collect_on_master(kpt.S_MM)
                if bd.rank != 0 or gd.rank != 0:
                    S_MM = np.empty((wfs.ksl.nao, wfs.ksl.nao),
                                    dtype=wfs.dtype)
                bd.comm.broadcast(S_MM, 0)
            else:
                S_MM = kpt.S_MM
            
            if gd.comm.rank == 0:
                S_mm = S_MM[Mvalues, :][:, Mvalues].copy()
                iS_mm = np.linalg.inv(S_mm).copy()
                CS_nm = np.dot(C_nM, S_MM[:, Mvalues])
                CSiS_nm = np.dot(CS_nm, iS_mm)
                w_n = (np.conj(CS_nm) * CSiS_nm).sum(axis=1) * kpt.weight
                w_skn[kpt.s, kpt.k, bd.beg:bd.end:bd.step] = w_n.real
                eps_skn[kpt.s, kpt.k, bd.beg:bd.end:bd.step] = kpt.eps_n
        
        for arr in [eps_skn, w_skn]:
            gd.comm.broadcast(arr, 0)
            bd.comm.sum(arr)
            kd.comm.sum(arr)

        if ravel:
            eps_n = eps_skn.ravel()
            w_n = w_skn.ravel()
            args = np.argsort(eps_n)
            eps_n = eps_n[args]
            w_n = w_n[args]
            return eps_n, w_n
        else:
            return eps_skn, w_skn


class RestartLCAODOS(LCAODOS):
    """Class for calculating LCAO subspace PDOS from a restarted calculator.

    Warning: This has side effects on the calculator.  The
    operation will allocate memory to diagonalize the Hamiltonian and
    set coefficients plus positions."""
    def __init__(self, calc):
        LCAODOS.__init__(self, calc)
        system = calc.get_atoms()
        calc.set_positions(system)
        calc.wfs.eigensolver.iterate(calc.hamiltonian, calc.wfs)
