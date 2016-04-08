
import numpy as np

from gpaw import debug, parsize_domain, parsize_bands
from gpaw import mpi
from gpaw.utilities.blas import axpy

"""
def dscf_find_lumo(paw,band):

    # http://trac.fysik.dtu.dk/projects/gpaw/browser/trunk/doc/documentation/dscf/lumo.py?format=raw

    assert band in [5,6]

    #Find band corresponding to lumo
    lumo = paw.get_pseudo_wave_function(band=band, kpt=0, spin=0)
    lumo = np.reshape(lumo, -1)

    wf1_k = [paw.get_pseudo_wave_function(band=5, kpt=k, spin=0) for k in range(paw.wfs.kd.nibzkpts)]
    wf2_k = [paw.get_pseudo_wave_function(band=6, kpt=k, spin=0) for k in range(paw.wfs.kd.nibzkpts)]

    band_k = []
    for k in range(paw.wfs.kd.nibzkpts):
        wf1 = np.reshape(wf1_k[k], -1)
        wf2 = np.reshape(wf2_k[k], -1)
        p1 = np.abs(np.dot(wf1, lumo))
        p2 = np.abs(np.dot(wf2, lumo))

        if p1 > p2:
            band_k.append(5)
        else:
            band_k.append(6)

    return band_k
"""

# -------------------------------------------------------------------

"""
def mpi_debug(text):
    if isinstance(text,list):
        for t in text:
            mpi_debug(t)
    else:
        print 'mpi.rank=%d, %s' % (mpi.rank,text)
"""

global msgcount
msgcount = 0

def mpi_debug(data, ordered=True):
    global msgcount

    if not isinstance(data, list):
        data = [data]

    if ordered:
        for i in range(mpi.rank):
            mpi.world.barrier()

    for txt in data:
        print('%02d-mpi%d, %s' % (msgcount, mpi.rank, txt))
        if ordered:
            msgcount += 1

    if ordered:
        for i in range(mpi.size-mpi.rank):
            mpi.world.barrier()

# -------------------------------------------------------------------

def dscf_find_atoms(atoms,symbol):
    chemsyms = atoms.get_chemical_symbols()
    return np.where(map(lambda s: s==symbol,chemsyms))[0]

# -------------------------------------------------------------------

# Helper function in case of hs_operators having ngroups > 1
def SliceGen(psit_nG, operator):
    assert psit_nG.ndim == 4
    assert psit_nG.shape[0] == operator.bd.mynbands
    assert np.all(psit_nG.shape[1:] == operator.gd.n_c)
    assert operator.bd.mynbands % operator.nblocks == 0
    M = operator.bd.mynbands // operator.nblocks
    for j in range(operator.nblocks):
        n1 = j * M
        n2 = n1 + M
        yield psit_nG[n1:n2]
    raise StopIteration

from gpaw.kpoint import GlobalKPoint
from gpaw.hs_operators import MatrixOperator

def dscf_kpoint_overlaps(paw, phasemod=True, broadcast=True):
    bd = paw.wfs.bd
    gd = paw.wfs.gd
    kd = paw.wfs.kd
    operator = MatrixOperator(paw.wfs.orthoksl, hermitian=False)
    atoms = paw.get_atoms()

    # Find the kpoint with lowest kpt.k_c (closest to gamma point)
    k0 = np.argmin(np.sum(paw.wfs.kd.ibzk_kc**2,axis=1)**0.5)

    # Maintain list of a single global reference kpoint for each spin
    kpt0_s = []
    for s0 in range(kd.nspins):
        q0 = k0 - kd.get_offset() % kd.nibzkpts
        kpt0 = GlobalKPoint(None, s0, k0, q0, None)
        kpt0.update(paw.wfs)
        kpt0_s.append(kpt0)

    if phasemod:
        # Scaled grid point positions used for exponential with ibzk_kc
        # cf. wavefunctions.py lines 90-91 rev 4500(ca)
        # phase_cd = np.exp(2j * np.pi * sdisp_cd * ibzk_kc[k, :, np.newaxis])
        r_cG = gd.empty(3)
        for c, r_G in enumerate(r_cG):
            slice_c2G = [np.newaxis, np.newaxis, np.newaxis]
            slice_c2G[c] = slice(None) #this means ':'
            r_G[:] = np.arange(gd.beg_c[c], gd.end_c[c], \
                               dtype=float)[slice_c2G] / gd.N_c[c]

    X_unn = np.empty((kd.mynks, bd.nbands, bd.nbands), dtype=paw.wfs.dtype)
    for myu, kpt in enumerate(paw.wfs.kpt_u):
        u = kd.global_index(myu)
        s, k = kd.what_is(u)
        kpt0 = kpt0_s[s]
        X_nn = X_unn[myu]

        if phasemod:
            assert paw.wfs.dtype == complex, 'Phase modification is complex!'

            k0_c = kd.ibzk_kc[k0]
            k_c = kd.ibzk_kc[k]
            eirk_G = np.exp(2j*np.pi*np.sum(r_cG*(k_c-k0_c)[:,np.newaxis,np.newaxis,np.newaxis], axis=0))
            psit0_nG = eirk_G[np.newaxis,...]*kpt0.psit_nG

            P0_ani = paw.wfs.pt.dict(bd.mynbands)
            spos_ac = atoms.get_scaled_positions() % 1.0
            for a, P0_ni in P0_ani.items():
                # Expanding the exponential exp(ikr)=exp(ikR)*exp(ik(r-R))
                # and neglecting the changed P_ani integral exp(ik(r-R))~1
                P0_ni[:] = np.exp(2j*np.pi*np.sum(spos_ac[a]*(k_c-k0_c), axis=0)) * kpt0.P_ani[a]

            ## NB: No exp(ikr) approximate here, but has a parallelization bug
            #kpt0_rank, myu0 = kd.get_rank_and_index(kpt0.s, kpt0.k)
            #if kd.comm.rank == kpt0_rank:
            #    paw.wfs.pt.integrate(psit0_nG, P0_ani, kpt0.q)
            #for a, P0_ni in P0_ani.items():
            #    kd.comm.broadcast(P0_ni, kpt0_rank)
        else:
            psit0_nG = kpt0.psit_nG
            P0_ani = kpt0.P_ani

        """
        if paw.wfs.world.size == 1:
            for n, psit_G in enumerate(kpt.psit_nG):
                for n0, psit0_G in enumerate(psit0_nG):
                    X_nn[n,n0] = np.vdot(psit_G, psit0_G)*gd.dv
            for a in range(len(paw.get_atoms())):
                P_ni, P0_ni, dO_ii = kpt.P_ani[a], P0_ani[a], paw.wfs.setups[a].dO_ii
                for n, P_i in enumerate(P_ni):
                    for n0, P0_i in enumerate(P0_ni):
                        X_nn[n,n0] += np.vdot(P_i, np.dot(dO_ii, P0_i))
        """
        X = lambda psit_nG, g=SliceGen(psit0_nG, operator): g.next()
        dX = lambda a, P_ni: np.dot(P0_ani[a], paw.wfs.setups[a].dO_ii)
        X_nn[:] = operator.calculate_matrix_elements(kpt.psit_nG, kpt.P_ani, X, dX).T

    if broadcast:
        if bd.comm.rank == 0:
            gd.comm.broadcast(X_unn, 0)
        bd.comm.broadcast(X_unn, 0)

    return kpt0_s, X_unn

# -------------------------------------------------------------------

def dscf_find_bands(paw,bands,data=None):
    """Entirely serial, but works regardless of parallelization. DOES NOT WORK WITH DOMAIN-DECOMPOSITION IN GPAW v0.5.2725 """ #TODO!

    raise DeprecationWarning('About to be replaced with something better.')

    if data is None:
        data = range(len(bands))
    else:
        assert len(data)==len(bands), 'Length mismatch.'

    k0 = 0 #TODO find kpt with lowest kpt.k_c (closest to gamma point)

    gamma_siG = []
    for s in range(paw.wfs.nspins):
        gamma_siG.append([paw.get_pseudo_wave_function(band=n,kpt=k0,spin=s).ravel() for n in bands]) #TODO! paw.get_pseudo fails with domain-decomposition from tar-file

    band_ui = []
    data_ui = []

    for u,kpt in enumerate(paw.wfs.kpt_u):
        band_i = []
        data_i = []

        for (i,n) in enumerate(bands):
            if kpt.k == k0:
                wf = gamma_siG[kpt.s][i]
            else:
                #wf = paw.get_pseudo_wave_function(n, kpt.k, kpt.s, pad=False).ravel()
                wf = paw.get_pseudo_wave_function(band=n,kpt=kpt.k,spin=kpt.s).ravel()

            overlaps = [np.abs(np.dot(wf,gamma_siG[kpt.s][i])) for i in range(len(bands))]
            if debug: mpi_debug('u=%d, i=%d, band=%d, overlaps=%s' % (u,i,n,str(overlaps)))
            p = np.argmax(overlaps)
            band_i.append(bands[p])
            data_i.append(data[p])

        assert len(np.unique(band_i))==len(np.sort(band_i)), 'Non-unique band range' #TODO!

        band_ui.append(band_i)
        data_ui.append(data_i)

    return (band_ui,data_ui,)

# -------------------------------------------------------------------

def dscf_linear_combination(paw, molecule, bands, coefficients):
    """Full parallelization over k-point - grid-decomposition parallelization needs heavy testing.""" #TODO!

    raise DeprecationWarning('About to be replaced with something better.')

    if debug: dumpkey = mpi.world.size == 1 and 'serial' or 'mpi'

    (band_ui,coeff_ui,) = dscf_find_bands(paw,bands,coefficients)

    if debug: mpi_debug('band_ui=%s, coeff_ui=%s' % (band_ui,coeff_ui))

    P_aui = {}
    for m,a in enumerate(molecule):
        if debug: mpi_debug('a=%d, paw.wfs.kd.nibzkpts=%d, len(paw.wfs.kpt_u)=%d, paw.wfs.setups[%d].ni=%d' % (a, paw.wfs.kd.nibzkpts, len(paw.wfs.kpt_u), a, paw.wfs.setups[a].ni))
        P_aui[m] = np.zeros((len(paw.wfs.kpt_u),paw.wfs.setups[a].ni),dtype=complex)

    for u,kpt in enumerate(paw.wfs.kpt_u):

        band_i = band_ui[u]
        coeff_i = coeff_ui[u]

        #if debug: mpi_debug(['paw.wfs.kpt_u[%d].P_ani[:,%d,:].shape=%s' % (u,n,str(paw.wfs.kpt_u[u].P_ani[:,n,:].shape)) for n in bands])

        for m,a in enumerate(molecule):
            """
            if debug:
                for n in bands:
                    print 'mpi.rank=%d, paw.nuclei[%d].P_uni[:,%d,:].shape=' % (mpi.rank,a,n), paw.nuclei[a].P_uni[:,n,:].shape

                print 'mpi.rank=%d, test.shape=' % mpi.rank, np.sum([c*paw.nuclei[a].P_uni[:,n,:] for (c,n) in zip(coefficients,bands)],axis=0).shape
            """

            #P_aui[m] += np.sum([c*paw.nuclei[a].P_uni[:,n,:] for (c,n) in zip(coefficients,bands)],axis=0)

            #if paw.nuclei[a].in_this_domain: #TODO what happened to this one in guc?
            if True:
                P_aui[m][u,:] += np.sum([c*kpt.P_ani[a][n,:] for (c,n) in zip(coeff_i,band_i)],axis=0)

                if debug: kpt.P_ani[a][:,:].dump('dscf_tool_P_ani_a%01d_k%01ds%01d_%s%02d.pickle' % (a,kpt.k,kpt.s,dumpkey,mpi.rank))

    #paw.wfs.gd.comm.sum(P_aui) #TODO HUH?!

    if debug: P_aui.dump('dscf_tool_P_aui_%s%02d.pickle' % (dumpkey,mpi.rank))

    """
    if debug and mpi.rank == 0:
        print 'P_aui.shape=',P_aui.shape

        for (a,P_ui) in enumerate(P_aui):
            print 'P_aui[%d].shape=' % a,P_ui.shape

        print 'P_aui=',P_aui

        print 'gd.Nc=',paw.wfs.gd.N_c
    """

    if debug: mpi_debug('P_aui.shape='+str(P_aui.shape))

    #wf_u = [np.sum([c*paw.wfs.kpt_u[u].psit_nG[n] for (c,n) in zip(coefficients,bands)],axis=0) for u in range(0,len(paw.wfs.kpt_u))]
    #wf_u = np.zeros((paw.wfs.kd.nibzkpts,paw.wfs.gd.N_c[0]-1,paw.wfs.gd.N_c[1]-1,paw.wfs.gd.N_c[2]-1))#,dtype=complex)
    wf_u = paw.wfs.gd.zeros(len(paw.wfs.kpt_u),dtype=complex)

    gd_slice = paw.wfs.gd.get_slice()

    if debug: mpi_debug('gd_slice='+str(gd_slice))

    for u,kpt in enumerate(paw.wfs.kpt_u):
        if debug: mpi_debug('u=%d, k=%d, s=%d, paw.wfs.kd.comm.rank=%d, paw.wfs.kd.comm.rank=%d, gd.shape=%s, psit.shape=%s' % (u, kpt.k, kpt.s, paw.wfs.kd.comm.rank, paw.wfs.kd.comm.rank, str(wf_u[0].shape), str(np.array(kpt.psit_nG[0])[gd_slice].shape)))

        #wf_u[u] += np.sum([c*np.array(kpt.psit_nG[n])[gd_slice] for (c,n) in zip(coefficients,bands)],axis=0)

        band_i = band_ui[u]
        coeff_i = coeff_ui[u]
        wf_u[u] += np.sum([c*np.array(kpt.psit_nG[n])[gd_slice] for (c,n) in zip(coeff_i,band_i)],axis=0)

    #paw.wfs.gd.comm.sum(wf_u)

    if debug: mpi_debug('|wf_u|^2=%s' % str([np.sum(np.abs(wf.flatten())**2) for wf in wf_u]))

    """
    if debug and mpi.rank == 0:
        print 'wf_u.shape=',wf_u.shape

        for (u,wf) in enumerate(wf_u):
            print 'wf[%d].shape=' % u,wf.shape

        for (u,wf) in enumerate(wf_u):
            print 'wf_u[%d].shape=' % u,wf.shape
    """

    return (P_aui,wf_u,)

# -------------------------------------------------------------------

def dscf_decompose_occupations(ft_mn):
    # We define a Hermitian matrix with elements f[m,n] = c[m]*c[n].conj()
    # and rewrite the coefficients in polar form as c[n] = r[n]*exp(1j*v[n]):
    # f[m,n] = r[m]*exp(1j*v[m])*r[n]*exp(-1j*v[n]) = r[m]*r[n]*exp(1j*(v[m]-v[n]))

    # Assuming c[i] is positive real for some i, i.e. v[i] = 0, then
    # c[m] = f[m,i]/r[i] = r[m]*exp(1j*v[m]) , where r[i] = sqrt(f[i,i])
    i = np.argmax(np.abs(ft_mn.diagonal()))
    c_n = ft_mn[:,i]/np.abs(ft_mn[i,i])**0.5

    if (np.abs(ft_mn-np.outer(c_n,c_n.conj()))>1e-12).any():
        raise RuntimeError('Hermitian matrix cannot be decomposed')

    # Note that c_n is only defined up to an arbitrary phase factor
    return c_n

# -------------------------------------------------------------------

from gpaw.utilities import unpack
from gpaw.utilities.tools import tri2full

def dscf_matrix_elements(paw, kpt, A, dA):
    operator = paw.wfs.overlap.operator
    A_nn = operator.calculate_matrix_elements(kpt.psit_nG, kpt.P_ani, \
        A, dA).T.copy() # transpose to get A_nn[m,n] = <m|A|n>
    tri2full(A_nn, 'U') # fill in from upper to lower...

    # Calculate <o|A|o'> = sum_nn' <o|n><n|A|n'><n'|o'> where c_on = <n|o>
    A_oo = np.dot(kpt.c_on.conj(), np.dot(A_nn, kpt.c_on.T))
    return A_oo, A_nn

def dscf_overlap_elements(paw, kpt):
    # Copy/paste from gpaw/overlap.py lines 83-86 rev. 4808
    operator = paw.wfs.overlap.operator
    assert operator.nblocks == 1

    S = lambda x: x
    dS_aii = dict([(a, paw.wfs.setups[a].dO_ii) for a in kpt.P_ani.keys()])
    return dscf_matrix_elements(paw, kpt, S, dS_aii)

def dscf_hamiltonian_elements(paw, kpt):
    # Copy/paste from gpaw/eigensolvers/eigensolver.py lines 155-170 rev. 4808
    operator = paw.wfs.overlap.operator
    assert operator.nblocks == 1
    Htpsit_xG = operator.suggest_temporary_buffer(kpt.psit_nG.dtype)

    def H(psit_xG):
        paw.wfs.kin.apply(psit_xG, Htpsit_xG, kpt.phase_cd)
        paw.hamiltonian.apply_local_potential(psit_xG, Htpsit_xG, kpt.s)
        paw.hamiltonian.xc.add_non_local_terms(psit_xG, Htpsit_xG, kpt)
        return Htpsit_xG

    dH_aii = dict([(a, unpack(dH_sp[kpt.s])) for a, dH_sp \
        in paw.hamiltonian.dH_asp.items()])

    H_oo, H_nn = dscf_matrix_elements(paw, kpt, H, dH_aii)
    eps_o = H_oo.real.diagonal()
    return eps_o, H_oo, H_nn

# -------------------------------------------------------------------

from ase.units import Hartree
from gpaw.io.tar import Writer, Reader, FileReference
from gpaw.occupations import FermiDirac
from gpaw.band_descriptor import BandDescriptor

def dscf_reconstruct_orbitals_k_point(paw, norbitals, mol, kpt):
    bd, gd = paw.wfs.bd, paw.wfs.gd
    if bd.comm.size != 1:
        raise NotImplementedError('Undefined action for band parallelization.')

    f_o = np.zeros(norbitals, dtype=float)
    eps_o = np.zeros(norbitals, dtype=float)
    wf_oG = gd.zeros(norbitals, dtype=complex)

    P_aoi = {}
    for a in mol:
        P_aoi[a] = np.zeros((norbitals,paw.wfs.setups[a].ni), dtype=complex)

    for o, c_n in enumerate(kpt.c_on):
        f_o[o] = np.dot(c_n.conj(), c_n)
        eps_o[o] = np.dot(np.abs(c_n)**2 / f_o[o], kpt.eps_n) # XXX use dscf_hamiltonian_elements for accuracy

        for n, psit_G in enumerate(np.asarray(kpt.psit_nG)):
            wf_oG[o] += c_n[n] / f_o[o]**0.5 * psit_G

        for a, P_oi in P_aoi.items():
            P_oi[o] += np.sum(c_n[:,np.newaxis] / f_o[o]**0.5 * kpt.P_ani[a], axis=0)

    return (f_o, eps_o, wf_oG, P_aoi,)


def dscf_save_band(filename, paw, n):
    """Extract and save all information for band `n` to a tar file."""
    world, bd, gd, kd = paw.wfs.world, paw.wfs.bd, paw.wfs.gd, paw.wfs.kd
    if world.rank == 0:
        # Minimal amount of information needed:
        w = Writer(filename)
        w.dimension('nspins', kd.nspins)
        w.dimension('nibzkpts', kd.nibzkpts)
        w.dimension('nproj', sum([setup.ni for setup in paw.wfs.setups]))
        ng = gd.get_size_of_global_array()
        w.dimension('ngptsx', ng[0])
        w.dimension('ngptsy', ng[1])
        w.dimension('ngptsz', ng[2])

    # Write projections:
    if world.rank == 0:
        w.add('Projection', ('nspins', 'nibzkpts', 'nproj'), dtype=paw.wfs.dtype)
    for s in range(kd.nspins):
        for k in range(kd.nibzkpts):
            all_P_ni = paw.wfs.collect_projections(k, s) # gets all bands
            if world.rank == 0:
                w.fill(all_P_ni[n])

    # Write wave functions:
    if world.rank == 0:
        w.add('PseudoWaveFunction', ('nspins', 'nibzkpts',
                                     'ngptsx', 'ngptsy', 'ngptsz'),
              dtype=paw.wfs.dtype)
    for s in range(kd.nspins):
        for k in range(kd.nibzkpts):
            psit_G = paw.wfs.get_wave_function_array(n, k, s)
            if world.rank == 0:
                w.fill(psit_G)

    if world.rank == 0:
        # Close the file here to ensure that the last wave function is
        # written to disk:
        w.close()

    # We don't want the slaves to start reading before the master has
    # finished writing:
    world.barrier()


def dscf_load_band(filename, paw, molecule=None):
    """Load and distribute all information for a band from a tar file."""
    if not paw.wfs:
        paw.initialize()
    world, bd, gd, kd = paw.wfs.world, paw.wfs.bd, paw.wfs.gd, paw.wfs.kd
    if bd.comm.size != 1:
        raise NotImplementedError('Undefined action for band parallelization.')

    r = Reader(filename)
    assert (r.dimension('nspins') == kd.nspins and \
        r.dimension('nibzkpts') == kd.nibzkpts), 'Incompatible spin/kpoints.'

    # Read wave function for every spin/kpoint owned by this rank
    psit_uG = gd.empty(kd.mynks, kd.dtype)
    for myu, psit_G in enumerate(psit_uG):
        u = kd.global_index(myu)
        s, k = kd.what_is(u)
        if gd.comm.rank == 0:
            big_psit_G = np.array(r.get('PseudoWaveFunction', s, k), paw.wfs.dtype)
        else:
            big_psit_G = None
        gd.distribute(big_psit_G, psit_G)

    # Find domain ranks for each atom
    atoms = paw.get_atoms()
    spos_ac = atoms.get_scaled_positions() % 1.0
    rank_a = gd.get_ranks_from_positions(spos_ac)
    #my_atom_indices = np.argwhere(rank_a == gd.comm.rank).ravel()
    #assert np.all(my_atom_indices == paw.wfs.pt.my_atom_indices)
    assert r.dimension('nproj') == sum([setup.ni for setup in paw.wfs.setups])

    if molecule is None:
        molecule = range(len(atoms))

    # Read projections for every spin/kpoint and atom owned by this rank
    P_uai = [{}] * kd.mynks #[paw.wfs.pt.dict() for myu in range(kd.mynks)]
    for myu, P_ai in enumerate(P_uai):
        u = kd.global_index(myu)
        s, k = kd.what_is(u)
        P_i = r.get('Projection', s, k)
        i1 = 0
        for a in molecule:
            setup = paw.wfs.setups[a]
            i2 = i1 + setup.ni
            if gd.comm.rank == rank_a[a]:
                P_ai[a] = np.array(P_i[i1:i2], paw.wfs.dtype)
            i1 = i2

    return psit_uG, P_uai


def dscf_collapse_orbitals(paw, nbands_max='occupied', f_tol=1e-4,
                           verify_density=True, nt_tol=1e-5, D_tol=1e-3):

    bd, gd, kd = paw.wfs.bd, paw.wfs.gd, paw.wfs.kd
    if bd.comm.size != 1:
        raise NotImplementedError('Undefined action for band parallelization.')

    f_skn = np.empty((kd.nspins, kd.nibzkpts, bd.nbands), dtype=float)
    for s, f_kn in enumerate(f_skn):
        for k, f_n in enumerate(f_kn):
            kpt_rank, myu = kd.get_rank_and_index(s, k)
            if kd.comm.rank == kpt_rank:
                f_n[:] = paw.wfs.kpt_u[myu].f_n
            kd.comm.broadcast(f_n, kpt_rank)

    # Find smallest band index, from which all bands have negligeble occupations
    n0 = np.argmax(f_skn<f_tol, axis=-1).max()
    assert np.all(f_skn[...,n0:]<f_tol) # XXX use f_skn[...,n0:].sum()<f_tol

    # Read the number of Delta-SCF orbitals
    norbitals = paw.occupations.norbitals
    if debug: mpi_debug('n0=%d, norbitals=%d, bd:%d, gd:%d, kd:%d' % (n0,norbitals,bd.comm.size,gd.comm.size,kd.comm.size))

    if nbands_max < 0:
        nbands_max = n0 + norbitals - nbands_max
    elif nbands_max == 'occupied':
        nbands_max = n0 + norbitals

    assert nbands_max >= n0 + norbitals, 'Too few bands to include occupations.'
    ncut = nbands_max-norbitals

    if debug: mpi_debug('nbands_max=%d' % nbands_max) 

    paw.wfs.initialize_wave_functions_from_restart_file() # hurts memmory

    for kpt in paw.wfs.kpt_u:
        mol = kpt.P_ani.keys() # XXX stupid
        (f_o, eps_o, wf_oG, P_aoi,) = dscf_reconstruct_orbitals_k_point(paw, norbitals, mol, kpt)

        assert np.abs(f_o-1).max() < 1e-9, 'Orbitals must be properly normalized.'
        f_o = kpt.ne_o # actual ocupatiion numbers

        # Crop band-data and inject data for Delta-SCF orbitals
        kpt.f_n = np.hstack((kpt.f_n[:n0], f_o, kpt.f_n[n0:ncut]))
        kpt.eps_n = np.hstack((kpt.eps_n[:n0], eps_o, kpt.eps_n[n0:ncut]))
        for a, P_ni in kpt.P_ani.items():
            kpt.P_ani[a] = np.vstack((P_ni[:n0], P_aoi[a], P_ni[n0:ncut]))

        old_psit_nG = kpt.psit_nG
        kpt.psit_nG = gd.empty(nbands_max, dtype=old_psit_nG.dtype)

        if isinstance(old_psit_nG, FileReference):
            assert old_psit_nG.shape[-3:] == wf_oG.shape[-3:], 'Shape mismatch!'

            # Read band-by-band to save memory as full psit_nG may be large
            for n,psit_G in enumerate(kpt.psit_nG):
                if n < n0:
                    full_psit_G = old_psit_nG[n]
                elif n in range(n0,n0+norbitals):
                    full_psit_G = wf_oG[n-n0]
                else:
                    full_psit_G = old_psit_nG[n-norbitals]
                gd.distribute(full_psit_G, psit_G)
        else:
            kpt.psit_nG[:n0] = old_psit_nG[:n0]
            kpt.psit_nG[n0:n0+norbitals] = wf_oG
            kpt.psit_nG[n0+norbitals:] = old_psit_nG[n0:ncut]

        del kpt.ne_o, kpt.c_on, old_psit_nG

    del paw.occupations.norbitals

    # Change various parameters related to new number of bands
    paw.wfs.bd = BandDescriptor(nbands_max, bd.comm, bd.strided)
    if paw.wfs.eigensolver:
        paw.wfs.eigensolver.initialized = False
    del bd

    # Crop convergence criteria nbands_converge to new number of bands
    par = paw.input_parameters
    if 'convergence' in par:
        cc = par['convergence']
        if 'bands' in cc:
            cc['bands'] = min(nbands_max, cc['bands'])

    # Replace occupations class with a fixed variant (gets the magmom right) XXX?!?
    fermilevel, magmom = paw.occupations.fermilevel, paw.occupations.magmom
    paw.occupations = FermiDirac(paw.occupations.width * Hartree, paw.occupations.fixmagmom)
    paw.occupations.set_fermi_level(fermilevel)
    paw.occupations.magmom = magmom
    del fermilevel, magmom

    # For good measure, self-consistency information should be destroyed
    paw.scf.reset()

    if verify_density:
        paw.initialize_positions()

        # Re-calculate pseudo density and watch for changes
        old_nt_sG = paw.density.nt_sG.copy()
        paw.density.calculate_pseudo_density(paw.wfs)
        if debug: mpi_debug('delta-density: %g' % np.abs(old_nt_sG-paw.density.nt_sG).max())
        assert np.abs(paw.density.nt_sG-old_nt_sG).max() < nt_tol, 'Density changed!'

        # Re-calculate atomic density matrices and watch for changes
        old_D_asp = {}
        for a,D_sp in paw.density.D_asp.items():
            old_D_asp[a] = D_sp.copy()
        paw.wfs.calculate_atomic_density_matrices(paw.density.D_asp)
        if debug: mpi_debug('delta-D_asp: %g' % max([0]+[np.abs(D_sp-old_D_asp[a]).max() for a,D_sp in paw.density.D_asp.items()]))
        for a,D_sp in paw.density.D_asp.items():
            assert np.abs(D_sp-old_D_asp[a]).max() < D_tol, 'Atom %d changed!' % a

