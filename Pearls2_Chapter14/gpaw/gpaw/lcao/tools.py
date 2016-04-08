from time import localtime
from numpy import linalg as la
import cPickle as pickle
import numpy as np

from gpaw.utilities import pack, unpack
from gpaw.utilities.tools import tri2full
from gpaw.utilities.blas import rk, gemm
from gpaw.basis_data import Basis
from gpaw.setup import types2atomtypes
from gpaw.coulomb import CoulombNEW as Coulomb
from gpaw.mpi import world, rank, MASTER, serial_comm
from gpaw import GPAW

from ase.units import Hartree
from ase.calculators.singlepoint import SinglePointCalculator


def get_bf_centers(atoms, basis=None):
    calc = atoms.get_calculator()
    if calc is None or isinstance(calc, SinglePointCalculator):
        symbols = atoms.get_chemical_symbols()
        basis_a = types2atomtypes(symbols, basis, 'dzp')
        nao_a = [Basis(symbol, type).nao
                 for symbol, type in zip(symbols, basis_a)]
    else:
        if not calc.initialized:
            calc.initialize(atoms)
        nao_a = [calc.wfs.setups[a].nao for a in range(len(atoms))]
    pos_ic = []
    for pos, nao in zip(atoms.get_positions(), nao_a):
        pos_ic.extend(pos[None].repeat(nao, 0))
    return np.array(pos_ic)


def get_bfi(calc, a_list):
    """basis function indices from a list of atom indices.
       a_list: atom indices
       Use: get_bfi(calc, [0, 4]) gives the functions indices
       corresponding to atom 0 and 4"""
    bfs_list = []
    for a in a_list:
        M = calc.wfs.basis_functions.M_a[a]
        bfs_list += range(M, M + calc.wfs.setups[a].nao)
    return bfs_list


def get_bfi2(symbols, basis, a_list):
    """Same as get_bfi, but does not require an LCAO calc"""
    basis = types2atomtypes(symbols, basis, default='dzp')
    bfs_list = []
    i = 0
    for a, symbol in enumerate(symbols):
        nao = Basis(symbol, basis[a]).nao
        if a in a_list:
            bfs_list += range(i, i + nao)
        i += nao
    return bfs_list
    

def get_mulliken(calc, a_list):
    """Mulliken charges from a list of atom indices (a_list)."""
    Q_a = {}
    for a in a_list:
        Q_a[a] = 0.0
    for kpt in calc.wfs.kpt_u:
        S_MM = calc.wfs.S_qMM[kpt.q]
        nao = S_MM.shape[0]
        rho_MM = np.empty((nao, nao), calc.wfs.dtype)
        calc.wfs.calculate_density_matrix(kpt.f_n, kpt.C_nM, rho_MM)
        Q_M = np.dot(rho_MM, S_MM).diagonal()
        for a in a_list:
            M1 = calc.wfs.basis_functions.M_a[a]
            M2 = M1 + calc.wfs.setups[a].nao
            Q_a[a] += np.sum(Q_M[M1:M2])
    return Q_a        


def get_realspace_hs(h_skmm, s_kmm, bzk_kc, weight_k,
                     R_c=(0, 0, 0), direction='x', 
                     symmetry={'enabled': False}):

    from gpaw.symmetry import Symmetry
    from ase.dft.kpoints import get_monkhorst_pack_size_and_offset, \
        monkhorst_pack
    
    if symmetry['point_group'] is True:
        raise NotImplementedError, 'Point group symmetry not implemented'

    nspins, nk, nbf = h_skmm.shape[:3]
    dir = 'xyz'.index(direction)
    transverse_dirs = np.delete([0, 1, 2], [dir])
    dtype = float
    if len(bzk_kc) > 1 or np.any(bzk_kc[0] != [0, 0, 0]):
        dtype = complex

    kpts_grid = get_monkhorst_pack_size_and_offset(bzk_kc)[0]

    # kpts in the transport direction
    nkpts_p = kpts_grid[dir]
    bzk_p_kc = monkhorst_pack((nkpts_p,1,1))[:, 0]
    weight_p_k = 1. / nkpts_p

   # kpts in the transverse directions
    bzk_t_kc = monkhorst_pack(tuple(kpts_grid[transverse_dirs]) + (1, ))
    if not 'time_reversal' in symmetry:
        symmetry['time_reversal'] = True
    if symmetry['time_reversal'] is True:
        #XXX a somewhat ugly hack:
        # By default GPAW reduces inversion sym in the z direction. The steps 
        # below assure reduction in the transverse dirs.
        # For now this part seems to do the job, but it may be written
        # in a smarter way in the future.
        symmetry = Symmetry([1], np.eye(3))
        symmetry.prune_symmetries_atoms(np.zeros((1, 3)))
        ibzk_kc, ibzweight_k = symmetry.reduce(bzk_kc)[:2]
        ibzk_t_kc, weights_t_k = symmetry.reduce(bzk_t_kc)[:2]
        ibzk_t_kc = ibzk_t_kc[:, :2]
        nkpts_t = len(ibzk_t_kc)
    else:
        ibzk_kc = bzk_kc.copy()
        ibzk_t_kc = bzk_t_kc
        nkpts_t = len(bzk_t_kc)
        weights_t_k = [1. / nkpts_t for k in range(nkpts_t)]

    h_skii = np.zeros((nspins, nkpts_t, nbf, nbf), dtype)
    if s_kmm is not None:
        s_kii = np.zeros((nkpts_t, nbf, nbf), dtype)

    tol = 7
    for j, k_t in enumerate(ibzk_t_kc):
        for k_p in bzk_p_kc:
            k = np.zeros((3,))
            k[dir] = k_p
            k[transverse_dirs] = k_t
            kpoint_list = [list(np.round(k_kc, tol)) for k_kc in ibzk_kc]
            if list(np.round(k, tol)) not in kpoint_list:
                k = -k # inversion
                index = kpoint_list.index(list(np.round(k,tol)))
                h = h_skmm[:, index].conjugate()
                if s_kmm is not None:
                    s = s_kmm[index].conjugate()
                k=-k
            else: # kpoint in the ibz
                index = kpoint_list.index(list(np.round(k, tol)))
                h = h_skmm[:, index]
                if s_kmm is not None:
                    s = s_kmm[index]

            c_k = np.exp(2.j * np.pi * np.dot(k, R_c)) * weight_p_k
            h_skii[:, j] += c_k * h
            if s_kmm is not None:
                s_kii[j] += c_k * s 
    
    if s_kmm is None:
        return ibzk_t_kc, weights_t_k, h_skii
    else:
        return ibzk_t_kc, weights_t_k, h_skii, s_kii


def remove_pbc(atoms, h, s=None, d=0, centers_ic=None, cutoff=None):
    if h.ndim > 2:
        raise KeyError('You have to run remove_pbc for each '
                       'spin/kpoint seperately.')
    L = atoms.cell[d, d]
    nao = len(h)
    dtype = h.dtype
    if centers_ic is None:
        centers_ic = get_bf_centers(atoms) # requires an attached LCAO calc
    ni = len(centers_ic)
    if nao != ni:
        assert nao == 2 * ni
        centers_ic = np.vstack((centers_ic, centers_ic))
        centers_ic[ni:, d] += L
        if cutoff is None:
            cutoff = L - 1e-3
    elif cutoff is None:
        cutoff = 0.5 * L - 1e-3
    pos_i = centers_ic[:, d]
    for i in range(nao):
        dpos_i = abs(pos_i - pos_i[i])
        mask_i = (dpos_i < cutoff).astype(dtype)
        h[i, :] *= mask_i
        h[:, i] *= mask_i
        if s is not None:
            s[i, :] *= mask_i
            s[:, i] *= mask_i


def dump_hamiltonian(filename, atoms, direction=None, Ef=None):
    h_skmm, s_kmm = get_hamiltonian(atoms)
    if direction is not None:
        d = 'xyz'.index(direction)
        for s in range(atoms.calc.nspins):
            for k in range(atoms.calc.nkpts):
                if s==0:
                    remove_pbc(atoms, h_skmm[s, k], s_kmm[k], d)
                else:
                    remove_pbc(atoms, h_skmm[s, k], None, d)
    
    if atoms.calc.master:
        fd = file(filename,'wb')
        pickle.dump((h_skmm, s_kmm), fd, 2)
        atoms_data = {'cell':atoms.cell, 'positions':atoms.positions,
                      'numbers':atoms.numbers, 'pbc':atoms.pbc}
        
        pickle.dump(atoms_data, fd, 2)
        calc_data ={'weight_k':atoms.calc.weight_k, 
                    'ibzk_kc':atoms.calc.ibzk_kc}
        
        pickle.dump(calc_data, fd, 2)
        fd.close()

    world.barrier()


def dump_hamiltonian_parallel(filename, atoms, direction=None, Ef=None):
    """
        Dump the lcao representation of H and S to file(s) beginning
        with filename. If direction is x, y or z, the periodic boundary
        conditions will be removed in the specified direction. 
        If the Fermi temperature is different from zero,  the
        energy zero-point is taken as the Fermi level.

        Note:
        H and S are parallized over spin and k-points and
        is for now dumped into a number of pickle files. This
        may be changed into a dump to a single file in the future.

    """
    if direction is not None:
        d = 'xyz'.index(direction)

    calc = atoms.calc
    wfs = calc.wfs
    nao = wfs.setups.nao
    nq = len(wfs.kpt_u) // wfs.nspins
    H_qMM = np.empty((wfs.nspins, nq, nao, nao), wfs.dtype)
    calc_data = {'k_q':{},
                 'skpt_qc':np.empty((nq, 3)), 
                 'weight_q':np.empty(nq)}

    S_qMM = wfs.S_qMM
   
    for kpt in wfs.kpt_u:
        calc_data['skpt_qc'][kpt.q] = calc.wfs.kd.ibzk_kc[kpt.k]
        calc_data['weight_q'][kpt.q] = calc.wfs.kd.weight_k[kpt.k]
        calc_data['k_q'][kpt.q] = kpt.k
##         print ('Calc. H matrix on proc. %i: '
##                '(rk, rd, q, k) = (%i, %i, %i, %i)') % (
##             wfs.world.rank, wfs.kd.comm.rank,
##             wfs.gd.domain.comm.rank, kpt.q, kpt.k)
        H_MM = wfs.eigensolver.calculate_hamiltonian_matrix(calc.hamiltonian,
                                                            wfs, 
                                                            kpt)

        H_qMM[kpt.s, kpt.q] = H_MM

        tri2full(H_qMM[kpt.s, kpt.q])
        if kpt.s==0:
            tri2full(S_qMM[kpt.q])
            if direction is not None:
                remove_pbc(atoms, H_qMM[kpt.s, kpt.q], S_qMM[kpt.q], d)
        else:
            if direction is not None:
                remove_pbc(atoms, H_qMM[kpt.s, kpt.q], None, d)
        if calc.occupations.width > 0:
            if Ef is None:
                Ef = calc.occupations.get_fermi_level()
            else:
                Ef = Ef / Hartree

            H_qMM[kpt.s, kpt.q] -= S_qMM[kpt.q] * Ef
    
    if wfs.gd.comm.rank == 0:
        fd = file(filename+'%i.pckl' % wfs.kd.comm.rank, 'wb')
        H_qMM *= Hartree
        pickle.dump((H_qMM, S_qMM),fd , 2)
        pickle.dump(calc_data, fd, 2) 
        fd.close()


def get_lcao_hamiltonian(calc):
    """Return H_skMM, S_kMM on master, (None, None) on slaves. H is in eV."""
    if calc.wfs.S_qMM is None:
        calc.wfs.set_positions(calc.get_atoms().get_scaled_positions() % 1)
    dtype = calc.wfs.dtype
    NM = calc.wfs.eigensolver.nao
    Nk = calc.wfs.kd.nibzkpts
    Ns = calc.wfs.nspins
    
    S_kMM = np.zeros((Nk, NM, NM), dtype)
    H_skMM = np.zeros((Ns, Nk, NM, NM), dtype)
    for kpt in calc.wfs.kpt_u:
        H_MM = calc.wfs.eigensolver.calculate_hamiltonian_matrix(
            calc.hamiltonian, calc.wfs, kpt)
        if kpt.s == 0:
            S_kMM[kpt.k] = calc.wfs.S_qMM[kpt.q]
            tri2full(S_kMM[kpt.k])
        H_skMM[kpt.s, kpt.k] = H_MM * Hartree
        tri2full(H_skMM[kpt.s, kpt.k])
    calc.wfs.kd.comm.sum(S_kMM, MASTER)
    calc.wfs.kd.comm.sum(H_skMM, MASTER)
    if rank == MASTER:
        return H_skMM, S_kMM
    else:
        return None, None


def get_lead_lcao_hamiltonian(calc, direction='x'):
    H_skMM, S_kMM = get_lcao_hamiltonian(calc)
    if rank == MASTER:
        return lead_kspace2realspace(H_skMM, S_kMM,
                                     bzk_kc=calc.wfs.kd.bzk_kc,
                                     weight_k=calc.wfs.kd.weight_k,
                                     direction=direction,
                                     symmetry=calc.input_parameters['symmetry'])
    else:
        return None, None, None, None


def lead_kspace2realspace(h_skmm, s_kmm, bzk_kc, weight_k, direction='x', 
                          symmetry={'point_group': False}):
    """Convert a k-dependent Hamiltonian to tight-binding onsite and coupling.

    For each transverse k-point:
    Convert k-dependent (in transport direction) Hamiltonian
    representing a lead to a real space tight-binding Hamiltonian
    of double size representing two principal layers and the coupling between.
    """

    dir = 'xyz'.index(direction)
    if symmetry['point_group'] is True:
        raise NotImplementedError

    R_c = [0, 0, 0]
    ibz_t_kc, weight_t_k, h_skii, s_kii = get_realspace_hs(
        h_skmm, s_kmm, bzk_kc, weight_k, R_c, direction, symmetry)

    R_c[dir] = 1
    h_skij, s_kij = get_realspace_hs(
        h_skmm, s_kmm, bzk_kc, weight_k, R_c, direction, symmetry)[-2:]

    nspins, nk, nbf = h_skii.shape[:-1]

    h_skmm = np.zeros((nspins, nk, 2 * nbf, 2 * nbf), h_skii.dtype)
    s_kmm = np.zeros((nk, 2 * nbf, 2 * nbf), h_skii.dtype)
    h_skmm[:, :, :nbf, :nbf] = h_skmm[:, :, nbf:, nbf:] = h_skii
    h_skmm[:, :, :nbf, nbf:] = h_skij
    h_skmm[:, :, nbf:, :nbf] = h_skij.swapaxes(2, 3).conj()

    s_kmm[:, :nbf, :nbf] = s_kmm[:, nbf:, nbf:] = s_kii
    s_kmm[:, :nbf, nbf:] = s_kij
    s_kmm[:, nbf:, :nbf] = s_kij.swapaxes(1,2).conj()

    return ibz_t_kc, weight_t_k, h_skmm, s_kmm


def zeta_pol(basis):
    """Get number of zeta func. and polarization func. indices in Basis."""
    zeta = []
    pol = []
    for bf in basis.bf_j:
        if 'polarization' in bf.type:
            pol.append(2 * bf.l + 1)
        else:
            zeta.append(2 * bf.l + 1)
    zeta = sum(zeta)
    pol = sum(pol)
    assert zeta + pol == basis.nao
    return zeta, pol


def basis_subset(symbol, largebasis, smallbasis):
    """Title.

    Determine which basis function indices from ``largebasis`` are also
    present in smallbasis.
    """
    blarge = Basis(symbol, largebasis)
    zeta_large, pol_large = zeta_pol(blarge)
    
    bsmall = Basis(symbol, smallbasis)
    zeta_small, pol_small = zeta_pol(bsmall)

    assert zeta_small <= zeta_large
    assert pol_small <= pol_large

    insmall = np.zeros(blarge.nao, bool)
    insmall[:zeta_small] = True
    insmall[zeta_large:zeta_large + pol_small] = True
    return insmall


def basis_subset2(symbols, largebasis='dzp', smallbasis='sz'):
    """Same as basis_subset, but for an entire list of atoms."""
    largebasis = types2atomtypes(symbols, largebasis, default='dzp')
    smallbasis = types2atomtypes(symbols, smallbasis, default='sz')
    mask = []
    for symbol, large, small in zip(symbols, largebasis, smallbasis):
        mask.extend(basis_subset(symbol, large, small))
    return np.asarray(mask, bool)


def collect_orbitals(a_xo, coords, root=0):
    """Collect array distributed over orbitals to root-CPU.

    Input matrix has last axis distributed amongst CPUs,
    return is None on slaves, and the collected array on root.

    The distribution can be uneven amongst CPUs. The list coords gives the
    number of values for each CPU.
    """
    a_xo = np.ascontiguousarray(a_xo)
    if world.size == 1:
        return a_xo

    # All slaves send their piece to ``root``:
    # There can be several sends before the corresponding receives
    # are posted, so use syncronous send here
    if world.rank != root:
        world.ssend(a_xo, root, 112)
        return None

    # On root, put the subdomains from the slaves into the big array
    # for the whole domain on root:
    xshape = a_xo.shape[:-1]
    Norb2 = sum(coords) # total number of orbital indices
    a_xO = np.empty(xshape + (Norb2,), a_xo.dtype)
    o = 0
    for rank, norb in enumerate(coords):
        if rank != root:
            tmp_xo = np.empty(xshape + (norb,), a_xo.dtype)
            world.receive(tmp_xo, rank, 112)
            a_xO[..., o:o + norb] = tmp_xo
        else:
            a_xO[..., o:o + norb] = a_xo
        o += norb
    return a_xO



def makeU(gpwfile='grid.gpw', orbitalfile='w_wG__P_awi.pckl',
          rotationfile='eps_q__U_pq.pckl', tolerance=1e-5,
          writeoptimizedpairs=False, dppname='D_pp.pckl', S_w=None):

    # S_w: None or diagonal of overlap matrix. In the latter case
    # the optimized and truncated pair orbitals are obtained from
    # normalized (to 1) orbitals.
    #     
    # Tolerance is used for truncation of optimized pairorbitals
    #calc = GPAW(gpwfile, txt=None)
    from gpaw import GPAW
    from gpaw.utilities import pack, unpack
    from gpaw.utilities.blas import rk, gemm
    from gpaw.mpi import world, MASTER
    calc = GPAW(gpwfile, txt='pairorb.txt') # XXX
    gd = calc.wfs.gd
    setups = calc.wfs.setups
    myatoms = calc.density.D_asp.keys()
    del calc

    # Load orbitals on master and distribute to slaves
    if world.rank == MASTER:
        wglobal_wG, P_awi = pickle.load(open(orbitalfile))
        Nw = len(wglobal_wG)
        print('Estimated total (serial) mem usage: %0.3f GB' % (
            np.prod(gd.N_c) * Nw**2 * 8 / 1024.**3))
    else:
        wglobal_wG = None
        Nw = 0
    Nw = gd.comm.sum(Nw) #distribute Nw to all nodes
    w_wG = gd.empty(n=Nw)
    gd.distribute(wglobal_wG, w_wG)
    del wglobal_wG
    
    # Make pairorbitals
    f_pG = gd.zeros(n=Nw**2)
    Np = len(f_pG)
    for p, (w1, w2) in enumerate(np.ndindex(Nw, Nw)):
        np.multiply(w_wG[w1], w_wG[w2], f_pG[p])
    del w_wG
    assert f_pG.flags.contiguous 
    # Make pairorbital overlap (lower triangle only)
    D_pp = np.zeros((Nw**2, Nw**2))
    rk(gd.dv, f_pG, 0., D_pp)

    # Add atomic corrections to pairorbital overlap
    for a in myatoms:
        if setups[a].type != 'ghost':
            P_pp = np.array([pack(np.outer(P_awi[a][w1], P_awi[a][w2]))
                             for w1, w2 in np.ndindex(Nw, Nw)])
            I4_pp = setups[a].four_phi_integrals()
            A = np.zeros((len(I4_pp), len(P_pp)))
            gemm(1.0, P_pp, I4_pp, 0.0, A, 't')
            gemm(1.0, A, P_pp, 1.0, D_pp)
            #D_pp += np.dot(P_pp, np.dot(I4_pp, P_pp.T))
   
    # Summ all contributions to master
    gd.comm.sum(D_pp, MASTER)

    if world.rank == MASTER:
        if S_w is not None:
            print('renormalizing pairorb overlap matrix (D_pp)')
            S2 = np.sqrt(S_w)
            for pa, (wa1, wa2) in enumerate(np.ndindex(Nw, Nw)):
                for pb, (wb1, wb2) in enumerate(np.ndindex(Nw, Nw)):
                    D_pp[pa, pb] /= S2[wa1] * S2[wa2] * S2[wb1] * S2[wb2]

        D_pp.dump(dppname) # XXX if the diagonalization below (on MASTER only)
                           # fails, then one can always restart the stuff
                           # below using only the stored D_pp matrix

        # Determine eigenvalues and vectors on master only
        eps_q, U_pq = np.linalg.eigh(D_pp, UPLO='L')
        del D_pp
        indices = np.argsort(-eps_q.real)
        eps_q = np.ascontiguousarray(eps_q.real[indices])
        U_pq = np.ascontiguousarray(U_pq[:, indices])

        # Truncate
        indices = eps_q > tolerance
        U_pq = np.ascontiguousarray(U_pq[:, indices])
        eps_q = np.ascontiguousarray(eps_q[indices])

        # Dump to file
        pickle.dump((eps_q, U_pq), open(rotationfile, 'wb'), 2)

    if writeoptimizedpairs is not False:
        assert world.size == 1 # works in parallel if U and eps are broadcast
        Uisq_qp = (U_pq / np.sqrt(eps_q)).T.copy()
        g_qG = gd.zeros(n=len(eps_q))
        gemm(1.0, f_pG, Uisq_qp, 0.0, g_qG)
        g_qG = gd.collect(g_qG)
        if world.rank == MASTER:
            P_app = dict([(a, np.array([pack(np.outer(P_wi[w1], P_wi[w2]),
                                             tolerance=1e3)
                                        for w1, w2 in np.ndindex(Nw, Nw)]))
                          for a, P_wi in P_awi.items()])
            P_aqp = dict([(a, np.dot(Uisq_qp, P_pp))
                          for a, P_pp in P_app.items()])
            pickle.dump((g_qG, P_aqp), open(writeoptimizedpairs, 'wb'), 2)


def makeV(gpwfile='grid.gpw', orbitalfile='w_wG__P_awi.pckl',
          rotationfile='eps_q__U_pq.pckl', coulombfile='V_qq.pckl',
          log='V_qq.log', fft=False):
    
    if isinstance(log, str) and world.rank == MASTER:
        log = open(log, 'w')

    # Extract data from files
    calc = GPAW(gpwfile, txt=None, communicator=serial_comm)
    spos_ac = calc.get_atoms().get_scaled_positions() % 1.0
    coulomb = Coulomb(calc.wfs.gd, calc.wfs.setups, spos_ac, fft)
    w_wG, P_awi = pickle.load(open(orbitalfile))
    eps_q, U_pq = pickle.load(open(rotationfile))
    del calc

    # Make rotation matrix divided by sqrt of norm
    Nq = len(eps_q)
    Np = len(U_pq)
    Ni = len(w_wG)
    Uisq_iqj = (U_pq/np.sqrt(eps_q)).reshape(Ni, Ni, Nq).swapaxes(1, 2).copy()
    del eps_q, U_pq

    # Determine number of opt. pairorb on each cpu
    Ncpu = world.size
    nq, R = divmod(Nq, Ncpu)
    nq_r = nq * np.ones(Ncpu, int)
    if R > 0:
        nq_r[-R:] += 1

    # Determine number of opt. pairorb on this cpu
    nq1 = nq_r[world.rank]
    q1end = nq_r[:world.rank + 1].sum()
    q1start = q1end - nq1
    V_qq = np.zeros((Nq, nq1), float)

    def make_optimized(qstart, qend):
        g_qG = np.zeros((qend - qstart,) + w_wG.shape[1:], float)
        P_aqp = {}
        for a, P_wi in P_awi.items():
            ni = P_wi.shape[1]
            nii = ni * (ni + 1) // 2
            P_aqp[a] = np.zeros((qend - qstart, nii), float)
        for w1, w1_G in enumerate(w_wG):
            U = Uisq_iqj[w1, qstart: qend].copy()
            gemm(1., w1_G * w_wG, U, 1.0, g_qG)
            for a, P_wi in P_awi.items():
                P_wp = np.array([pack(np.outer(P_wi[w1], P_wi[w2])) 
                                for w2 in range(Ni)])
                gemm(1., P_wp, U, 1.0, P_aqp[a])
        return g_qG, P_aqp

    g1_qG, P1_aqp = make_optimized(q1start, q1end)
    for block, nq2 in enumerate(nq_r):
        if block == world.rank:
            g2_qG, P2_aqp = g1_qG, P1_aqp
            q2start, q2end = q1start, q1end
        else:
            q2end = nq_r[:block + 1].sum()
            q2start = q2end - nq2
            g2_qG, P2_aqp = make_optimized(q2start, q2end)

        for q1, q2 in np.ndindex(nq1, nq2):
            P1_ap = dict([(a, P_qp[q1]) for a, P_qp in P1_aqp.items()])
            P2_ap = dict([(a, P_qp[q2]) for a, P_qp in P2_aqp.items()])
            V_qq[q2 + q2start, q1] = coulomb.calculate(g1_qG[q1], g2_qG[q2],
                                                       P1_ap, P2_ap)
            if q2 == 0 and world.rank == MASTER:
                T = localtime()
                log.write(
                    'Block %i/%i is %4.1f percent done at %02i:%02i:%02i\n' % (
                    block + 1, world.size, 100.0 * q1 / nq1, T[3], T[4], T[5]))
                log.flush()

    # Collect V_qq array on master node
    if world.rank == MASTER:
        T = localtime()
        log.write('Starting collect at %02i:%02i:%02i\n' % (
            T[3], T[4], T[5]))
        log.flush()

    V_qq = collect_orbitals(V_qq, coords=nq_r, root=MASTER)
    if world.rank == MASTER:
        # V can be slightly asymmetric due to numerics
        V_qq = 0.5 * (V_qq + V_qq.T)
        V_qq.dump(coulombfile)

        T = localtime()
        log.write('Finished at %02i:%02i:%02i\n' % (
            T[3], T[4], T[5]))
        log.flush()
