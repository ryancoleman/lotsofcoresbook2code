import copy
import cPickle
import numpy as np

from ase.utils.timing import Timer
from ase.units import Hartree

from gpaw.utilities import unpack
from gpaw.mpi import world, rank, send_string, receive_string, broadcast_string
from gpaw.utilities.blas import gemm
from gpaw.utilities.lapack import inverse_general
import _gpaw

def tw(mat, filename):
    fd = file(filename, 'wb')
    cPickle.dump(mat, fd, 2)
    fd.close()

def tr(filename):
    fd = file(filename, 'r')
    mat = cPickle.load(fd)
    fd.close()
    return mat

def write(filename, name, data, dimension, dtype=float):
    import gpaw.io.tar as io
    if world.rank == 0:
        w = io.Writer(filename)
        dim = ()
        for i in range(len(dimension)):
            w.dimension(str(i), dimension[i])
            dim += (str(i),)
        w.add(name, dim, dtype=dtype)
        w.fill(data)
        w.close()

def fermidistribution(energy, kt):
    #fermi level is fixed to zero
    return 1.0 / (1.0 + np.exp(energy / kt))

def zeroTFermi(energy):
    if np.real(energy) > 0:
        return 0.
    else:
        return 1.

def get_tri_type(mat):
    #mat is lower triangular or upper triangular matrix
    tol = 1e-10
    mat = abs(mat)
    dim = mat.shape[-1]
    sum = [0, 0]
    for i in range(dim):
        sum[0] += np.trace(mat, -j)
        sum[1] += np.trace(mat, j)
    diff = sum[0] - sum[1]
    if diff >= 0:
        ans = 'L'
    elif diff < 0:
        ans = 'U'
    if abs(diff) < tol:
        print('Warning: can not define the triangular matrix')
    return ans
    
def tri2full(M,UL='L'):
    """UP='L' => fill upper triangle from lower triangle
       such that M=M^d"""
    nbf = len(M)
    if UL=='L':
        for i in range(nbf-1):
            M[i,i:] = M[i:,i].conjugate()
    elif UL=='U':
        for i in range(nbf-1):
            M[i:,i] = M[i,i:].conjugate()

def dagger(matrix):
    return np.conj(matrix.T)

def get_matrix_index(ind1, ind2=None):
    if ind2 is None:
        dim1 = len(ind1)
        return np.resize(ind1, (dim1, dim1))
    else:
        dim1 = len(ind1)
        dim2 = len(ind2)
    return np.resize(ind1, (dim2, dim1)).T, np.resize(ind2, (dim1, dim2))
    
def aa1d(a, d=2):
    # array average in one dimension
    dim = a.shape
    b = [np.sum(np.take(a, [i], axis=d)) for i in range(dim[d])]
    b = np.array(b)
    b = (b * dim[d]) / np.product(dim)
    return b
    
def aa2d(a, d=0):
    # array average in two dimensions
    b = np.sum(a, axis=d) / a.shape[d]
    return b

def k2r_hs(h_skmm, s_kmm, ibzk_kc, weight_k, R_c=(0,0,0), magnet=None):
    phase_k = np.dot(2 * np.pi * ibzk_kc, R_c)
    c_k = np.exp(1.0j * phase_k) * weight_k
    c_k.shape = (len(ibzk_kc),1,1)

    if h_skmm is not None:
        nbf = h_skmm.shape[-1]
        nspins = len(h_skmm)
        h_smm = np.empty((nspins, nbf, nbf),complex)
        for s in range(nspins):
            h_smm[s] = np.sum((h_skmm[s] * c_k), axis=0)
    if s_kmm is not None:
        nbf = s_kmm.shape[-1]
        s_mm = np.empty((nbf, nbf),complex)
        s_mm[:] = np.sum((s_kmm * c_k), axis=0)     
    #if magnet is not None:
    #    MM = magnet.trans_matrix(diag=True)
    #    assert np.sum(R_c) == R_c[2]
    #    MM = MM ** R_c[2]
    #    if h_skmm is not None:
    #        for s in range(nspins):
    #            h_smm[s] *= MM
    #    if s_kmm is not None:  
    #        s_mm *= MM    
    if h_skmm is not None and s_kmm is not None:
        return h_smm, s_mm
    elif h_skmm is None:
        return s_mm
    elif s_kmm is None:
        return h_smm

def r2k_hs(h_srmm, s_rmm, R_vector, kvector=(0,0,0), magnet=None):
    phase_k = np.dot(2 * np.pi * R_vector, kvector)
    c_k = np.exp(-1.0j * phase_k)
    c_k.shape = (len(R_vector), 1, 1)
   
    if h_srmm is not None:
        nbf = h_srmm.shape[-1]
        nspins = len(h_srmm)
        h_smm = np.empty((nspins, nbf, nbf), complex)
        for s in range(nspins):
            h_smm[s] = np.sum((h_srmm[s] * c_k), axis=0)
    if s_rmm is not None:
        nbf = s_rmm.shape[-1]
        s_mm = np.empty((nbf, nbf), complex)
        s_mm[:] = np.sum((s_rmm * c_k), axis=0)
#    if magnet is not None:
#        MM = magnet.trans_matrix(diag=True)
#       assert np.sum(R_vector) == R_vector[2]
#        MM = MM ** (-R_vector[2])
#       if h_srmm is not None:
#            for s in range(nspins):
#               h_smm[s] *= MM
#       if s_mm is not None:    
#           s_mm *= MM    
    if h_srmm is not None and s_rmm is not None:   
        return h_smm, s_mm
    elif h_srmm is None:
        return s_mm
    elif s_rmm is None:
        return h_smm

def r2k2(s_rmm, R_vector, kvector=(0,0,0), symmetrize=True):
    phase_k = np.dot(2 * np.pi * R_vector, kvector)
    c_k = np.exp(-1.0j * phase_k)
    c_k.shape = (len(R_vector), 1, 1)
    nbf = s_rmm.shape[-1]
    s_mm = np.zeros((nbf, nbf), complex)
    for i in range(len(R_vector)):
        tmp = s_rmm[i,:,:]*c_k[i]
        if symmetrize and R_vector[i,:].any():
            s_mm = s_mm + tmp + dagger(tmp)
        else:
            s_mm = s_mm + tmp
    return s_mm

def collect_lead_mat(lead_hsd, lead_couple_hsd, s, pk, flag='S'):
    diag_h = []
    upc_h = []
    dwnc_h = []
    for i, hsd, c_hsd in zip(range(len(lead_hsd)), lead_hsd, lead_couple_hsd):
        if flag == 'S':
            band_mat, cp_mat = hsd.S[pk], c_hsd.S[pk]
        elif flag == 'H':
            band_mat, cp_mat = hsd.H[s][pk], c_hsd.H[s][pk]
        else:
            band_mat, cp_mat = hsd.D[s][pk], c_hsd.D[s][pk]
        diag_h.append(copy.deepcopy(band_mat))
        upc_h.append(cp_mat.recover('c'))
        dwnc_h.append(cp_mat.recover('n'))
    return diag_h, upc_h, dwnc_h        
        
def get_hs(atoms):
    """Calculate the Hamiltonian and overlap matrix."""
    calc = atoms.calc
    wfs = calc.wfs
    wfs.gd.comm.broadcast(wfs.S_qMM, 0)    
    Ef = calc.get_fermi_level()
    eigensolver = wfs.eigensolver
    ham = calc.hamiltonian
    S_qMM = wfs.S_qMM.copy()
    for S_MM in S_qMM:
        tri2full(S_MM)
    H_sqMM = np.empty((wfs.nspins,) + S_qMM.shape, wfs.dtype)
    for kpt in wfs.kpt_u:
        H_MM = eigensolver.calculate_hamiltonian_matrix(ham, wfs, kpt)
        tri2full(H_MM)
        H_MM *= Hartree
        #H_MM -= Ef * S_qMM[kpt.q]
        H_sqMM[kpt.s, kpt.q] = H_MM
    wfs.gd.comm.broadcast(H_sqMM, 0)        
    return H_sqMM, S_qMM

def substract_pk(d, npk, ntk, kpts, k_mm, hors='s', position=[0, 0, 0], magnet=None):
    weight = np.array([1.0 / ntk] * ntk )
    if hors not in 'hs':
        raise KeyError('hors should be h or s!')
    if hors == 'h':
        dim = k_mm.shape[:]
        dim = (dim[0],) + (dim[1] // ntk,) + dim[2:]
        pk_mm = np.empty(dim, k_mm.dtype)
        dim = (dim[0],) + (ntk,) + dim[2:]
        tk_mm = np.empty(dim, k_mm.dtype)
    elif hors == 's':
        dim = k_mm.shape[:]
        dim = (dim[0] // ntk,) + dim[1:]
        pk_mm = np.empty(dim, k_mm.dtype)
        dim = (ntk,) + dim[1:]
        tk_mm = np.empty(dim, k_mm.dtype)

    tkpts = pick_out_tkpts(d, npk, ntk, kpts)
    for i in range(npk):
        n = i * ntk
        for j in range(ntk):
            if hors == 'h':
                tk_mm[:, j] = np.copy(k_mm[:, n + j])
            elif hors == 's':
                tk_mm[j] = np.copy(k_mm[n + j])
        if hors == 'h':
            pk_mm[:, i] = k2r_hs(tk_mm, None, tkpts, weight, position, magnet)
        elif hors == 's':
            pk_mm[i] = k2r_hs(None, tk_mm, tkpts, weight, position, magnet)
    return pk_mm   

def pick_out_tkpts(d, npk, ntk, kpts):
    tkpts = np.zeros([ntk, 3])
    for i in range(ntk):
        tkpts[i, d] = kpts[i, d]
    return tkpts

def count_tkpts_num(d, kpts):
    tol = 1e-6
    tkpts = [kpts[0]]
    for kpt in kpts:
        flag = False
        for tkpt in tkpts:
            if abs(kpt[d] - tkpt[d]) < tol:
                flag = True
        if not flag:
            tkpts.append(kpt)
    return len(tkpts)
    
def dot(a, b, transa='n'):
    assert len(a.shape) == 2 and a.shape[1] == b.shape[0]
    dtype = complex
    if a.dtype == complex and b.dtype == complex:
        c = a
        d = b
    elif a.dtype == float and b.dtype == complex:
        c = np.array(a, complex)
        d = b
    elif a.dtype == complex and b.dtype == float:
        d = np.array(b, complex)
        c = a
    else:
        dtype = float
        c = a
        d = b
    e = np.zeros([c.shape[0], d.shape[1]], dtype)
    gemm(1.0, np.ascontiguousarray(d), np.ascontiguousarray(c), 0.0, e, transa)
    return e

def gcd(m,n):
    while n:
        m,n=n,m%n
    return m

def plot_diag(mtx, ind=1):
    import pylab
    dim = mtx.shape
    if len(dim) != 2:
        print('Warning! check the dimenstion of the matrix')
    if dim[0] != dim[1]:
        print('Warinng! check if the matrix is square')
    diag_element = np.diag(mtx)
    y_data = pick(diag_element, ind)
    x_data = range(len(y_data))
    pylab.plot(x_data, y_data,'b-o')
    pylab.show()

def get_atom_indices(subatoms, setups):
    basis_list = [setup.nao for setup in setups]
    index = []
    for j, lj  in zip(subatoms, range(len(subatoms))):
        begin = np.sum(np.array(basis_list[:j], int))
        for n in range(basis_list[j]):
            index.append(begin + n) 
    return np.array(index, int)    

def mp_distribution(e, kt, n=1):
    x = e / kt
    re = 0.5 * error_function(x)
    for i in range(n):
        re += coff_function(i + 1) * hermite_poly(2 * i + 1, x) * np.exp(-x**2) 
    return re        

def coff_function(n):
    return (-1)**n / (np.product(np.arange(1, n + 1)) * 4.** n * np.sqrt(np.pi))
    
def hermite_poly(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2 * x
    else:
        return 2 * x * hermite_poly(n - 1, x) \
                                      - 2 * (n - 1) * hermite_poly(n - 2 , x)

def error_function(x):
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * np.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
        t*(.09678418+t*(-.18628806+t*(.27886807+
        t*(-1.13520398+t*(1.48851587+t*(-.82215223+
        t*.17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r

def sum_by_unit(x, unit):
    dim = x.shape[0]
    dim1 = int(np.ceil(dim / unit))
    y = np.empty([dim1], dtype=x.dtype)
    for i in range(dim1 - 1):
        y[i] = np.sum(x[i * unit: (i + 1) * unit]) / unit
    y[0] = y[1]
    y[-1] = y[-2]
    return y

def diag_cell(cell):
    if len(cell.shape) == 2:
        cell = np.diag(cell)
    return cell
    
def get_pk_hsd(d, ntk, kpts, hl_skmm, sl_kmm, dl_skmm, txt=None,
                                                  dtype=complex, direction=0, magnet=None):
    npk = len(kpts) // ntk
    position = [0, 0, 0]
    hl_spkmm = substract_pk(d, npk, ntk, kpts, hl_skmm, hors='h', magnet=magnet)
    dl_spkmm = substract_pk(d, npk, ntk, kpts, dl_skmm, hors='h', magnet=magnet)
    sl_pkmm = substract_pk(d, npk, ntk, kpts, sl_kmm, hors='s', magnet=magnet)
    
    if direction==0:
        position[d] = 1
    else:
        position[d] = -1
    
    hl_spkcmm = substract_pk(d, npk, ntk, kpts, hl_skmm, 'h', position, magnet=magnet)
    dl_spkcmm = substract_pk(d, npk, ntk, kpts, dl_skmm, 'h', position, magnet=magnet)
    sl_pkcmm = substract_pk(d, npk, ntk, kpts, sl_kmm, 's', position, magnet=magnet)
    
    tol = 1e-6
    position[d] = 2.0
    s_test = substract_pk(d, npk, ntk, kpts, sl_kmm, 's', position, magnet=magnet)
    
    matmax = np.max(abs(s_test))
    if matmax > tol:
        if txt is not None:
            txt('Warning*: the principle layer should be larger, \
                                                      matmax=%f' % matmax)
        else:
            print('Warning*: the principle layer should be larger, \
                                                      matmax=%f' % matmax)
    if dtype == float:
        hl_spkmm = np.real(hl_spkmm).copy()
        sl_pkmm = np.real(sl_pkmm).copy()
        dl_spkmm = np.real(dl_spkmm).copy()
        hl_spkcmm = np.real(hl_spkcmm).copy()
        sl_pkcmm = np.real(sl_pkcmm).copy()
        dl_spkcmm = np.real(dl_spkcmm).copy()
    return hl_spkmm, sl_pkmm, dl_spkmm * ntk, hl_spkcmm, \
                                                    sl_pkcmm, dl_spkcmm * ntk
   
def get_lcao_density_matrix(calc):
    wfs = calc.wfs
    kpts = wfs.kd.ibzk_qc
    nq = len(kpts)
    my_ns = len(wfs.kpt_u) // nq
    nao = wfs.setups.nao
    
    # calculate_density_matrix involves gemm and doesn't work well with empty()
    d_skmm = np.zeros([my_ns, nq, nao, nao], wfs.dtype)
    for kpt in wfs.kpt_u:
        if my_ns == 1:
            wfs.calculate_density_matrix(kpt.f_n, kpt.C_nM, d_skmm[0, kpt.q])
        else:
            wfs.calculate_density_matrix(kpt.f_n, kpt.C_nM, d_skmm[kpt.s, kpt.q])            
    return d_skmm

def collect_atomic_matrices(asp, setups, ns, comm, rank_a):
    all_asp = []
    for a, setup in enumerate(setups):
        sp = asp.get(a)
        if sp is None:
            ni = setup.ni
            sp = np.empty((ns, ni * (ni + 1) // 2))
        if comm.size > 1:
            comm.broadcast(sp, rank_a[a])
        all_asp.append(sp)      
    return all_asp

def distribute_atomic_matrices(all_asp, asp, setups):
    for a in range(len(setups)):
        if asp.get(a) is not None:
            asp[a] = all_asp[a]    

def collect_and_distribute_atomic_matrices(D_ap, setups, setups0, rank_a, comm, keys):
    gD_ap = []
    D_ap0 = [None] * len(keys)
    for a, setup in enumerate(setups):
        if a not in keys:
            ni = setup.ni
            sp = np.empty((ni * (ni + 1) // 2))
        else:
            sp = D_ap[keys.index(a)]
        if comm.size > 1:
            comm.broadcast(sp, rank_a[a])
        gD_ap.append(sp)
    for a in range(len(setups0)):
        if a in keys:
            D_ap0[keys.index(a)] = gD_ap[a]
    return D_ap0            

def generate_selfenergy_database(atoms, ntk, filename, direction=0, kt=0.1,
                                 bias=[-3,3], depth=3, comm=None):
    from gpaw.transport.sparse_matrix import Banded_Sparse_HSD, CP_Sparse_HSD, Se_Sparse_Matrix
    from gpaw.transport.selfenergy import LeadSelfEnergy
    from gpaw.transport.contour import Contour
    hl_skmm, sl_kmm = get_hs(atoms)
    dl_skmm = get_lcao_density_matrix(atoms.calc)
    fermi = atoms.calc.get_fermi_level()
    wfs = atoms.calc.wfs
    hl_spkmm, sl_pkmm, dl_spkmm,  \
    hl_spkcmm, sl_pkcmm, dl_spkcmm = get_pk_hsd(2, ntk,
                                                wfs.kd.ibzk_qc,
                                                hl_skmm, sl_kmm, dl_skmm,
                                                None, wfs.dtype,
                                                direction=direction)    
    my_npk = len(wfs.kd.ibzk_qc) / ntk
    my_nspins = len(wfs.kpt_u) / ( my_npk * ntk)
    
    lead_hsd = Banded_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    lead_couple_hsd = CP_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    for pk in range(my_npk):
        lead_hsd.reset(0, pk, sl_pkmm[pk], 'S', init=True)
        lead_couple_hsd.reset(0, pk, sl_pkcmm[pk], 'S', init=True)
        for s in range(my_nspins):
            lead_hsd.reset(s, pk, hl_spkmm[s, pk], 'H', init=True)     
            lead_hsd.reset(s, pk, dl_spkmm[s, pk], 'D', init=True)
            lead_couple_hsd.reset(s, pk, hl_spkcmm[s, pk], 'H', init=True)     
            lead_couple_hsd.reset(s, pk, dl_spkcmm[s, pk], 'D', init=True)          
    lead_se = LeadSelfEnergy(lead_hsd, lead_couple_hsd)
    contour = Contour(kt, [fermi] * 2, bias, depth, comm=comm)    
    path = contour.get_plot_path(ex=True)
    for nid, energy in zip(path.my_nids, path.my_energies):
        for kpt in wfs.kpt_u:
                if kpt.q % ntk == 0:
                    flag = str(kpt.k // ntk) + '_' + str(nid)
                    lead_se.s = kpt.s
                    lead_se.pk = kpt.q // ntk
                    data = lead_se(energy)
                    fd = file(flag, 'w')    
                    cPickle.dump(data, fd, 2)
                    fd.close()    

def test_selfenergy_interpolation(atoms, ntk, filename, begin, end, base, scale, direction=0):
    from gpaw.transport.sparse_matrix import Banded_Sparse_HSD, CP_Sparse_HSD, Se_Sparse_Matrix
    from gpaw.transport.selfenergy import LeadSelfEnergy
    from gpaw.transport.contour import Contour
    hl_skmm, sl_kmm = get_hs(atoms)
    dl_skmm = get_lcao_density_matrix(atoms.calc)
    fermi = atoms.calc.get_fermi_level()
    wfs = atoms.calc.wfs
    hl_spkmm, sl_pkmm, dl_spkmm,  \
    hl_spkcmm, sl_pkcmm, dl_spkcmm = get_pk_hsd(2, ntk,
                                                wfs.kd.ibzk_qc,
                                                hl_skmm, sl_kmm, dl_skmm,
                                                None, wfs.dtype,
                                                direction=direction)    
    my_npk = len(wfs.kd.ibzk_qc) / ntk
    my_nspins = len(wfs.kpt_u) / ( my_npk * ntk)
    
    lead_hsd = Banded_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    lead_couple_hsd = CP_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    for pk in range(my_npk):
        lead_hsd.reset(0, pk, sl_pkmm[pk], 'S', init=True)
        lead_couple_hsd.reset(0, pk, sl_pkcmm[pk], 'S', init=True)
        for s in range(my_nspins):
            lead_hsd.reset(s, pk, hl_spkmm[s, pk], 'H', init=True)     
            lead_hsd.reset(s, pk, dl_spkmm[s, pk], 'D', init=True)
            lead_couple_hsd.reset(s, pk, hl_spkcmm[s, pk], 'H', init=True)     
            lead_couple_hsd.reset(s, pk, dl_spkcmm[s, pk], 'D', init=True)          
    lead_se = LeadSelfEnergy(lead_hsd, lead_couple_hsd)
    begin += fermi
    end += fermi
    
    ee = np.linspace(begin, end, base)
    cmp_ee = np.linspace(begin, end, base * scale)
  
    se = []
    cmp_se = []
    from scipy import interpolate

    for e in ee:
        se.append(lead_se(e).recover())
    se = np.array(se)
    ne, ny, nz= se.shape
    nie = len(cmp_ee)
    data = np.zeros([nie, ny, nz], se.dtype)
    for yy in range(ny):
        for zz in range(nz):
            ydata = se[:, yy, zz]
            f = interpolate.interp1d(ee, ydata)
            data[:, yy, zz] = f(cmp_ee)
    inter_se_linear = data
    
    for e in cmp_ee:
        cmp_se.append(lead_se(e).recover())
    
    fd = file(filename, 'w')
    cPickle.dump((cmp_se, inter_se_linear, ee, cmp_ee), fd, 2)
    fd.close()
    
    for i,e in enumerate(cmp_ee):
        print(e, np.max(abs(cmp_se[i] - inter_se_linear[i])), 'linear', np.max(abs(cmp_se[i])))


def path_selfenergy(atoms, ntk, filename, begin, end, num= 257, direction=0):
    from gpaw.transport.sparse_matrix import Banded_Sparse_HSD, CP_Sparse_HSD, Se_Sparse_Matrix
    from gpaw.transport.selfenergy import LeadSelfEnergy
    from gpaw.transport.contour import Contour
    hl_skmm, sl_kmm = get_hs(atoms)
    dl_skmm = get_lcao_density_matrix(atoms.calc)
    fermi = atoms.calc.get_fermi_level()
    wfs = atoms.calc.wfs
    hl_spkmm, sl_pkmm, dl_spkmm,  \
    hl_spkcmm, sl_pkcmm, dl_spkcmm = get_pk_hsd(2, ntk,
                                                wfs.kd.ibzk_qc,
                                                hl_skmm, sl_kmm, dl_skmm,
                                                None, wfs.dtype,
                                                direction=direction)    
    my_npk = len(wfs.kd.ibzk_qc) / ntk
    my_nspins = len(wfs.kpt_u) / ( my_npk * ntk)
    
    lead_hsd = Banded_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    lead_couple_hsd = CP_Sparse_HSD(wfs.dtype, my_nspins, my_npk)
    for pk in range(my_npk):
        lead_hsd.reset(0, pk, sl_pkmm[pk], 'S', init=True)
        lead_couple_hsd.reset(0, pk, sl_pkcmm[pk], 'S', init=True)
        for s in range(my_nspins):
            lead_hsd.reset(s, pk, hl_spkmm[s, pk], 'H', init=True)     
            lead_hsd.reset(s, pk, dl_spkmm[s, pk], 'D', init=True)
            lead_couple_hsd.reset(s, pk, hl_spkcmm[s, pk], 'H', init=True)     
            lead_couple_hsd.reset(s, pk, dl_spkcmm[s, pk], 'D', init=True)          
    lead_se = LeadSelfEnergy(lead_hsd, lead_couple_hsd)
    begin += fermi
    end += fermi
    
    ee = np.linspace(begin, end, num)
    se = []

    for e in ee:
        se.append(lead_se(e).recover())
    se = np.array(se)
    
    fd = file(filename + '_' + str(world.rank), 'w')
    cPickle.dump((se, ee), fd, 2)
    fd.close()

def sort_atoms(atoms):
    pos = atoms.positions.copy()
    ind2 = fuzzy_sort(pos[:,2])
    ind1 = fuzzy_sort(pos[:,1])
    ind0 = fuzzy_sort(pos[:,0])
    tmp = ind2 * 1e12 + ind1 * 1e6 + ind0
    indices = np.argsort(tmp)
    atoms.positions = atoms.positions[indices]
    atoms.numbers = atoms.numbers[indices]
    if 'magmoms' in atoms.arrays:
        atoms.arrays['magmoms'] = atoms.arrays['magmoms'][indices]

def fuzzy_sort(seq0, tol=1e-6):
    seq = seq0.copy()
    ind0 = np.zeros_like(seq)
    ind1 = []
    n = 0
    while len(ind1) < len(seq):
        ind = []
        am = np.argmin(seq)
        tmp = seq - seq[am]        
        for j, i in enumerate(tmp):
            if abs(i) < tol:
                ind.append(j)
        seq[ind] = 1e19
        ind0[ind] = n
        ind1 += ind
        n += 1
    return ind0

def cubicing(atoms):
    cell = atoms._cell
    positions = atoms.positions
    print('cubicing only ok to [a,0,0][a/2, b, 0],[0,0,c] type ')
    tol = 1e-6
    if abs(cell[1,0]*2 - cell[0,0]) < tol:
        print('ok, possible to get a cubic structure')
        natoms = len(positions)
        new_pos = np.empty([natoms * 2, 3])
        for pos, i in zip(positions, range(natoms)):
            new_pos[i] = pos
            new_pos[i + natoms] = pos + cell[1]
        atoms += atoms
        atoms.positions = new_pos
        sort_atoms(atoms)
        cell[1, 0] = 0
        cell[1, 1] *= 2
        atoms.set_cell(cell)
    return atoms

class P_info:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0
        self.Pxsign = 1
        self.Pysign = 1
        self.Pzsign = 1
        self.N = 0
class D_info:
    def __init__(self):
        self.xy = 0
        self.xz = 0
        self.yz = 0
        self.x2y2 = 0
        self.z2r2 = 0
        self.N = 0

def egodic(nums):
    if len(nums)==1:
        return np.array(nums)
    else:
        rows = np.product(np.arange(1, len(nums) + 1))
        cols = len(nums)
        all = np.zeros([rows, cols])
        
        for i, n in enumerate(nums):
            subrows = np.product(np.arange(1, len(nums)))            
            all[i*subrows: (i+1)*subrows, 0] = n
            left_nums = nums[:]
            left_nums.remove(n)
            all[i*subrows: (i+1)*subrows, 1:] = egodic(left_nums)
        return all
    
#PPP = egodic(range(3))

def PutP(index, X, P, T):
    if P.N == 2:
        P.x = index
    if P.N == 0:
        P.y = index
    if P.N == 1:
        P.z = index
    P.N += 1
    
    if P.N == 3:
        bs = np.array([P.x, P.y, P.z])
        #c = np.array([P.Pxsign, P.Pysign, P.Pzsign])
        #c = np.resize(c, [3, 3])
        #cf = c / c.T
        #ind = np.resize(bs, [3, 3])
        ind = get_matrix_index(bs)
        T[ind.T, ind] = X
        #T[ind.T, ind] = X * cf
        P.__init__()

#DDD = egodic(range(5))

def PutD(index, X, D, T):
    if D.N == 0:
        D.xy = index
    if D.N == 3:
        D.xz = index
    if D.N == 1:
        D.yz = index
    if D.N == 4:
        D.x2y2 = index
    if D.N == 2:
        D.z2r2 = index
        
    D.N += 1
    if D.N == 5:
        sqrt = np.sqrt
        Dxy = np.array([[0, 1, 0],
                        [1, 0, 0],
                        [0, 0, 0]])
        D2xy = np.dot(X, Dxy)
        D2xy = np.dot(D2xy, X.T)
        
        Dxz = np.array([[0, 0, 1],
                        [0, 0, 0],
                        [1, 0, 0]])
        D2xz = np.dot(X, Dxz)
        D2xz = np.dot(D2xz, X.T)
        
        Dyz = np.array([[0, 0, 0],
                        [0, 0, 1],
                        [0, 1, 0]])
        D2yz = np.dot(X, Dyz)
        D2yz = np.dot(D2yz, X.T)

        Dx2y2 = np.array([[1, 0 , 0],
                          [0, -1, 0],
                          [0, 0,  0]])
        D2x2y2 = np.dot(X, Dx2y2)
        D2x2y2 = np.dot(D2x2y2, X.T)
        
        Dz2r2 = np.array([[-1, 0, 0],
                          [0, -1, 0],
                          [0,  0, 2]]) / sqrt(3)
        D2z2r2 = np.dot(X, Dz2r2)
        D2z2r2 = np.dot(D2z2r2, X.T)
        
        T[D.xy, D.xy] = D2xy[0, 1]               
        T[D.xz, D.xy] = D2xy[0, 2]               
        T[D.yz, D.xy] = D2xy[1, 2]               
        T[D.x2y2, D.xy] = (D2xy[0, 0] - D2xy[1, 1]) / 2 
        T[D.z2r2, D.xy] = sqrt(3) / 2 * D2xy[2, 2]     

        T[D.xy, D.xz] = D2xz[0, 1]               
        T[D.xz, D.xz] = D2xz[0, 2]               
        T[D.yz, D.xz] = D2xz[1, 2]               
        T[D.x2y2, D.xz] = (D2xz[0, 0] - D2xz[1, 1]) / 2 
        T[D.z2r2, D.xz] = sqrt(3) / 2 * D2xz[2,2];     

        T[D.xy , D.yz] = D2yz[0, 1]               
        T[D.xz , D.yz] = D2yz[0, 2]               
        T[D.yz , D.yz] = D2yz[1, 2]               
        T[D.x2y2, D.yz] = (D2yz[0, 0] - D2yz[1, 1]) / 2 
        T[D.z2r2, D.yz] = sqrt(3) / 2 * D2yz[2, 2]     

        T[D.xy , D.x2y2] = D2x2y2[0, 1]               
        T[D.xz , D.x2y2] = D2x2y2[0, 2]               
        T[D.yz , D.x2y2] = D2x2y2[1, 2]               
        T[D.x2y2, D.x2y2] = (D2x2y2[0, 0] - D2x2y2[1, 1]) / 2 
        T[D.z2r2, D.x2y2] = sqrt(3) / 2 * D2x2y2[2, 2]     

        T[D.xy, D.z2r2] = D2z2r2[0, 1]               
        T[D.xz, D.z2r2] = D2z2r2[0, 2]               
        T[D.yz, D.z2r2] = D2z2r2[1, 2]               
        T[D.x2y2, D.z2r2] = (D2z2r2[0, 0] - D2z2r2[1, 1]) / 2 
        T[D.z2r2, D.z2r2] = sqrt(3) / 2 * D2z2r2[2, 2]     
        
        D.__init__()      
        
def orbital_matrix_rotate_transformation(X, orbital_indices):
    nb = orbital_indices.shape[0]
    assert len(X) == 3 
    T = np.zeros([nb, nb])
    P = P_info()
    D = D_info()
    for i in range(nb):
        if orbital_indices[i, 1] == 0:
            T[i, i] = 1
        elif orbital_indices[i, 1] == 1:
            PutP(i, X, P, T)
        elif orbital_indices[i, 1] == 2:
            PutD(i, X, D, T)
        else:
            raise NotImplementError('undown shell name')
    return T

def normalize(r):
    return r/np.sqrt(np.sum(r*r))

def vector_to_paramid(r):
    r = normalize(r)
    x, y, z = r
    if z!=0:
        a1, b1, c1 = 0, 1, -y/z
        a2, b2, c2 = 1, -x*y/(y**2 + z**2), -x*z/(y**2 + z**2)
    elif y!=0:
        a1, b1, c1 = 1, -x/y, 0
        a2, b2, c2 = -x*z/(y**2 + x**2), -y*z/(y**2 + x**2), 1
    elif x!=0:
        a1, b1, c1 = -z/x, 0, 1
        a2, b2, c2 = -y*z/(x**2 + z**2), 1, -x*y/(x**2 + z**2)
    else:
        raise RuntimeError('The input vector is zero!')
    r1 = np.array([a1, b1, c1])
    r2 = np.array([a2, b2, c2])
    r1 = normalize(r1)
    r2 = normalize(r2)
    R1 = r + r1
    R2 = r - r1 / 2. + r2 / 2.
    R3 = r - r1 / 2. - r2 / 2. 
    return np.array([R1, R2, R3])
  
def transform_3d(rs1, rs2):
    assert rs1.shape == rs2.shape
    if rs1.shape[0] == 2:
        r1 = rs1[1] - rs1[0]
        r2 = rs2[1] - rs2[0]
        RS1 = vector_to_paramid(r1)
        RS2 = vector_to_paramid(r2)
    elif rs1.shape[0] == 4:
        RS1 = rs1[1:] - rs1[0]
        RS2 = rs2[1:] - rs2[0]
    else:
        raise RuntimeError('Transform atoms indices wrong!')
    X = np.dot(RS2.T, np.linalg.inv(RS1.T))
    return X

def interpolate_2d(mat):
    from gpaw.grid_descriptor import GridDescriptor
    from gpaw.transformers import Transformer
    nn = 10
    N_c = np.zeros([3], dtype=int)
    N_c[1:] = mat.shape[:2]
    N_c[0] = nn
    bmat = np.resize(mat, N_c)
    gd = GridDescriptor(N_c, N_c)
    finegd = GridDescriptor(N_c * 2, N_c)
    interpolator = Transformer(gd, finegd, 3)
    fine_bmat = finegd.zeros()
    interpolator.apply(bmat, fine_bmat)
    return fine_bmat[0]
    
def interpolate_array(array, gd, h, di='+'):
    try:
        from scipy import interpolate
        ip = True
    except ImportError:
        ip = False
    if not ip:
        return array

    dim = len(array.shape)
    assert dim == 3 or dim == 4
    spin_relate = dim == 4
    if h <= gd.h_cv[2, 2]:
        if di == '+':
            x = np.arange(gd.N_c[2]) * gd.h_cv[2, 2]
            xnew = np.arange(gd.N_c[2]) * h
        else:
            x = np.arange(-gd.N_c[2], 0) * gd.h_cv[2, 2]
            xnew = np.arange(-gd.N_c[2], 0) * h            
    else:
        if di == '+':
            x = np.arange(gd.N_c[2] * 2) * gd.h_cv[2, 2]
            xnew = np.arange(gd.N_c[2]) * h
        else:
            x = np.arange(-gd.N_c[2] * 2, 0) * gd.h_cv[2, 2]
            xnew = np.arange(-gd.N_c[2], 0) * h         
        
    if spin_relate:
        ns, nx, ny, nz = array.shape
        array.shape = (ns * nx * ny, nz)
        new_array = gd.zeros(ns, global_array=True)
        new_array.shape = (ns * nx * ny, nz)
    else:
        nx, ny, nz = array.shape
        array.shape = (nx * ny, nz)
        new_array = gd.zeros(global_array=True)
        new_array.shape = (nx * ny, nz)
      
    if h > gd.h_cv[2, 2]:
        array = np.append(array, array, 1)
        
    for i, line in enumerate(array):
        tck = interpolate.splrep(x, line, s=0)
        new_array[i] = interpolate.splev(xnew, tck, der=0)
    
    if spin_relate:
        #array.shape = (ns, nx, ny, nz)
        new_array.shape = (ns, nx, ny, nz)
    else:
        #array.shape = (nx, ny, nz)
        new_array.shape = (nx, ny, nz)
    
    return new_array
        
def eig_states_norm(orbital, s_mm):
    #normalize orbital to satisfy orbital.T.conj()*SM*orbital=unit
    norm_error = 1e-10
    ortho_error = 1e-8
    nstates = orbital.shape[1]
    d_mm = np.dot(orbital.T.conj(), s_mm)
    d_mm = np.dot(d_mm, orbital)
    for i in range(1, nstates):
        for j in range(i):
            if abs(d_mm[j ,i]) > ortho_error:
                orbital[:, i] -= orbital[:, j] * d_mm[j, i] / d_mm[j ,j]
                d_mm = np.dot(orbital.T.conj(), s_mm)
                d_mm = np.dot(d_mm, orbital)
    for i in range(nstates):
        orbital[:, i] /= np.sqrt(d_mm[i, i])

    if orbital.shape[-1] == 0:
        error = 0
    else:
        error = np.max(np.dot(np.dot(orbital.T.conj(), s_mm), orbital) -
                   np.eye(nstates)) / nstates
  
    if  abs(error) > norm_error:
        print('Warning! Normalization error %f' % error)
    return orbital

def shtm(l):
    #Spherical Harmonics transformation matrix, from the complex to real
    #The harmonics should be arranged in sequence -m, -m+1 ... m-1, m
    n = 2 * l + 1
    mtx = np.zeros([n,n], complex)
    for i in range(n):
        if i < l:
            #mtx[i, i] = 1.j*(-1.)**(l - i) / np.sqrt(2)
            #mtx[i, i + 2*(l-i)] = -1.j / np.sqrt(2)
            # The change is due to Condon-Shortley phase
            mtx[i, i] = 1.j / np.sqrt(2)
            mtx[i, i + 2*(l-i)] = -1.j*(-1.)**(l - i) / np.sqrt(2)
        elif i == l:
            mtx[i, i] = 1.
        else:
            #mtx[i, i] = 1./np.sqrt(2)
            #mtx[i, i + 2*(l-i)] = (-1.)**(i - l) / np.sqrt(2)
            mtx[i, i] = (-1.)**(i - l) /np.sqrt(2)
            mtx[i, i + 2*(l-i)] = 1. / np.sqrt(2)
    return mtx.T            

def construct_spherical_transformation_matrix(l_list):
    #construct a transformation matrix from complex harmonics to real for
    #lcao basis
    nao = np.sum([2 * l + 1 for l in l_list])
    mtx = np.zeros([nao, nao], complex)
    start = 0
    for l in l_list:
        n = 2 * l + 1
        mtx[start:start+n, start:start+n] = shtm(l)
        start += n
    return mtx

def aml(ss, l, direction):
    #calculat angular momentum matrix(complex spherical harmonics)
    #elements based on the overlap 
    amss = np.zeros(ss.shape, complex)
    n = 2 * l + 1
    for i in range(n):
        m = i - l
        # x direction=0, lx=(l+ + l-) /2
        if direction == 0:
            a1 = 0.5 * np.sqrt(l*(l+1.)-m*(m+1.))
            a2 = 0.5 * np.sqrt(l*(l+1.)-m*(m-1.))
        # y direction=1, ly=(l+ - l-) /2i
        if direction == 1:    
            a1 = -0.5 * 1.j * np.sqrt(l*(l+1.)-m*(m+1.))
            a2 = 0.5 * 1.j * np.sqrt(l*(l+1.)-m*(m-1.))
        if direction == 0 or direction == 1:  
            if m + 1 <= l:
                amss[i] += a1 * ss[i + 1]
            if m - 1 >= -l:
                amss[i] += a2 * ss[i - 1]
        elif direction == 2: #z direction=2
            amss[i] = m * ss[i]
        else: 
            raise RuntimeError('unknown direction %d' % direction)
    return amss

def angular_momentum_slice(overlap_slice, l, direction):
    #given a overlap matrix slice <l m| l_i,m> for a fixed l_i
    #return the angular momentum matrix slice <l m| l |l_i, m>
    #the slice has a shape (nao, 2*l_i + 1)
    #the magnetic number m is in the increasing sequence
    nao, n = overlap_slice.shape
    am_slice = np.zeros([nao, n], complex)
    for i in range(nao):
        ss = overlap_slice[i]
        am_slice[i] = aml(ss, l, direction)
    return am_slice     

def cut_grids_side(array, gd, gd0):
    #abstract the grid value from a including-buffer-layer calculation
    #the vaccum buffer layer is fixed on the right side
    #Assume the buffer system has the same domain spliting with the original one
    from scipy import interpolate
    global_array = gd.collect(array)
    nx, ny, nz = gd.N_c
    if gd.comm.rank == 0:
        global_array.shape = (nx * ny, nz)
    new_array = gd0.zeros()
    global_new_array = gd0.collect(new_array)
    x = np.arange(gd.N_c[2]) * gd.h_cv[2, 2]
    xnew = np.arange(gd0.N_c[2]) * gd0.h_cv[2, 2]
    nz0 = gd0.N_c[2]
    if gd0.comm.rank == 0:
        global_new_array.shape = (nx * ny, nz0)
        for i, line in enumerate(global_array):
            tck = interpolate.splrep(x, line, s=0)
            global_new_array[i] = interpolate.splev(xnew, tck, der=0)
        global_new_array.shape = (nx, ny, nz0)
    gd0.distribute(global_new_array, new_array)
    return new_array

def save_bias_data_file(Lead1, Lead2, Device):
    import pickle
    ham = Device.calc.hamiltonian
    density = Device.calc.density
    hamL = Lead1.calc.hamiltonian
    hamR = Lead2.calc.hamiltonian
    Ef = Device.calc.get_fermi_level()
    Ef_L = Lead1.calc.get_fermi_level()
    Ef_R = Lead2.calc.get_fermi_level()
    vt_sG = ham.gd.collect(ham.vt_sG) 
    vt_sG_L = hamL.gd.collect(hamL.vt_sG)
    vt_sG_R = hamR.gd.collect(hamR.vt_sG)
    vt_sG_L += (Ef - Ef_L) / Hartree
    vt_sG_R += (Ef - Ef_R) / Hartree
    vt_sG=np.append(vt_sG_L, vt_sG,axis=3)
    vt_sG=np.append(vt_sG,vt_sG_R,axis=3)
    dH_asp = collect_atomic_matrices(ham.dH_asp, ham.setups,
                                     ham.nspins, ham.gd.comm,
                                     density.rank_a)
    pickle.dump(([0.0,0,0], vt_sG, dH_asp), file('bias_data1','wb'),2)

def find(condition, flag=0):
    if flag == 1: # return an int
        return np.int(np.nonzero(condition)[0])
    else: # return an array     
        return np.nonzero(condition)[0]
   
def gather_ndarray_list(data, comm):
    #data is a numpy array, maybe has different shape in different cpus
    #this function gather them to the all_data in master, all_data is
    # a list with the lenth world.size, all_data[i] = data {on i}
    all_data = []
    dim = len(data.shape)
    shape_array = np.zeros([comm.size, dim], int)
    shape_array[comm.rank] = data.shape
    comm.sum(shape_array)
    
    if comm.rank == 0:
        all_data.append(data)
        for i in range(1, comm.size):
            tmp = np.zeros(shape_array[i], dtype=data.dtype)
            comm.receive(tmp, i, 546)
            all_data.append(tmp[:])
    else:
        comm.ssend(data, 0, 546)
    return all_data            
        
def gather_ndarray_dict(data, comm, broadcast=False):
    #data is dict of a numpy array, maybe has different shape in different cpus
    #this function gather them to the all_data in master, all_data is
    # a dict with the lenth world.size
    all_data = {}
    data_len = np.zeros([comm.size], int)
    data_len[comm.rank] = len(data)
    comm.sum(data_len)
    
    info = np.zeros([np.sum(data_len), 3], int)
    dtypes = [int, float, complex]
    for i, name in enumerate(data):
        base = np.sum(data_len[:comm.rank])
        info[base + i, 0] = len(data[name].shape)
        info[base + i, 1] = comm.rank
        info[base + i, 2] = dtypes.index(data[name].dtype)
    comm.sum(info)
    
    if comm.rank == 0:
        for name in data:
            all_data[name] = data[name]
 
        for i in range(1, comm.size):
            base = np.sum(data_len[:i])
            for j in range(data_len[i]):
                shape = np.zeros([info[base + j, 0]], int)
                dtype = dtypes[info[base + j, 2]]
                name = receive_string(i, comm)
                comm.receive(shape, i, 123)
                tmp = np.zeros(shape, dtype)
                comm.receive(tmp, i, 546)
                all_data[name] = tmp

    else:
        for name in data:
            send_string(name, 0, comm)
            shape = np.array(data[name].shape, int)
            comm.ssend(shape, 0, 123)
            comm.ssend(data[name], 0, 546)
    
    if broadcast:
        num = np.zeros([1], int)
        if comm.rank == 0:
            num[0] = len(all_data)
            comm.broadcast(num, 0)
            for name in all_data:
                broadcast_string(name, 0, comm)
                shape = np.array(all_data[name].shape, int)
                comm.broadcast(shape, 0)
                comm.broadcast(all_data[name], 0)
        else:
            comm.broadcast(num, 0)              
            for i in range(num):
                name = broadcast_string(None, 0, comm)
                shape = np.zeros([info[i, 0]], int)
                dtype = dtypes[info[i, 2]]
                comm.broadcast(shape, 0)
                tmp = np.zeros(shape, dtype)
                comm.broadcast(tmp, 0)
                all_data[name] = tmp
    return all_data
    
    
