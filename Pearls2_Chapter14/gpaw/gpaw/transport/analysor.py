from gpaw.transport.selfenergy import LeadSelfEnergy
from gpaw.transport.tools import get_matrix_index, aa1d, aa2d, sum_by_unit, \
                                dot, fermidistribution, eig_states_norm, \
                                find, get_atom_indices, dagger, \
                                write, gather_ndarray_dict
from gpaw.transport.sparse_matrix import Tp_Sparse_Matrix
from gpaw.mpi import world
from gpaw import GPAW, Mixer, MixerDif, PoissonSolver
from ase.units import Hartree, Bohr
from gpaw.utilities.memory import maxrss
import numpy as np
import copy
import cPickle
import os

class Transport_Analysor:
    def __init__(self, tp, restart=False):
        self.restart = restart
        self.data = {}
        self.ele_steps = []
        self.bias_steps = []
        self.ion_steps = []
        self.n_ele_step = 0
        self.n_bias_step = tp.n_bias_step
        self.n_ion_step = tp.n_ion_step
        self.matrix_foot_print = False
        self.reset = False
        self.scattering_states_initialized = False
        self.initialize(tp)
       
    def initialize(self, tp):
        kw = tp.analysis_parameters
        if self.restart:
            self.initialize_selfenergy_and_green_function(tp)
        else:
            self.selfenergies = tp.selfenergies
        p = self.set_default_analysis_parameters()
        for key in kw:
            if key in ['energies', 'lead_pairs', 'dos_project_atoms',
                       'project_equal_atoms', 'project_molecular_levels',
                       'isolate_atoms', 'dos_project_orbital',
                       'trans_project_orbital', 'eig_trans_channel_energies',
                        'eig_trans_channel_num', 'dos_realspace_energies']:
                p[key] = kw[key]
        for key in p:
            vars(self)[key] = p[key]       

        ef = tp.lead_fermi[0]
        path = tp.contour.get_plot_path()
        self.energies = path.energies
        self.my_energies = path.my_energies
        self.weights = path.weights
        self.my_weights = path.my_weights
        self.nids = path.nids
        self.my_nids = path.my_nids
        
        for variable in [self.eig_trans_channel_energies,
                         self.dos_realspace_energies]:
            if variable is not None:
                variable = np.array(variable) + ef
        ecomm = tp.contour.comm
        
        if self.eig_trans_channel_energies is not None:
            self.my_eig_trans_channel_energies = np.array_split(
                self.eig_trans_channel_energies, ecomm.size)[ecomm.rank]
            self.my_eig_trans_channel_nids = np.array_split(
                 np.arange(len(self.eig_trans_channel_energies)),
                 ecomm.size)[ecomm.rank]
        
        if self.dos_realspace_energies is not None:
            self.my_dos_realspace_energies = np.array_split(
                 self.dos_realspace_energies, ecomm.size)[ecomm.rank]
            self.my_dos_realspace_nids = np.array_split(
                np.arange(len(self.dos_realspace_energies)),
                ecomm.size)[ecomm.rank]

        setups = tp.inner_setups
        self.project_atoms_in_device = self.project_equal_atoms[0]
        self.project_atoms_in_molecule = self.project_equal_atoms[1]
        self.project_basis_in_device = get_atom_indices(
                                  self.project_atoms_in_device, setups)
        if self.isolate_atoms is not None:
            self.calculate_isolate_molecular_levels(tp)
        self.overhead_data_saved = False
        if world.rank == 0:
            if not os.access('analysis_data', os.F_OK):
                os.mkdir('analysis_data')
            if not os.access('analysis_data/ionic_step_0', os.F_OK):            
                os.mkdir('analysis_data/ionic_step_0')
        world.barrier()
        self.data = {}

    def save_overhead_data(self, tp):
        tp.log('-> save_overhead_data()')
        contour_parameters = {}
        cp = contour_parameters
        cp['neintmethod'] = tp.neintmethod
        cp['neintstep'] = tp.neintstep
        cp['lead_fermi'] = tp.lead_fermi
        cp['kt'] = tp.occupations.width * Hartree
        if not tp.non_sc:
            cp['eqinttol'] = tp.contour.eqinttol
        
        basis_information = {}
        bi = basis_information
        bi['orbital_indices'] = tp.orbital_indices
        bi['lead_orbital_indices'] = tp.lead_orbital_indices
        tp.log('   tp.hsd.S[0] : {0}'.format(tp.hsd.S[0]))
        bi['ll_index'] = tp.hsd.S[0].ll_index
        bi['ex_ll_index'] = tp.hsd.S[0].ex_ll_index
          
        atoms = tp.atoms.copy()
        if world.rank == 0:
            fd = file('analysis_overhead', 'wb')
            tp.log('   cPickle.dump()')
            cPickle.dump((atoms, basis_information, contour_parameters), fd, 2)
            fd.close()
        tp.log('<- save_overhead_data()')
           
    def set_default_analysis_parameters(self):
        p = {}
        p['energies'] = np.linspace(-5., 5., 201) 
        p['lead_pairs'] = [[0,1]]
        p['dos_project_atoms'] = None
        p['project_molecular_levels'] = None
        p['project_equal_atoms'] = [[], []]
        p['isolate_atoms'] = None
        p['dos_project_orbital'] = None
        p['trans_project_orbital'] = None
        p['eig_trans_channel_energies'] = None
        p['eig_trans_channel_num'] = 0
        p['dos_realspace_energies'] = None
        return p      
                
    def initialize_selfenergy_and_green_function(self, tp):
        self.selfenergies = []
        if tp.use_lead:
            for i in range(tp.lead_num):
                self.selfenergies.append(LeadSelfEnergy(tp.lead_hsd[i],
                                                     tp.lead_couple_hsd[i]))
    
                self.selfenergies[i].set_bias(tp.bias[i])
            
    def reset_selfenergy(self, tp, s, k):
        for i in range(tp.lead_num):
            self.selfenergies[i].s = s
            self.selfenergies[i].pk = k
    
    def reset_green_function(self, tp, s, k): 
        tp.hsd.s = s
        tp.hsd.pk = k
      
    def calculate_green_function_of_k_point(self, tp, s, k, energy, sigma,
                                                                 full=False):
        self.reset_green_function(tp, s, k)
        return tp.hsd.calculate_eq_green_function(energy, sigma,
                                                       ex=False, full=full)
 
    def calculate_sigma(self, tp, s, k, energy, nid_flag=None):
        sigma = []
        for i in range(tp.lead_num):
            sigma.append(self.selfenergies[i](energy, nid_flag))        
        return sigma

    def get_gamma(self, tp, sigma):
        gamma = []
        for i in range(tp.lead_num):
            gamma.append(1.j * (sigma[i].recover() -
                                                sigma[i].recover().T.conj()))
        return gamma
        
    def calculate_transmission(self, tp, s, k, energy, nid_flag=None):
        tp.log('-> calculate_transmission()')
        self.reset_selfenergy(tp, s, k)
        self.reset_green_function(tp, s, k)
        tp.log('  calculate_sigma()')
        sigma = self.calculate_sigma(tp, s, k, energy, nid_flag)
        tp.log('  calculate_gamma()')
        gamma = self.get_gamma(tp, sigma)    
        trans_coff = []
        for i, lead_pair in enumerate(self.lead_pairs):
            l1, l2 = lead_pair
            if i == 0:
                tp.log('  tp.hsd.abstract_sub_green_matrix()')
                gr_sub, inv_mat = tp.hsd.abstract_sub_green_matrix(
                                                    energy, sigma, l1, l2)
            else:
                tp.log('  tp.hsd.abstract_sub_green_matrix()')
                gr_sub = tp.hsd.abstract_sub_green_matrix(energy,
                                                    sigma, l1, l2, inv_mat)                    
            transmission =  dot(dot(gamma[l1], gr_sub),
                                              dot(gamma[l2], gr_sub.T.conj()))
            trans_coff.append(np.real(np.trace(transmission)))        
        trans_coff = np.array(trans_coff)
        tp.log('<- calculate_transmission()')
        return trans_coff
    
    def calculate_dos(self, tp, s, k, energy, nid_flag=None):
        self.reset_selfenergy(tp, s, k)
        self.reset_green_function(tp, s, k)
        sigma = self.calculate_sigma(tp, s, k, energy, nid_flag)
        gr = self.calculate_green_function_of_k_point(tp, s, k, energy, sigma)        
        dos = - np.imag(np.diag(dot(gr, tp.hsd.S[k].recover()))) / np.pi         
        return dos
    
    def calculate_eigen_transport_channel(self, tp, s, q, energy):
        t, vc, k, vl = self.central_scattering_states(tp, energy, s, q)
        weights = []
        velocities = []
        vectors = []
        for i in range(tp.lead_num):
            weights.append([])
            velocities.append([])
            vectors.append([])
            for j in range(tp.lead_num):
                zeta = t[i][j]
                zeta2 = np.dot(zeta.T.conj(), zeta)
                if np.product(zeta2.shape) != 0:
                    w, v = np.linalg.eig(zeta2)
                    weights[i].append(w)
                    velocities[i].append(v)
                    vectors[i].append(np.dot(vc[i], v))
        weights = np.array(weights)
        velocities = np.array(velocities)
        vectors = np.array(vectors)
        return weights, velocities, vectors
  
    def calculate_charge_distribution(self, tp, s, q):
        d_mm = tp.hsd.D[s][q].recover(True)
        s_mm = tp.hsd.S[q].recover(True)
        q_mm = np.dot(d_mm, s_mm)
        return np.diag(q_mm).real

    def calculate_isolate_molecular_levels(self, tp):
        atoms = self.isolate_atoms
        atoms.pbc = True
        p = tp.gpw_kwargs.copy()
        p['nbands'] = None
        if 'mixer' in p:
            if not tp.spinpol:
                p['mixer'] = Mixer(0.1, 5, weight=100.0)
            else:
                p['mixer'] = MixerDif(0.1, 5, weight=100.0)
        p['poissonsolver'] = PoissonSolver(nn=2)
        if type(p['basis']) is dict and len(p['basis']) == len(tp.atoms):
            p['basis'] = 'dzp'
            raise Warning('the dict basis is not surpported in isolate atoms')
        if 'txt' in p and p['txt'] != '-':
            p['txt'] = 'isolate_' + p['txt']
        atoms.set_calculator(GPAW(**p))
        atoms.get_potential_energy()
        setups = atoms.calc.wfs.setups
        self.project_basis_in_molecule = get_atom_indices(
                                      self.project_atoms_in_molecule, setups)         
        kpt = atoms.calc.wfs.kpt_u[0]
        s_mm = atoms.calc.wfs.S_qMM[0]
        c_nm = eig_states_norm(kpt.C_nM, s_mm)
        self.isolate_eigen_values = kpt.eps_n
        self.isolate_eigen_vectors = c_nm
        self.isolate_s_mm = s_mm

    def calculate_project_transmission(self, tp, s, q, energy):
        eps_n, c_nm, s_mm = self.isolate_eigen_values, \
                                 self.isolate_eigen_vectors, self.isolate_s_mm
        nl = tp.lead_num
        T0 = np.zeros(c_nm.shape[0])
        ind1 = self.project_basis_in_molecule
        ind2 = self.project_basis_in_device
        t, vc, k, vl = self.central_scattering_states(tp, energy, s, q)
        project_transmission = []
        for j in range(nl):
            project_transmission.append([])
            for k in range(nl):
                vs = vc[k][ind2]
                vm = np.dot(np.dot(c_nm.T.conj(), s_mm)[:, ind1], vs)
                t0 = t[j][k]
                if len(t0) > 0:
                    pt = vm * vm.conj() * \
                             np.diag(np.dot(t0.T.conj(), t0)) \
                                  / np.diag(np.dot(vm.T.conj(), vm))
                else:
                    pt = T0
                project_transmission[j].append(pt)
        project_transmission = np.array(project_transmission)
        return project_transmission

    def calculate_realspace_wave_functions(self, tp, C_nm, q):
        #nl number of molecular levels
        wfs = tp.wfs
        nao = wfs.setups.nao
        nb, nl = C_nm.shape
        if wfs.dtype == float:
            C_nm = C_nm.real.copy()
        #extended_C_nm = np.zeros([nl, nao], wfs.dtype)
        total_psi_g = []
        for i in range(nl):
            psi_g = tp.gd.zeros(nl, dtype=wfs.dtype)
            c_nm = C_nm.reshape(1, -1)
            psi_g = psi_g.reshape(1, -1)
            wfs.basis_functions.lcao_to_grid(c_nm, psi_g, q)
            psi_g.shape = tp.gd.n_c
            global_psi_g = tp.gd.collect(psi_g / Bohr**1.5)
            total_psi_g.append(global_psi_g)
        return np.array(total_psi_g)

    def get_left_channels(self, tp, energy, s, k):
        # to get the left scattering channel from lead to scattering region
        sigma = self.calculate_sigma(tp, s, k, energy)
        g_s_ii = self.calculate_green_function_of_k_point(tp, s, k, energy,
                                                          sigma, full=True)
        nb = g_s_ii.shape[-1]
        dtype = g_s_ii.dtype
        lambda_l_ii = np.zeros([nb, nb], dtype)
        lambda_r_ii = np.zeros([nb, nb], dtype)
        ind = get_matrix_index(tp.hsd.S[0].ll_index[0][-1])
        lambda_l_ii[ind.T, ind] = 1.j * (sigma[0].recover() -
                                                  sigma[0].recover().T.conj())
        ind = get_matrix_index(tp.hsd.S[0].ll_index[1][-1])        
        lambda_r_ii[ind.T, ind] = 1.j * (sigma[1].recover() -
                                                  sigma[1].recover().T.conj())
        s_mm = tp.hsd.S[k].recover()
        s_s_i, s_s_ii = np.linalg.eig(s_mm)
        s_s_i = np.abs(s_s_i)
        s_s_sqrt_i = np.sqrt(s_s_i) # sqrt of eigenvalues  
        s_s_sqrt_ii = np.dot(s_s_ii * s_s_sqrt_i, dagger(s_s_ii))
        s_s_isqrt_ii = np.dot(s_s_ii / s_s_sqrt_i, dagger(s_s_ii))

        lambdab_r_ii = np.dot(np.dot(s_s_isqrt_ii, lambda_r_ii),s_s_isqrt_ii)
        a_l_ii = np.dot(np.dot(g_s_ii, lambda_l_ii), dagger(g_s_ii))
        ab_l_ii = np.dot(np.dot(s_s_sqrt_ii, a_l_ii), s_s_sqrt_ii)
        lambda_i, u_ii = np.linalg.eig(ab_l_ii)
        ut_ii = np.sqrt(lambda_i / (2.0 * np.pi)) * u_ii
        m_ii = 2 * np.pi * np.dot(np.dot(dagger(ut_ii), lambdab_r_ii),ut_ii)
        T_i,c_in = np.linalg.eig(m_ii)
        T_i = np.abs(T_i)
        channels = np.argsort(-T_i)
        c_in = np.take(c_in, channels, axis=1)
        T_n = np.take(T_i, channels)
        v_in = np.dot(np.dot(s_s_isqrt_ii, ut_ii), c_in)
        return T_n, v_in
    
    def calculate_realspace_dos(self, tp, energy):
        wfs = tp.extended_calc.wfs
        dosg = wfs.gd.zeros(tp.nspins)
        for kpt in wfs.kpt_u:
            s = kpt.s
            q = kpt.q
            sigma = self.calculate_sigma(tp, s, q, energy)

            self.reset_green_function(tp, s, q)
            tp.hsd.G.reset_from_others(tp.hsd.S[q],
                                       tp.hsd.H[s][q],
                                       energy, -1, init=True)
            
            for i in range(tp.hsd.G.lead_num):
                tp.hsd.G.diag_h[i][-1].reset_minus(sigma[i])
            
            gr = np.linalg.inv(tp.hsd.G.recover(True))
            dos_mm = np.dot(gr, tp.hsd.S[q].recover(True))
        
            if wfs.dtype == float:
                dos_mm = np.real(dos_mm).copy()
            wfs.basis_functions.construct_density(dos_mm, dosg[kpt.s], kpt.q)
        wfs.kd.comm.sum(dosg)
        wfs.band_comm.sum(dosg)
        global_dosg = wfs.gd.collect(dosg)        
        return global_dosg
       
    def lead_k_matrix(self, tp, l, s, pk, k_vec, hors='S'):
        if hors == 'S':
            h00 = tp.lead_hsd[l].S[pk].recover()
            h01 = tp.lead_couple_hsd[l].S[pk].recover()
        else:
            h00 = tp.lead_hsd[l].H[s][pk].recover()
            h01 = tp.lead_couple_hsd[l].H[s][pk].recover()
        h01 *=  np.exp(2 * np.pi * 1.j * k_vec)
        return h00 + h01 + h01.T.conj()
      
    def lead_scattering_states(self, tp, energy, l, s, q):
        #Calculating the scattering states in electrodes
        #l ---- index of electrode
        #s ---- index of spin
        #q ---- index of local k point
        
        #if it is multi-terminal system, should add a part that can rotate
        # the lead hamiltonian
        MaxLambda = 1e2
        MinErr = 1e-8
        energy += tp.bias[l]
        hes00 = tp.lead_hsd[l].H[s][q].recover() - \
                                    tp.lead_hsd[l].S[q].recover() * energy
        hes01 = tp.lead_couple_hsd[l].H[s][q].recover() - \
                             tp.lead_couple_hsd[l].S[q].recover() * energy
        nb = hes00.shape[-1]
        dtype = hes00.dtype
        A = np.zeros([2*nb, 2*nb], dtype)
        B = np.zeros([2*nb, 2*nb], dtype)
            
        A[:nb, nb:2*nb] = np.eye(nb)
        A[nb:2*nb, :nb] = hes01.T.conj()
        A[nb:2*nb, nb:2*nb] = hes00
            
        B[:nb, :nb] = np.eye(nb)
        B[nb:2*nb, nb:2*nb] = -hes01
            
        from scipy.linalg import eig
        D, V = eig(A, B)
        index = np.argsort(abs(D))
        D = D[index]
        V = V[:, index]
        
        #delete NaN
        index = find(np.abs(D) >= 0)
        D = D[index]
        V = V[:, index]
        
        #delete some unreasonable solutions
        index = find(abs(D) > MaxLambda)
        cutlen = len(D) - index[0]
        index = np.arange(cutlen, len(D) - cutlen)
        D = D[index]
        V = V[:, index]
            
        k = np.log(D) * (-1.j) / 2 / np.pi
        Vk = V[:nb]
            
        #sort scattering states 
        proindex = find(abs(k.imag) < MinErr)
            
        if len(proindex) > 0:
            k_sort = np.sort(k[proindex].real)
            index = np.argsort(k[proindex].real)
            k[proindex] = k[proindex[index]].real
            Vk[:, proindex] = Vk[:, proindex[index]]
                
            #normalization the scattering states
            j = proindex[0]
            while j <= proindex[-1]:
                same_k_index = find(abs((k - k[j]).real) < MinErr)
                sk = self.lead_k_matrix(tp, l, s, q, k[j])
                Vk[:, same_k_index] = eig_states_norm(Vk[:, same_k_index],
                                                                         sk)
                j += len(same_k_index)
        return np.array(k[proindex]), np.array(Vk[:, proindex])
                   
    def central_scattering_states(self, tp, energy, s, q):
        #To get the scattering states corresponding to a Bloch vector in lead
        MaxLambda = 1e2
        MinErr = 1e-8
        bc = tp.inner_mol_index
        nc = len(bc)
        molhes = tp.hsd.H[s][q].recover(True) - \
                              tp.hsd.S[q].recover(True) * energy
        blead = []
        lead_hes = []
        lead_couple_hes = []
        for i in range(tp.lead_num):
            blead.append(tp.lead_layer_index[i][-1])
            lead_hes.append(tp.lead_hsd[i].H[s][q].recover() -
                                       tp.lead_hsd[i].S[q].recover() * energy)
            lead_couple_hes.append(tp.lead_couple_hsd[i].H[s][q].recover() -
                                tp.lead_couple_hsd[i].S[q].recover() * energy)

        ex_ll_index = tp.hsd.S[0].ex_ll_index    
        
        total_k = []
        total_vk = []
        total_lam = []
        total_v = []
        total_pro_right_index = []
        total_pro_left_index = []
        total_left_index = []
        total_right_index = []
        total_kr = []
        total_kt = []
        total_lambdar = []
        total_lambdat = []
        total_vkr = []
        total_vkt = []
        total_vr = []
        total_vt = []
        total_len_k = []
        total_nblead = []
        total_bA1 = []
        total_bA2 = []
        total_bproin = []
        total_nbB2 = []
        total_bB2 = []
                
        for i in range(tp.lead_num):
            k, vk = self.lead_scattering_states(tp, energy, i, s, q)
            total_k.append(k)
            total_vk.append(vk)
            
            lam = np.exp(2 * np.pi * k * 1.j)
            total_lam.append(lam)
            
            #calculating v = dE/dk
            de2 = 1e-8
            k2, vk2 = self.lead_scattering_states(tp, energy + de2, i, s, q)
            v = de2 / (k2 - k) / 2 / np.pi 
            total_v.append(v)
            
            #seperating left scaterring states and right scattering states
            #left scattering: lead->mol, right scattering: mol->lead
            proindex = find(abs(k.imag) < MinErr)
            pro_left_index = proindex[find(v[proindex].real < 0)]
            pro_left_index = pro_left_index[-np.arange(len(pro_left_index))
                                                                          - 1]
            
            pro_right_index = proindex[find(v[proindex].real > 0)]
               
            left_index = find(k.imag > MinErr)
            left_index = np.append(left_index, pro_left_index)
            
            right_index = pro_right_index.copy()
            right_index = np.append(right_index, find(k.imag < -MinErr))
            
            total_pro_left_index.append(pro_left_index)
            total_pro_right_index.append(pro_right_index)
            total_left_index.append(left_index)
            total_right_index.append(right_index)
           
            kr = k[left_index]
            kt = k[right_index]
            total_kr.append(kr)
            total_kt.append(kt)
            
            lambdar = np.diag(np.exp(2 * np.pi * kr * 1.j))
            lambdat = np.diag(np.exp(2 * np.pi * kt * 1.j))
            total_lambdar.append(lambdar)
            total_lambdat.append(lambdat)
            
            vkr = np.take(vk, left_index, axis=1)
            vkt = np.take(vk, right_index, axis=1)
            vr = v[pro_left_index]
            vt = v[pro_right_index]
            total_vkr.append(vkr)
            total_vkt.append(vkt)
            total_vr.append(vr)
            total_vt.append(vt)
            
            #abstract basis information
            len_k = len(right_index)
            total_len_k.append(len_k)
            
            #lead i basis seqeunce in whole matrix
            nblead = len(ex_ll_index[i][-1])
            total_nblead.append(nblead)
            bA1 = nc + int(np.sum(total_nblead[:i])) + np.arange(nblead)
            #sequence of the incident wave   
            bA2 = nc + int(np.sum(total_len_k[:i])) + np.arange(len_k)
            # the first n in vkt are scattering waves, the rest are the decaying ones
            # the first n in vkr are decaying waves, the rest are the scattering ones
            total_bA1.append(bA1)
            total_bA2.append(bA2)  
           
            bproin = np.arange(len(pro_right_index)) + \
                                          len(kr) - len(pro_right_index)
            ### this line need to check it...
            total_bproin.append(bproin)
            nbB2 = len(bproin)
            total_nbB2.append(nbB2)
            bB2 = int(np.sum(total_nbB2[:i])) + np.arange(nbB2)
            total_bB2.append(bB2)

        ind = get_matrix_index(bc)              
        Acc = molhes[ind.T, ind]
           
        total_Alc = []
        total_hcl = []
        total_All = []
        total_Acl = []
        total_Bc = []
        total_Bl = []    
            
        for i in range(tp.lead_num):
            ind1, ind2 = get_matrix_index(blead[i], bc)         
            Alc = molhes[ind1, ind2]
            ind1, ind2 = get_matrix_index(bc, blead[i])         
            hcl = molhes[ind1, ind2]
            vkt = total_vkt[i]
            All = np.dot(tp.lead_hsd[i].H[s][q].recover(), vkt) + \
                    np.dot(np.dot(tp.lead_couple_hsd[i].H[s][q].recover(), \
                                                     vkt), total_lambdat[i])
            Acl = np.dot(hcl, vkt)
            bpi = total_bproin[i]
            vkr = total_vkr[i]
            vkr_bpi = np.take(vkr, bpi, axis=1)
            Bc = -np.dot(hcl, vkr_bpi)
            ind = get_matrix_index(bpi)
            Bl = -np.dot(tp.lead_hsd[i].H[s][q].recover(), vkr_bpi) + \
                    np.dot(np.dot(tp.lead_couple_hsd[i].H[s][q].recover(),
                        vkr_bpi) ,total_lambdar[i][ind.T, ind])
            total_Alc.append(Alc)
            total_Acl.append(Acl)
            total_hcl.append(hcl)
            total_All.append(All)
            total_Bc.append(Bc)
            total_Bl.append(Bl)
        total_bc = molhes.shape[-1]
        MatA = np.zeros([total_bc, nc + np.sum(total_len_k)], complex)    
        MatB = np.zeros([total_bc, np.sum(total_nbB2)], complex)
        ind = get_matrix_index(np.arange(nc))
        MatA[ind.T, ind] = Acc
        for i in range(tp.lead_num):
            ind1, ind2 = get_matrix_index(total_bA1[i], total_bA2[i])
            MatA[ind1, ind2] = total_All[i]
            ind1, ind2 = get_matrix_index(np.arange(nc), total_bA2[i])
            MatA[ind1, ind2] = total_Acl[i]
            ind1, ind2 = get_matrix_index(total_bA1[i], np.arange(nc))
            MatA[ind1, ind2] = total_Alc[i]
            ind1, ind2 = get_matrix_index(np.arange(nc), total_bB2[i])
            MatB[ind1, ind2] = total_Bc[i]
            ind1, ind2 = get_matrix_index(total_bA1[i], total_bB2[i])
            MatB[ind1, ind2] = total_Bl[i]
        Vx, residues, rank, singular = np.linalg.lstsq(MatA, MatB)
        
        total_vl = []
        for i in range(tp.lead_num):
            total_k[i] = total_k[i][total_pro_right_index[i]]
            total_vl.append(np.take(total_vk[i],
                                    total_pro_right_index[i], axis=1))
         
        total_vc = []
        t = []
        for i in range(tp.lead_num):
            t.append([])
            for j in range(tp.lead_num):
                t[i].append([])
        for i in range(tp.lead_num):
            ind1, ind2 = get_matrix_index(np.arange(nc), total_bB2[i])
            total_vc.append(Vx[ind1, ind2])
            for j in range(tp.lead_num):
                bx = total_bA2[j]
                bx = bx[:len(total_pro_right_index[j])]
                ind1, ind2 = get_matrix_index(bx, total_bB2[i])
                #t[j].append(Vx[ind1, ind2])
                t[j][i] = Vx[ind1, ind2]
                for m in range(len(total_pro_left_index[i])):
                    for n in range(len(total_pro_right_index[j])):
                        t[j][i][n, m] *= np.sqrt(abs(total_vt[i][n]/
                                                          total_vr[j][m]))
        return np.array(t), np.array(total_vc), \
                                         np.array(total_k), np.array(total_vl)
            
    def cal_sstates(self, tp, energies, s, q):
        #To see how much the transmision from one orbital(local basis) in lead a
        # to the orbital in the lead b, the transmission is decomposited to
        # Bloch wave in leads
        MinErr = 1e-8
        MaxLambda = 1 + MinErr
        
        energies = np.real(energies)
        ne = len(energies)
        t_all = np.zeros([ne, tp.lead_num, tp.lead_num])
        t_lead_all = []
        for i in range(tp.lead_num):
            t_lead_all.append([])
            for j in range(tp.lead_num):
                t_lead_all[i].append(np.zeros([tp.nblead[i],
                                                     tp.nblead[j], ne]))
        total_t = []
        total_vc = []
        total_k = []
        total_vl = []
        sk1 = []
        sk2 = []
        for n, energy in enumerate(energies):
            t, vc, k, vl = self.central_scattering_states(tp, energy, s, q)
            total_t.append(t)
            total_vc.append(vc)
            total_k.append(k)
            total_vl.append(vl)
            for i in range(tp.lead_num):
                for j in range(tp.lead_num):
                    t0 = t[i][j]
                    t_all[n, i, j] = np.sum(np.dot(t0.T.conj(), t0))
                    for m1 in range(t0.shape[0]):
                        sk1.append(self.lead_k_matrix(tp, i, s, q, k[i][m1], 'S'))
                    for m2 in range(t0.shape[1]):
                        sk2.append(self.lead_k_matrix(tp, i, s, q, k[j][m2], 'S'))
                    for m1 in range(t0.shape[0]):
                        for m2 in range(t0.shape[1]):
                            w1 = (vl[j][:, m1].conj() * np.dot(sk1[m1],
                                                        vl[j][:, m1])).real
                            w1[find(w1 < 0)] = 0
                            w1 /= np.sum(w1)
                            
                            w2 = (vl[i][:, m2].conj() * np.dot(sk2[m2],
                                                        vl[i][:, m2])).real
                            w2[find(w2 < 0)] = 0                        
                            w2 /= np.sum(w2)
                            
                            t_lead_all[i][j][:, :, n] += abs(t0[m1,
                                          m2]) ** 2 * np.dot(w2, w1.T.conj())
        return np.array(t_lead_all), np.array(t_all)

    
    def save_ele_step(self, tp):
        if not self.matrix_foot_print:
            fd = file('matrix_sample', 'wb')
            sample = Tp_Sparse_Matrix(complex, tp.hsd.ll_index)
            sample.reset_from_others(tp.hsd.S[0], tp.hsd.H[0][0], 1.,
                                                      -1.j, init=True)
            cPickle.dump(sample, fd, 2)
            fd.close()
            self.matrix_foot_print = True
        #selfconsistent calculation: extended_calc ---> extended_atoms
        #non_selfconsistent calculation: extended_calc ---> original_atoms        
        
        calc = tp.extended_calc
        gd = calc.gd
        finegd = calc.hamiltonian.finegd

        if not tp.use_qzk_boundary:
            nt_sG = tp.density.nt_sG
            nt_sg = tp.density.nt_sg
        else:
            nt_sG = calc.density.nt_sG            
            nt_sg = calc.density.nt_sg
        vt_sG = calc.hamiltonian.vt_sG
        vt_sg = calc.hamiltonian.vt_sg
        data = self.data
        flag = 'ele_' + str(self.n_ele_step) + '_'
        nt = []
        vt = []
        ntg = []
        vtg = []
        for s in range(tp.nspins): 
            nts = aa1d(nt_sG[s]) 
            vts = aa1d(vt_sG[s]) * Hartree
            ntgs = aa1d(nt_sg[s]) 
            vtgs = aa1d(vt_sg[s]) * Hartree
            nt.append(nts)
            vt.append(vts)
            ntg.append(ntgs)
            vtg.append(vtgs)

        data[flag + 'nt'] = np.array(nt)
        data[flag + 'vt'] = np.array(vt)
        data[flag + 'ntg'] = np.array(ntg)
        data[flag + 'vtg'] = np.array(vtg)
        data[flag + 'vbarg'] = aa1d(calc.hamiltonian.vbar_g)
        data[flag + 'df'] = np.diag(tp.hsd.H[0][0].recover(True))
        data[flag + 'dd'] = np.diag(tp.hsd.D[0][0].recover(True))            
        
        if not tp.use_qzk_boundary: 
            gd = tp.finegd
            rhot_g = tp.density.rhot_g
            rho = aa1d(rhot_g)
            data[flag + 'rho'] = np.array(rho)
        else:
            gd = calc.finegd
            rhot_g = calc.density.rhot_g
            rho = aa1d(rhot_g)                
            data[flag + 'rho'] = np.array(rho)
       
        gd = finegd
        vHt_g = calc.hamiltonian.vHt_g
        vHt = aa1d(vHt_g) * Hartree
        data[flag + 'vHt'] = vHt
        self.n_ele_step += 1
      
    def save_bias_step(self, tp):
        tp.log('-> save_bias_step()')
        tp.text('----- saving bias step {0} ------'.format(self.n_bias_step))
        if not self.overhead_data_saved:
            tp.log('   save_overhead_data()')
            self.save_overhead_data(tp)
            self.overhead_data_saved = True
        ks_map = tp.my_ks_map

        if 'tc' in tp.analysis_data_list:
            tp.log('   collect_transmission_and_dos()')
            tc, dos = self.collect_transmission_and_dos(tp)
            if not tp.non_sc:
                tp.log('   calculate_current()')
                current = self.calculate_current(tp, tc)
            else:
                current = np.array(0)
            self.data['tc'] = tc
            self.data['dos'] = dos
            self.data['current'] = current
        tp.log('   abstract_d_and_v()')
        nt, vt, ntx, vtx, nty, vty = self.abstract_d_and_v(tp)
        for name in ['nt', 'vt', 'ntx', 'vtx', 'nty', 'vty']:
                self.data[name] = eval(name)
        if tp.non_sc or tp.analysis_mode:
            force = None
            contour = None
        else:
            if 'force' in tp.analysis_data_list:
                if not tp.use_qzk_boundary:
                    force = tp.calculate_force() * Hartree / Bohr
                else:
                    force = tp.extended_calc.get_forces(tp.extended_atoms
                                                    )[:len(tp.atoms)]
            else:
                force = None
            tp.F_av = None
            tp.log('   collect_contour()')
            contour = self.collect_contour(tp)            
    
        tp.log('   collect_charge()')
        charge = self.collect_charge(tp)
        magmom = tp.occupations.magmom
        local_magmom = tp.get_magnetic_moments()   

        lead_fermi = np.array(tp.lead_fermi)
        lead_pairs = np.array(self.lead_pairs)
        bias = np.array(tp.bias)
        gate = np.array(tp.gate)
        magmom = np.array(magmom)
        local_magmom = np.array(local_magmom)
        for name in ['lead_fermi', 'lead_pairs', 'bias', 'gate', 
                                      'charge', 'magmom', 'local_magmom']:
            self.data[name] = eval(name)
        # do not include contour now because it is a dict, not a array able to
        # collect, but will do it at last
        #self.data = gather_ndarray_dict(self.data, tp.contour.comm)
        self.data['contour'] = contour
        self.data['force'] = force
        self.data['kpt_rank'] = tp.wfs.kd.comm.rank
        self.data['kpt_size'] = tp.wfs.kd.comm.size
        self.data['domain_rank'] = tp.gd.comm.rank
        self.data['domain_parpos'] = tp.gd.parpos_c
        self.data['domain_parsize'] = tp.gd.parsize_c
        ns, npk, nlp, ne = tp.nspins, tp.npk, len(self.lead_pairs), len(self.energies)
        self.data['transmission_dimension'] = (ns, npk, nlp, ne)
        nb = tp.nbmol_inner
        self.data['dos_dimension'] = (ns, npk, nb, ne)

        if tp.non_sc:
            self.data['total_energy'] = tp.guess_total_energy
        for condition, obj, name in tp.special_datas:
            if eval(condition):
                self.data[name] = eval(obj)
       
        tp.log('   tio.save_data()')
        tp.tio.save_data(self, option='Analysis', bias_step=self.n_bias_step)
        tp.log('   tio.arrange_analysis_data()')
        tp.tio.arrange_analysis_data(self.n_bias_step, 
                                            self.n_ion_step, tp.analysis_mode)
        #if world.rank == 0:
        #    if tp.analysis_mode:
        #        filename = '/abias_step_' + str(self.n_bias_step)
        #    else:
        #        filename = '/bias_step_' + str(self.n_bias_step)
        #    fd = file('analysis_data/ionic_step_' + str(self.n_ion_step)
        #              + filename, 'wb')
        #    cPickle.dump(self.data, fd, 2)
        #    fd.close()
        self.data = {}
        self.n_ele_step = 0
        self.n_bias_step += 1
        tp.log('<- save_bias_step()')
        

    def collect_transmission_and_dos(self, tp, energies=None, nids=None):
        tp.log('-> collect_transmission_and_dos()')
        if energies is None:
            energies = self.my_energies
        if nids is None:
            nids = self.my_nids
      
        nlp = len(self.lead_pairs)
        ne = len(energies)
        nbmol = tp.nbmol_inner
        ns, npk = tp.my_nspins, tp.my_npk
        local_tc_array = np.empty([ns, npk, nlp, ne], float)
        local_dos_array = np.empty([ns, npk, nbmol, ne], float)
        
        tp.log('  nr. energies: {0}'.format(len(energies)))       
        tp.log('  nids: {0}'.format(nids))
        for s in range(ns):
            for q in range(npk):
                for e, energy, nid in zip(range(len(energies)), energies, nids):
                    tp.log('  calculate_dos() {0} {1} {2}'.format(s,q,e))
                    local_dos_array[s, q, :, e] = self.calculate_dos(tp, s, q,
                                                               energy, nid)
                    tp.log('  calculate_transmission() {0} {1} {2}'.format(s,q,e))
                    local_tc_array[s, q, :, e] =  self.calculate_transmission(tp, s,
                                                            q, energy, nid)

        #kpt_comm = tp.wfs.kd.comm
        #ns, npk = tp.nspins, tp.npk
        #if kpt_comm.rank == 0:
        #    tc_array = np.empty([ns, npk, nlp, ne], float)
        #    dos_array = np.empty([ns, npk, nbmol, ne], float)
        #else:
        #    tc_array = None
        #    dos_array = None
        #kpt_comm.gather(local_tc_array, 0, tc_array)
        #kpt_comm.gather(local_dos_array, 0, dos_array)                    
        tp.log('<- collect_transmission_and_dos()')
        return local_tc_array, local_dos_array
    
    def collect_charge(self, tp):
        nbmol = tp.nbmol + np.sum(tp.nblead)
        ns, npk = tp.my_nspins, tp.my_npk
        local_charge_array = np.zeros([ns, npk, nbmol])
        for s in range(ns):
            for q in range(npk):
                    local_charge_array[s, q] = \
                                     self.calculate_charge_distribution(tp, s, q)
        charge_array = np.sum(local_charge_array, axis=1)
        return charge_array

    def collect_contour(self, tp):
        my_eq_contour = {}
        my_ne_contour = {}
        my_loc_contour = {}
        num = 0
        for s in range(tp.my_nspins):
            for q in range(tp.my_npk):
                flag = str(tp.my_ks_map[num, 0]) + str(tp.my_ks_map[num, 1])
                my_eq_contour[flag] = np.array(tp.eqpathinfo[s][q].energy)
                my_ne_contour[flag] = np.array(tp.nepathinfo[s][q].energy)
                if not tp.ground:
                    my_loc_contour[flag] = np.array(tp.locpathinfo[s][q].energy)        
        #eq_contour = gather_ndarray_dict(my_eq_contour, tp.wfs.kd.comm)
        #ne_contour = gather_ndarray_dict(my_ne_contour, tp.wfs.kd.comm)        
        #if not tp.ground:
        #    loc_contour = gather_ndarray_dict(my_loc_contour, tp.wfs.kd.comm)
        #else:
        #    loc_contour = None
        contour = {'eq': my_eq_contour, 'ne': my_ne_contour, 'loc': my_loc_contour}
        return contour
 
    def save_ion_step(self, tp):
        if world.rank == 0:
            fd = file('analysis_data/ionic_step_' +
                  str(self.n_ion_step) +'/positions', 'wb')
            cPickle.dump(tp.atoms.positions, fd, 2)
            fd.close()
        self.n_bias_step = 0
        self.n_ion_step += 1
        if world.rank == 0:
            dirname = 'analysis_data/ionic_step_' + str(self.n_ion_step)
            if not os.access(dirname, os.F_OK):            
                os.mkdir(dirname)
        world.barrier()
 
    def abstract_d_and_v(self, tp):
        calc = tp.extended_calc
        gd = calc.gd
        if not tp.use_qzk_boundary:
            nt_sG = tp.density.nt_sG
        else:
            nt_sG = calc.density.nt_sG            
        vt_sG = calc.hamiltonian.vt_sG
        nt = []
        vt = []
        ntx = []
        nty = []
        vtx = []
        vty = []
        for s in range(tp.nspins):
            nts = aa1d(nt_sG[s])
            vts = aa1d(vt_sG[s])
            ntsx = aa2d(nt_sG[s], 0)            
            vtsx = aa2d(vt_sG[s], 0)
            ntsy = aa2d(nt_sG[s], 1)            
            vtsy = aa2d(vt_sG[s], 1)                
            nt.append(nts)
            vt.append(vts)
            ntx.append(ntsx)
            vtx.append(vtsx)
            nty.append(ntsy)
            vty.append(vtsy)
        nt = np.array(nt)
        vt = np.array(vt)
        ntx = np.array(ntx)
        vtx = np.array(vtx)
        nty = np.array(nty)
        vty = np.array(vty)
        return nt, vt, ntx, vtx, nty, vty
    
    def calculate_current(self, tp, tc_array, lead_pair_index=0):             
        current = np.array([0, 0], complex)
        contour = tp.contour
        kt = 0.02
        fd = fermidistribution        
        lead_ef1 = contour.leadfermi[self.lead_pairs[lead_pair_index][0]]
        lead_ef2 = contour.leadfermi[self.lead_pairs[lead_pair_index][1]]
        if lead_ef2 > lead_ef1:
            lead_ef1, lead_ef2 = lead_ef2, lead_ef1
        interval = np.real(self.my_energies[1] - self.my_energies[0])
        tc_all = np.sum(tc_array, axis=1) / tp.npk 
        fermi_factor = fd(self.my_energies - lead_ef1, kt) - fd(
                                         self.my_energies - lead_ef2, kt)
        for s in range(tp.my_nspins):
            current[s] = np.sum(tc_all[s, lead_pair_index] * fermi_factor *
                                                          self.my_weights)
        current = np.array(current)
        #tp.contour.comm.sum(current)
        return current 
         
class Transport_Plotter:
    def __init__(self, mode=0):
        self.my_options = False
        self.ion_step = 0
        self.bias_step = 0
        self.ele_step = 0
        self.plot_setup()
        self.initialize()
        self.mode = mode
    
    def initialize(self):
        self.xlabels = {}
        self.ylabels = {}
        names = ['tc', 'dos', 'nt', 'vt']
        xls = ['Energy(eV)', 'Energy(eV)', 'Transport Axis', 'Transport Axis']
        yls = ['Transmission Coefficient', 'Density of States(Electron/eV)',
               'Electron Density', 'Effective Potential']
        for name, xl, yl in zip(names, xls, yls):
            self.xlabels[name] = xl
            self.ylabels[name] = yl

    def read_overhead(self):
        fd = file('analysis_overhead', 'r')
        atoms, basis_information, contour_information = cPickle.load(fd)
        fd.close()
        self.atoms = atoms
        self.basis = basis_information
        self.contour = contour_information
 
    def plot_setup(self):
        from matplotlib import rcParams
        rcParams['xtick.labelsize'] = 18
        rcParams['ytick.labelsize'] = 18
        rcParams['legend.fontsize'] = 18
        rcParams['axes.titlesize'] = 18
        rcParams['axes.labelsize'] = 18
        rcParams['font.size'] = 18
        rcParams['figure.subplot.left'] = 0.15

    def get_data(self, bias_step, ion_step):
        if self.mode == 0:
            fd = file('analysis_data/ionic_step_' + str(ion_step) +
                  '/bias_step_' + str(bias_step), 'r')
        else:
            fd = file('analysis_data/ionic_step_' + str(ion_step) +
                  '/abias_step_' + str(bias_step), 'r')            
        data = cPickle.load(fd)
        fd.close()        
        return data
    
    def get_info(self, name, bias_step, ion_step=0):
        data = self.get_data(bias_step, ion_step)
        if 'newform' in data and data['newform'] == True:
            info = self.collect_info(data, name)
        else:
            if name in ['tc', 'dos']:
                info = data['ER_0_' + name]        
                data_name = 'ER_1_' + name
                n = 1
                while data_name in data:
                    info = np.append(info, data[data_name], axis=-1)
                    n += 1
                    data_name = 'ER_' + str(n) + '_' + name        
            else:
                info = data[name]
        return info

    def collect_info(self, data, name):
        kpt_size = data['kpt_size']
        domain_parsize = data['domain_parsize']
        domain_size = np.product(domain_parsize)
        if name == 'vt' or name == 'nt':
            vt = []
            for i in range(domain_parsize[2]):
                vt.append([])
            for i in range(domain_size):
                d1, d2, d3 =  (i // (domain_parsize[1] * domain_parsize[2]),
                               i % (domain_parsize[1] * 
                               domain_parsize[2]) // domain_parsize[2], 
                               i % domain_parsize[2])        
                lvt = data['K_0D_' + str(i) + '_' + name]
                #lvt = np.sum(lvt, axis=1) / lvt.shape[1]
                #lvt = np.sum(lvt, axis=1) / lvt.shape[1]
                vt[d3].append(lvt) 
            for i in range(domain_parsize[2]):
                vt[i] = np.sum(vt[i],axis=0) / len(vt[i])
            output = np.concatenate(vt, axis=-1)
        elif name =='vty' or name =='nty':
            vty = []
            for i in range(domain_parsize[0]):
                vty.append([])
                for j in range(domain_parsize[2]):
                    vty[i].append([])
            for i in range(domain_size):
                d1, d2, d3 =  (i // (domain_parsize[1] * domain_parsize[2]),
                               i % (domain_parsize[1] * 
                               domain_parsize[2]) // domain_parsize[2], 
                               i % domain_parsize[2])        
                lvty = data['K_0D_' + str(i) + '_' + name]
                vty[d1][d3].append(lvty) 
            for i in range(domain_parsize[0]):
                for j in range(domain_parsize[2]):
                    vty[i][j] = np.sum(vty[i][j],axis=0) / len(vty[i][j])
            output = np.concatenate(vty, axis=-1)           
            output = np.concatenate(output, axis=-1)        
        elif name =='vtx' or name =='ntx':
            vtx = []
            for i in range(domain_parsize[1]):
                vtx.append([])
                for j in range(domain_parsize[2]):
                    vtx[i].append([])
            for i in range(domain_size):
                d1, d2, d3 =  (i // (domain_parsize[1] * domain_parsize[2]),
                               i % (domain_parsize[1] * 
                               domain_parsize[2]) // domain_parsize[2], 
                               i % domain_parsize[2])        
                lvtx = data['K_0D_' + str(i) + '_' + name]
                vtx[d2][d3].append(lvtx) 
            for i in range(domain_parsize[1]):
                for j in range(domain_parsize[2]):
                    vtx[i][j] = np.sum(vtx[i][j],axis=0) / len(vtx[i][j])
            output = np.concatenate(vtx, axis=-1)           
            output = np.concatenate(output, axis=-1)        
        elif name == 'vHt' or name == 'rho':
            vHt = []
            for i in range(domain_parsize[2]):
                vHt.append([])
            for i in range(domain_size):
                d1, d2, d3 =  (i // (domain_parsize[1] * domain_parsize[2]),
                               i % (domain_parsize[1] * 
                               domain_parsize[2]) // domain_parsize[2], 
                               i % domain_parsize[2])        
                lvHt = data['K_0D_' + str(i) + '_' + name]
                lvHt = np.sum(lvHt, axis=0) / lvHt.shape[0]
                lvHt = np.sum(lvHt, axis=0) / lvHt.shape[0]
                vHt[d3].append(lvHt) 
            for i in range(domain_parsize[2]):
                vHt[i] = np.sum(vHt[i],axis=0) / vHt[i].shape[0]
            output = np.concatenate(vHt, axis=-1)
        elif name == 'tc':
            ns, npk, npl, ne = data['transmission_dimension']
            tc = []
            for i in range(kpt_size):
                ltc = []
                for j in range(domain_size):
                    ltc.append(data['K_' + str(i) + 'D_' + str(j) + '_tc'])
                ltc = np.concatenate(ltc, axis=-1)
                tc.append(ltc)
            output = np.array(tc)
            output.shape = (ns, npk, npl, ne)
        elif name == 'dos':
            ns, npk, nb, ne = data['dos_dimension']
            dos = []
            for i in range(kpt_size):
                ldos = []
                for j in range(domain_size):
                    ldos.append(data['K_' + str(i) + 'D_' + str(j) + '_dos'])
                ldos = np.concatenate(ldos, axis=-1)
                dos.append(ldos)
            output = np.array(dos)
            output.shape = (ns, npk, nb, ne)
        elif name == 'charge':
            data = cPickle.load(file('temperary_data/KC_0_DC_0AD_bias_0','r'))
            data1 = cPickle.load(file('temperary_data/KC_1_DC_0AD_bias_0','r'))
            output = [data['charge'], data1['charge']]
            output=np.sum(output,axis=1)
        elif name == 'current':
            # temperary lines
            # ! THa: Hack, load energy scale from saved file
            try:
                ee  = np.load('plot_energy_range.npy')
                enp = len(ee)
            except IOError:
                print('WARNING: could not find file plot_energy_range.npy,\
                        fall back to default: np.linspoace(-5,5,201)')
                ee = np.linspace(-5,5,201)
            enp = len(ee)
            interval = np.abs(ee[0] - ee[-1]) / enp
            # ! THa: Hack end
            kt = 0.02
            lead_ef1 = data['bias'][0]
            lead_ef2 = data['bias'][1]
            ns, npk, npl, ne = data['transmission_dimension']
            tc = []
            for i in range(kpt_size):
                ltc = []
                for j in range(domain_size):
                    ltc.append(data['K_' + str(i) + 'D_' + str(j) + '_tc'])
                ltc = np.concatenate(ltc, axis=-1)
                tc.append(ltc)
            tc = np.array(tc)
            tc.shape = (ns, npk,  ne)
            tc = np.sum(tc, axis=1) / npk
            fd = fermidistribution
            fermi_factor = fd(ee - lead_ef1, kt) - fd(ee - lead_ef2, kt)
            output = np.sum(tc * fermi_factor * interval, axis=-1)
        else:
            output = data[name]
        return output       

    def get_positions(self, ion_step=0):
        fd = file('analysis_data/ionic_step_' + str(ion_step) + '/positions', 'r')
        positions = cPickle.load(fd)
        fd.close()
        return positions

    def process(self, info, s=0, k=None, lp=None):
        if s is None:
            info = np.sum(info, axis=0) / info.shape[0]
        else:
            info = info[s]
        if k is None:
            info = np.sum(info, axis=0) / info.shape[0]
        else:
            info = info[k]
        if lp is not None:
            info = info[lp]
        return info
        
    def tc(self, bias_step, ion_step=0, lp=0, s=0, k=None):
        info = self.get_info('tc', bias_step, ion_step)
        return self.process(info, s, k, lp)
    
    def dos_array(self, bias_step, ion_step=0, s=0, k=None):
        info = self.get_info('dos', bias_step, ion_step)
        return self.process(info, s, k)
    
    def dos(self, bias_step, ion_step=0, s=0, k=None):
        info = self.get_info('dos', bias_step, ion_step)
        info = self.process(info, s, k)
        return np.sum(info, axis=-2)
 
    def partial_dos(self, bias_step, ion_step=0, s=0, k=None,
                    atom_indices=None, orbital_type=None, direction=None):
        self.read_overhead()
        dos_array = self.dos_array(bias_step, ion_step, s, k)
        orbital_indices = self.basis['orbital_indices']
        orbital_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3}
        direction_map = {'x': 2, 'y': 0, 'z': 1, 'xy': 0, 'yz': 1,
                         'z2r2': 2, 'xz': 3, 'x2y2': 4}        
        if direction is None:
            direction_index = np.zeros([orbital_indices.shape[0]]) + 1
        else:
            directions = []
            num = -1
            oio = 0
            for oi in orbital_indices:
                if oi[1] != oio:
                    num = -1
                oio = oi[1]      
                if oi[1] != orbital_map[orbital_type]:
                    directions.append(0)
                else:
                    num += 1
                    if num == direction_map[direction]:
                        directions.append(1)
                    else:
                        directions.append(0)
            direction_index = np.array(directions)
            
        if orbital_type is None:
            orbital_index = np.zeros([orbital_indices.shape[0]]) + 1
        else:
            orbital_index = orbital_indices[:, 1] - orbital_map[
                                                      orbital_type] ==  0
        atom_index = np.zeros([orbital_indices.shape[0]]) + 1
        if atom_indices is not None:
            atom_index -= 1
            for i in atom_indices:
                atom_index += orbital_indices[:, 0] - i == 0
        pdos = []
        for i in range(dos_array.shape[1]):
            pdos.append(np.sum(dos_array[:,i] * orbital_index * atom_index * direction_index))
        pdos = np.array(pdos)
        return pdos

    def charge(self, bias_step, ion_step=0, s=0, atom_indices=None,
               orbital_type=None, direction=None):
        self.read_overhead()
        charge_array = self.get_info('charge', bias_step, ion_step)
        orbital_indices = self.basis['orbital_indices']
        orbital_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3}
        direction_map = {'x': 2, 'y': 0, 'z': 1, 'xy': 0, 'yz': 1,
                         'z2r2': 2, 'xz': 3, 'x2y2': 4}        
        if direction is None:
            direction_index = np.zeros([orbital_indices.shape[0]]) + 1
        else:
            directions = []
            num = -1
            oio = 0
            for oi in orbital_indices:
                if oi[1] != oio:
                    num = -1
                oio = oi[1]      
                if oi[1] != orbital_map[orbital_type]:
                    directions.append(0)
                else:
                    num += 1
                    if num == direction_map[direction]:
                        directions.append(1)
                    else:
                        directions.append(0)
            direction_index = np.array(directions)

        if orbital_type is None:
            orbital_index = np.zeros([orbital_indices.shape[0]]) + 1
        else:
            orbital_index = orbital_indices[:, 1] - orbital_map[
                                                      orbital_type] ==  0
        n_atom_basis = len(orbital_index)
        charge_array = charge_array[s, :n_atom_basis]
        atom_index = np.zeros([orbital_indices.shape[0]]) + 1
        if atom_indices is not None:
            atom_index -= 1
            for i in atom_indices:
                atom_index += orbital_indices[:, 0] - i == 0
        charge = np.sum(charge_array * atom_index * orbital_index * direction_index)
        return charge
       
    def iv(self, nsteps=16, spinpol=False):
        current = []
        bias = []
        for i in range(nsteps):
            bias_list = self.get_info('bias', i)
            bias.append(bias_list[0]-bias_list[1])
            current.append(self.get_info('current', i))
        unit = 6.624 * 1e3 
        current = np.array(current) * unit / (Hartree * 2 * np.pi)
        if not spinpol:
            current *= 2
        bias = np.array(bias)
        return bias, current

    def diamond(self, bias_step, ionic_step=0, dense_level=1, spinpol=False):
        tc = self.tc(bias_step, ionic_step)
        if dense_level > 1:
            from scipy import interpolate
            ee = np.linspace(-5, 5, 201)
            new_ee = np.linspace(-5, 5, dense_level * 201)
            tck = interpolate.splrep(ee, tc, s=0)
            new_tc = interpolate.splev(new_ee, tck, der=0)
            tc = new_tc
        de = 0.05 / dense_level
        nb = 30 * dense_level
        ng = 140 * dense_level
        currents = np.zeros([nb, ng])
        for i in range(ng):
            for j in range(nb):
                begin = i + nb - j
                end = i + nb + j
                currents[j, i] = np.sum(tc[begin:end])
        unit = 6.624 * 1e3 
        currents = np.array(currents) * de * unit / (Hartree * 2 * np.pi)
        if not spinpol:
            currents *= 2
        return currents
           
    def tvs(self, nsteps=16, spinpol=True):
        bias, current = self.iv(nsteps,  spinpol)
        current *= 1e-6 # switch unit to Ampier
        ydata = np.log(abs(current) / (bias ** 2))
        xdata = 1 / bias
        return xdata, ydata



