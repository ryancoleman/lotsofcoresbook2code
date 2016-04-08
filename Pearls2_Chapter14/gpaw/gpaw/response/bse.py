from __future__ import print_function

import pickle
from math import pi
from time import time, ctime

import numpy as np
from ase.io import write
from ase.units import Hartree
from ase.utils import opencew

from gpaw.response.df0 import DF
from gpaw.io.tar import Writer, Reader
from gpaw.response.base import BASECHI
from gpaw.utilities.memory import maxrss
from gpaw.blacs import BlacsGrid, Redistributor
from gpaw.mpi import world, size, rank, serial_comm
from gpaw.response.kernel import calculate_Kc, calculate_Kc_q
from gpaw.response.parallel import parallel_partition, gatherv


class BSE(BASECHI):
    """This class defines Bethe-Salpeter equations."""

    def __init__(self,
                 calc=None,
                 nbands=None,
                 nc=None,
                 nv=None,
                 w=None,
                 q=None,
                 eshift=None,
                 ecut=10.,
                 eta=0.2,
                 gw_skn=None, # GW QP energies in Hartree
                 rpad=None,
                 vcut=None,   # Coulomb cutoff only 2D works
                 ftol=1e-5,
                 txt=None,
                 optical_limit=None,
                 integrate_coulomb=None,
                 print_coulomb=False,
                 coupling=False,  # False : use Tamm-Dancoff Approx
                 mode='BSE',      # BSE, TDHF or RPA
                 kernel_file=None,#'W_qGG',
                 qsymm=True): 

        BASECHI.__init__(self, calc=calc, nbands=nbands, w=w, q=q,
                         eshift=eshift, ecut=ecut, eta=eta, rpad=rpad,
                         ftol=ftol, txt=txt, optical_limit=optical_limit)

        assert mode in ['RPA', 'TDHF', 'BSE']

        self.epsilon_w = None
        self.coupling = coupling
        self.vcut = vcut
        self.nc = nc # conduction band index
        self.nv = nv # valence band index
        self.gw_skn = gw_skn
        self.mode = mode
        self.integrate_coulomb = integrate_coulomb
        self.print_coulomb = print_coulomb
        self.kernel_file = kernel_file
        self.qsymm = qsymm

    def initialize(self):

        self.printtxt('----------------------------------------')
        self.printtxt('Bethe-Salpeter Equation calculation')
        self.printtxt('----------------------------------------')
        self.printtxt('Started at:  %s' % ctime())
        self.printtxt('')
        BASECHI.initialize(self)
        assert self.nspins == 1
        
        calc = self.calc
        self.kd = kd = calc.wfs.kd

        # frequency points init
        self.dw = self.w_w[1] - self.w_w[0]
        assert ((self.w_w[1:] - self.w_w[:-1] - self.dw) < 1e-10).all() # make sure its linear w grid
        assert self.w_w.max() == self.w_w[-1]

        self.dw /= Hartree
        self.w_w  /= Hartree
        self.wmax = self.w_w[-1] 
        self.Nw  = int(self.wmax / self.dw) + 1

        # band init
        if self.nc is None:
            nv = self.nvalence / 2 - 1
            self.nv = np.array([nv, nv+1]) # conduction band start / end
            self.nc = np.array([nv+1, nv+2]) # valence band start / end

        self.printtxt('')
        self.printtxt('Number of electrons      : %d' % (self.nvalence))
        self.printtxt('Valence band included    : (band %d to band %d)' %(self.nv[0], self.nv[1]-1))
        self.printtxt('Conduction band included : (band %d to band %d)' %(self.nc[0], self.nc[1]-1))
        if self.eshift is not None:
            self.printtxt('Scissors operator        : %2.3f eV' % self.eshift)
        self.printtxt('')
            
        # find the pair index and initialized pair energy (e_i - e_j) and occupation(f_i-f_j)
        self.e_S = {}
        focc_s = {}
        self.Sindex_S3 = {}
        iS = 0
        kq_k = self.kq_k
        for k1 in range(self.kd.nbzkpts):
            ibzkpt1 = kd.bz2ibz_k[k1]
            ibzkpt2 = kd.bz2ibz_k[kq_k[k1]]
            for n1 in range(self.nv[0], self.nv[1]): 
                for m1 in range(self.nc[0], self.nc[1]): 
                    focc = self.f_skn[0][ibzkpt1,n1] - self.f_skn[0][ibzkpt2,m1]
                    if self.coupling: # Dont use Tamm-Dancoff Approx.
                        check_ftol = np.abs(focc) > self.ftol
                    else:
                        check_ftol = focc > self.ftol
                    if check_ftol:
                        if self.gw_skn is None:
                            self.e_S[iS] = self.e_skn[0][ibzkpt2,m1] - self.e_skn[0][ibzkpt1,n1]
                        else:
                            self.e_S[iS] = self.gw_skn[0][ibzkpt2,m1] - self.gw_skn[0][ibzkpt1,n1]
                            
                        focc_s[iS] = focc
                        self.Sindex_S3[iS] = (k1, n1, m1)
                        iS += 1
        self.nS = iS
        self.focc_S = np.zeros(self.nS)
        for iS in range(self.nS):
            self.focc_S[iS] = focc_s[iS]

        # q points init
        self.bzq_qc = kd.get_bz_q_points()
        if not self.qsymm:
            self.ibzq_qc = self.bzq_qc
        else:
            (self.ibzq_qc, self.ibzq_q, self.iop_q,
             self.timerev_q, self.diff_qc) = kd.get_ibz_q_points(self.bzq_qc,
                                                                 calc.wfs.kd.symmetry.op_scc)
            if np.abs(self.bzq_qc - kd.bzk_kc).sum() < 1e-8:
                assert np.abs(self.ibzq_qc - kd.ibzk_kc).sum() < 1e-8
        self.nibzq = len(self.ibzq_qc)

        # Parallel initialization 
        # kcomm and wScomm is only to be used when wavefunctions are distributed in parallel.
        self.comm = self.Scomm = world
        self.kcomm = world
        self.wScomm = serial_comm
        self.nS, self.nS_local, self.nS_start, self.nS_end = parallel_partition(
            self.nS, world.rank, world.size, reshape=False)

        self.print_bse()

        if calc.input_parameters['mode'] == 'lcao':
            calc.initialize_positions()

        # Coulomb interaction at q=0 for Hartree coupling
        ### 2D z direction only !!!!!!!!!!!!!!!!!!!!!!
        if self.integrate_coulomb is None:
            self.integrate_coulomb = []
            if self.vcut is None:
                pass
            elif self.vcut == '2D':
                for iG in range(len(self.Gvec_Gc)):
                    if self.Gvec_Gc[iG, 0] == 0 and self.Gvec_Gc[iG, 1] == 0:
                        self.integrate_coulomb.append(iG)
            else:
                raise NotImplementedError
        elif type(self.integrate_coulomb) is int:
            self.integrate_coulomb = range(self.integrate_coulomb)
        elif self.integrate_coulomb == 'all':
            self.integrate_coulomb = range(len(self.Gvec_Gc))
        elif type(self.integrate_coulomb) is list:
            pass
        else:
            raise 'Invalid option for integrate_coulomb'
        
        self.printtxt('')
        self.printtxt('Calculating bare Coulomb kernel')
        if not len(self.integrate_coulomb) == 0:
            self.printtxt('Integrating Coulomb kernel at %s reciprocal lattice vector(s)' % len(self.integrate_coulomb))
        
        # Coulomb interaction at problematic G's for exchange coupling
        if len(self.integrate_coulomb) != 0:
            self.vint_Gq = []
            for iG in self.integrate_coulomb:
                v_q, v0_q = calculate_Kc_q(self.acell_cv,
                                           self.bcell_cv,
                                           self.pbc,
                                           self.kd.N_c,
                                           vcut=self.vcut,
                                           Gvec_c=self.Gvec_Gc[iG],
                                           q_qc=self.ibzq_qc.copy())
                self.vint_Gq.append(v_q)                    
                if self.print_coulomb:
                    self.printtxt('')
                    self.printtxt('Average kernel relative to bare kernel - \int v(q)dq / v(q0): ')
                    self.printtxt('  G: % s' % self.Gvec_Gc[iG])
                    for iq in range(len(v_q)):
                        q_s = '    q = [%1.2f,  %1.2f,  %1.2f]: ' % (self.ibzq_qc[iq,0],
                                                                     self.ibzq_qc[iq,1],
                                                                     self.ibzq_qc[iq,2])
                        v_rel = v_q[iq] / v0_q[iq]
                        self.printtxt(q_s + '%1.3f' % v_rel)
        self.printtxt('')

        self.V_qGG = self.full_bare_interaction()

    def calculate(self):
        calc = self.calc
        focc_S = self.focc_S
        e_S = self.e_S
        op_scc = calc.wfs.kd.symmetry.op_scc

        # Get phi_qaGp
        if self.mode == 'RPA':
            self.phi_aGp = self.get_phi_aGp()
        else:
            fd = opencew('phi_qaGp')
            if fd is None:
                self.reader = Reader('phi_qaGp')
                tmp = self.load_phi_aGp(self.reader, 0)[0]
                assert len(tmp) == self.npw
                self.printtxt('Finished reading phi_aGp')
            else:
                self.printtxt('Calculating phi_qaGp')
                self.get_phi_qaGp()
                world.barrier()
                self.reader = Reader('phi_qaGp')                
            self.printtxt('Memory used %f M' % (maxrss() / 1024.**2))
            self.printtxt('')

        if self.optical_limit:
            iq = np.where(np.sum(abs(self.ibzq_qc), axis=1) < 1e-5)[0][0]
        else:
            iq = np.where(np.sum(abs(self.ibzq_qc - self.q_c), axis=1) < 1e-5)[0][0]
        kc_G = np.array([self.V_qGG[iq, iG, iG] for iG in range(self.npw)])
        if self.optical_limit:
            kc_G[0] = 0.

        # Get screened Coulomb kernel
        if self.mode == 'BSE':
            try:
                # Read
                data = pickle.load(open(self.kernel_file+'.pckl'))
                W_qGG = data['W_qGG']
                assert np.shape(W_qGG) == np.shape(self.V_qGG)
                self.printtxt('Finished reading screening interaction kernel')
            except:
                # Calculate from scratch
                self.printtxt('Calculating screening interaction kernel.')
                W_qGG = self.full_static_screened_interaction()
            self.printtxt('')
        else:
            W_qGG = self.V_qGG
 
        t0 = time()
        self.printtxt('Calculating %s matrix elements' % self.mode)

        # Calculate full kernel
        K_SS = np.zeros((self.nS_local, self.nS), dtype=complex)
        self.rhoG0_S = np.zeros(self.nS, dtype=complex)

        #noGmap = 0
        for iS in range(self.nS_start, self.nS_end):
            k1, n1, m1 = self.Sindex_S3[iS]
            rho1_G = self.density_matrix(n1,m1,k1)
            self.rhoG0_S[iS] = rho1_G[0]
            for jS in range(self.nS):
                k2, n2, m2 = self.Sindex_S3[jS]
                rho2_G = self.density_matrix(n2,m2,k2)
                K_SS[iS-self.nS_start, jS] = np.sum(rho1_G.conj() * rho2_G * kc_G)

                if not self.mode == 'RPA':
                    rho3_G = self.density_matrix(n1,n2,k1,k2)
                    rho4_G = self.density_matrix(m1,m2,self.kq_k[k1],
                                                 self.kq_k[k2])

                    q_c = self.kd.bzk_kc[k2] - self.kd.bzk_kc[k1]
                    q_c[np.where(q_c > 0.501)] -= 1.
                    q_c[np.where(q_c < -0.499)] += 1.
                    iq = self.kd.where_is_q(q_c, self.bzq_qc)
                    
                    if not self.qsymm:    
                        W_GG = W_qGG[iq]
                    else:
                        ibzq = self.ibzq_q[iq]
                        W_GG_tmp = W_qGG[ibzq]

                        iop = self.iop_q[iq]
                        timerev = self.timerev_q[iq]
                        diff_c = self.diff_qc[iq]
                        invop = np.linalg.inv(op_scc[iop])
                        Gindex = np.zeros(self.npw, dtype=int)
                        for iG in range(self.npw):
                            G_c = self.Gvec_Gc[iG]
                            if timerev:
                                RotG_c = -np.int8(np.dot(invop, G_c+diff_c).round())
                            else:
                                RotG_c = np.int8(np.dot(invop, G_c+diff_c).round())
                            tmp_G = np.abs(self.Gvec_Gc - RotG_c).sum(axis=1)
                            try:
                                Gindex[iG] = np.where(tmp_G < 1e-5)[0][0]
                            except:
                                #noGmap += 1
                                Gindex[iG] = -1
    
                        W_GG = np.zeros_like(W_GG_tmp)
                        for iG in range(self.npw):
                            for jG in range(self.npw):
                                if Gindex[iG] == -1 or Gindex[jG] == -1:
                                    W_GG[iG, jG] = 0
                                else:
                                    W_GG[iG, jG] = W_GG_tmp[Gindex[iG], Gindex[jG]]
                                    
                    if self.mode == 'BSE':
                        tmp_GG = np.outer(rho3_G.conj(), rho4_G) * W_GG
                        K_SS[iS-self.nS_start, jS] -= 0.5 * np.sum(tmp_GG)
                    else:
                        tmp_G = rho3_G.conj() * rho4_G * np.diag(W_GG)
                        K_SS[iS-self.nS_start, jS] -= 0.5 * np.sum(tmp_G)
            self.timing(iS, t0, self.nS_local, 'pair orbital') 
 
        K_SS /= self.vol

        world.sum(self.rhoG0_S)
        #self.printtxt('Number of G indices outside the Gvec_Gc: %d' % noGmap)

        # Get and solve Hamiltonian
        H_sS = np.zeros_like(K_SS)
        for iS in range(self.nS_start, self.nS_end):
            H_sS[iS-self.nS_start,iS] = e_S[iS]
            for jS in range(self.nS):
                H_sS[iS-self.nS_start,jS] += focc_S[iS] * K_SS[iS-self.nS_start,jS]
  
        # Force matrix to be Hermitian
        if not self.coupling:
            if world.size > 1:
                H_Ss = self.redistribute_H(H_sS)
            else:
                H_Ss = H_sS
            H_sS = (np.real(H_sS) + np.real(H_Ss.T)) / 2. + 1j * (np.imag(H_sS) - np.imag(H_Ss.T)) /2.

        # Save H_sS matrix
        self.par_save('H_SS','H_SS', H_sS)

        return H_sS

    def diagonalize(self, H_sS):
        if self.coupling: # Non-Hermitian matrix can only use linalg.eig
            self.printtxt('Use numpy.linalg.eig')
            H_SS = np.zeros((self.nS, self.nS), dtype=complex)
            if self.nS % world.size == 0:
                world.all_gather(H_sS, H_SS)
            else:
                H_SS = gatherv(H_sS)

            self.w_S, self.v_SS = np.linalg.eig(H_SS)
            self.par_save('v_SS', 'v_SS', self.v_SS[self.nS_start:self.nS_end, :].copy())
        else:
            if world.size == 1:
                self.printtxt('Use lapack.')
                from gpaw.utilities.lapack import diagonalize
                self.w_S = np.zeros(self.nS)
                H_SS = H_sS
                diagonalize(H_SS, self.w_S)
                self.v_SS = H_SS.conj() # eigenvectors in the rows, transposed later
            else:
                self.printtxt('Use scalapack')
                self.w_S, self.v_sS = self.scalapack_diagonalize(H_sS)
                self.v_SS = self.v_sS # just use the same name
            self.par_save('v_SS', 'v_SS', self.v_SS)
        
        return 


    def par_save(self,filename, name, A_sS):
        from gpaw.io import open 

        nS = self.nS
        
        if rank == 0:
            w = open(filename, 'w', world)
            w.dimension('nS', nS)
            
            if name == 'v_SS':
                w.add('w_S', ('nS',), dtype=self.w_S.dtype)
                w.fill(self.w_S)
            w.add('rhoG0_S', ('nS',), dtype=complex)
            w.fill(self.rhoG0_S)
            w.add(name, ('nS', 'nS'), dtype=complex)

            tmp = np.zeros_like(A_sS)

        # Assumes that H_SS is written in order from rank 0 - rank N
        for irank in range(size):
            if irank == 0:
                if rank == 0:
                    w.fill(A_sS)
            else:
                if rank == irank:
                    world.send(A_sS, 0, irank+100)
                if rank == 0:
                    world.receive(tmp, irank, irank+100)
                    w.fill(tmp)
        if rank == 0:
            w.close()
        world.barrier()

    def par_load(self,filename, name):
        from gpaw.io import open 
        r = open(filename, 'r')
        nS = r.dimension('nS')

        if name == 'v_SS':
            self.w_S = r.get('w_S')

        A_SS = np.zeros((self.nS_local, nS), dtype=complex)
        for iS in range(self.nS_start, self.nS_end):   
            A_SS[iS-self.nS_start,:] = r.get(name, iS)

        self.rhoG0_S = r.get('rhoG0_S')

        r.close()
        
        return A_SS
    

    def full_static_screened_interaction(self):
        """Calcuate W_GG(q)"""
        W_qGG = np.zeros((self.nibzq, self.npw, self.npw), dtype=complex)

        t0 = time()
        for iq in range(self.nibzq):
            q = self.ibzq_qc[iq]
            optical_limit = False
            if np.abs(q).sum() < 1e-8:
                q = self.q_c.copy()
                optical_limit = True
            df = DF(calc=self.calc,
                    q=q,
                    w=(0.,),
                    optical_limit=optical_limit,
                    nbands=self.nbands,
                    hilbert_trans=False,
                    eta=0.0001,
                    ecut=self.ecut*Hartree,
                    xc='RPA',
                    txt='df.out')
            df.initialize()
            df.calculate()

            if optical_limit:
                K_GG = self.V_qGG[iq].copy()
                K0 = calculate_Kc(q,
                                  self.Gvec_Gc,
                                  self.acell_cv,
                                  self.bcell_cv,
                                  self.pbc,
                                  vcut=self.vcut)[0,0]

                for iG in range(1,self.npw):
                    K_GG[0, iG] = self.V_qGG[iq, iG, iG]**0.5 * K0**0.5
                    K_GG[iG, 0] = self.V_qGG[iq, iG, iG]**0.5 * K0**0.5
                K_GG[0,0] = K0
                df_GG = np.eye(self.npw, self.npw) - K_GG*df.chi0_wGG[0]
            else:
                df_GG = np.eye(self.npw, self.npw) - self.V_qGG[iq]*df.chi0_wGG[0]
            dfinv_GG = np.linalg.inv(df_GG)
            
            if optical_limit:
                eps = 1/dfinv_GG[0,0]
                self.printtxt('    RPA macroscopic dielectric constant is: %3.3f' %  eps.real)
            W_qGG[iq] = dfinv_GG * self.V_qGG[iq]
            self.timing(iq, t0, self.nibzq, 'iq')
            
        if rank == 0:
            if self.kernel_file is not None:
                data = {'W_qGG': W_qGG}
                name = self.kernel_file+'.pckl'
                pickle.dump(data, open(name, 'w'), -1)
        
        return W_qGG


    def full_bare_interaction(self):
        """Calculate V_GG(q)"""

        V_qGG = np.zeros((self.nibzq, self.npw, self.npw), dtype=complex)
        
        for iq in range(self.nibzq):
            q = self.ibzq_qc[iq]

            Vq_G = np.diag(calculate_Kc(q,
                                        self.Gvec_Gc,
                                        self.acell_cv,
                                        self.bcell_cv,
                                        self.pbc,
                                        integrate_gamma=True,
                                        N_k=self.kd.N_c,
                                        vcut=self.vcut))**0.5
            
            for i, iG in enumerate(self.integrate_coulomb):
                Vq_G[iG] = self.vint_Gq[i][iq]**0.5

            V_qGG[iq] = np.outer(Vq_G, Vq_G)

        return V_qGG

                                          
    def print_bse(self):

        printtxt = self.printtxt

        if not self.mode == 'RPA':
            printtxt('Number of q points           : %d' %(self.nibzq))
        printtxt('Number of frequency points   : %d' %(self.Nw) )
        printtxt('Number of pair orbitals      : %d' %(self.nS) )
        printtxt('')
        printtxt('Parallelization scheme:')
        printtxt('   Total cpus         : %d' %(world.size))
        printtxt('   pair orb parsize   : %d' %(self.Scomm.size))        
        
        return


    def get_phi_qaGp(self):

        N1_max = 0
        N2_max = 0
        natoms = len(self.calc.wfs.setups)
        for id in range(natoms):
            N1 = self.npw
            N2 = self.calc.wfs.setups[id].ni**2
            if N1 > N1_max:
                N1_max = N1
            if N2 > N2_max:
                N2_max = N2
        
        nbzq = self.kd.nbzkpts
        nbzq, nq_local, q_start, q_end = parallel_partition(
                                   nbzq, world.rank, world.size, reshape=False)
        phimax_qaGp = np.zeros((nq_local, natoms, N1_max, N2_max), dtype=complex)
        #phimax_qaGp = np.zeros((nbzq, natoms, N1_max, N2_max), dtype=complex)

        t0 = time()
        for iq in range(nq_local):
            q_c = self.bzq_qc[iq + q_start]
            tmp_aGp = self.get_phi_aGp(q_c, parallel=False)
            for id in range(natoms):
                N1, N2 = tmp_aGp[id].shape
                phimax_qaGp[iq, id, :N1, :N2] = tmp_aGp[id]
            self.timing(iq*world.size, t0, nq_local, 'iq')
        world.barrier()

        # Write to disk
        filename = 'phi_qaGp'
        if world.rank == 0:
            w = Writer(filename)
            w.dimension('nbzq', nbzq)
            w.dimension('natoms', natoms)
            w.dimension('nG', N1_max)
            w.dimension('nii', N2_max)
            w.add('phi_qaGp', ('nbzq', 'natoms', 'nG', 'nii',), dtype=complex)

        for q in range(nbzq):
            residual = nbzq % size
            N_local = nbzq // size
            if q < residual * (N_local + 1):
                qrank = q // (N_local + 1)
            else:
                qrank = (q - residual * (N_local + 1)) // N_local + residual
                
            if qrank == 0:
                if world.rank == 0:
                    phi_aGp = phimax_qaGp[q - q_start]
            else:
                if world.rank == qrank:
                    phi_aGp = phimax_qaGp[q - q_start]
                    world.send(phi_aGp, 0, q)
                elif world.rank == 0:
                    world.receive(phi_aGp, qrank, q)
            if world.rank == 0:
                w.fill(phi_aGp)
        if world.rank == 0:
            w.close()
        world.barrier()

    def load_phi_aGp(self, reader, iq):

        phimax_aGp = np.array(reader.get('phi_qaGp', iq), complex)

        phi_aGp = {}
        natoms = len(phimax_aGp)
        for a in range(natoms):
            N1 = self.npw
            N2 = self.calc.wfs.setups[a].ni**2
            phi_aGp[a] = phimax_aGp[a, :N1, :N2]

        return phi_aGp

    def get_dielectric_function(self, filename='df.dat', readfile=None):

        if self.epsilon_w is None:
            self.initialize()

            if readfile is None:
                H_sS = self.calculate()
                self.printtxt('Diagonalizing %s matrix.' % self.mode)
                self.diagonalize(H_sS)
                self.printtxt('Calculating dielectric function.')
            elif readfile == 'H_SS':
                H_sS = self.par_load('H_SS', 'H_SS')
                self.printtxt('Finished reading H_SS.gpw')
                self.diagonalize(H_sS)
                self.printtxt('Finished diagonalizing BSE matrix')
            elif readfile == 'v_SS':
                self.v_SS = self.par_load('v_SS', 'v_SS')
                self.printtxt('Finished reading v_SS.gpw')
            else:
                1 / 0

            w_S = self.w_S
            if not self.coupling:
                v_SS = self.v_SS.T # v_SS[:,lamda]
            else:
                v_SS = self.v_SS
            rhoG0_S = self.rhoG0_S
            focc_S = self.focc_S

            # get overlap matrix
            if self.coupling:
                tmp = np.dot(v_SS.conj().T, v_SS )
                overlap_SS = np.linalg.inv(tmp)
    
            # get chi
            epsilon_w = np.zeros(self.Nw, dtype=complex)

            A_S = np.dot(rhoG0_S, v_SS)
            B_S = np.dot(rhoG0_S*focc_S, v_SS)
            if self.coupling:
                C_S = np.dot(B_S.conj(), overlap_SS.T) * A_S
            else:
                if world.size == 1:
                    C_S = B_S.conj() * A_S
                else:
                    tmp = B_S.conj() * A_S
                    C_S = gatherv(tmp, self.nS)

            for iw in range(self.Nw):
                tmp_S = 1. / (iw*self.dw - w_S + 1j*self.eta)
                epsilon_w[iw] += np.dot(tmp_S, C_S)
    
            epsilon_w *=  - 4 * pi / np.inner(self.qq_v, self.qq_v) / self.vol
            epsilon_w += 1        

            self.epsilon_w = epsilon_w
    
        if rank == 0:
            f = open(filename, 'w')
            for iw in range(self.Nw):
                energy = iw * self.dw * Hartree
                print(energy, np.real(epsilon_w[iw]), np.imag(epsilon_w[iw]),
                      file=f)
            f.close()
            #g.close()
        # Wait for I/O to finish
        world.barrier()

        """Check f-sum rule."""
        N1 = 0
        for iw in range(self.Nw):
            w = iw * self.dw
            N1 += np.imag(epsilon_w[iw]) * w
        N1 *= self.dw * self.vol / (2 * pi**2)

        self.printtxt('')
        self.printtxt('Sum rule:')
        nv = self.nvalence
        self.printtxt('N1 = %f, %f  %% error' %(N1, (N1 - nv) / nv * 100) )

        return epsilon_w


    def get_e_h_density(self, lamda=None, filename=None):

        if filename is not None:
            self.load(filename)
            self.initialize()
            
        gd = self.gd
        v_SS = self.v_SS
        A_S = v_SS[:, lamda]
        kq_k = self.kq_k
        kd = self.kd

        # Electron density
        nte_R = gd.zeros()
        
        for iS in range(self.nS_start, self.nS_end):
            print('electron density:', iS)
            k1, n1, m1 = self.Sindex_S3[iS]
            ibzkpt1 = kd.bz2ibz_k[k1]
            psitold_g = self.get_wavefunction(ibzkpt1, n1)
            psit1_g = kd.transform_wave_function(psitold_g, k1)

            for jS in range(self.nS):
                k2, n2, m2 = self.Sindex_S3[jS]
                if m1 == m2 and k1 == k2:
                    psitold_g = self.get_wavefunction(ibzkpt1, n2)
                    psit2_g = kd.transform_wave_function(psitold_g, k1)

                    nte_R += A_S[iS] * A_S[jS].conj() * psit1_g.conj() * psit2_g

        # Hole density
        nth_R = gd.zeros()
        
        for iS in range(self.nS_start, self.nS_end):
            print('hole density:', iS)
            k1, n1, m1 = self.Sindex_S3[iS]
            ibzkpt1 = kd.bz2ibz_k[kq_k[k1]]
            psitold_g = self.get_wavefunction(ibzkpt1, m1)
            psit1_g = kd.transform_wave_function(psitold_g, kq_k[k1])

            for jS in range(self.nS):
                k2, n2, m2 = self.Sindex_S3[jS]
                if n1 == n2 and k1 == k2:
                    psitold_g = self.get_wavefunction(ibzkpt1, m2)
                    psit2_g = kd.transform_wave_function(psitold_g, kq_k[k1])

                    nth_R += A_S[iS] * A_S[jS].conj() * psit1_g * psit2_g.conj()
                    
        self.Scomm.sum(nte_R)
        self.Scomm.sum(nth_R)


        if rank == 0:
            write('rho_e.cube',self.calc.atoms, format='cube', data=nte_R)
            write('rho_h.cube',self.calc.atoms, format='cube', data=nth_R)
            
        world.barrier()
        
        return 

    def get_excitation_wavefunction(self, lamda=None,filename=None, re_c=None, rh_c=None):
        """ garbage at the moment. come back later"""
        if filename is not None:
            self.load(filename)
            self.initialize()
            
        gd = self.gd
        v_SS = self.v_SS
        A_S = v_SS[:, lamda]
        kq_k = self.kq_k
        kd = self.kd

        nx, ny, nz = self.gd.N_c
        nR = 9
        nR2 = (nR - 1 ) // 2
        if re_c is not None:
            psith_R = gd.zeros(dtype=complex)
            psith2_R = np.zeros((nR*nx, nR*ny, nz), dtype=complex)
            
        elif rh_c is not None:
            psite_R = gd.zeros(dtype=complex)
            psite2_R = np.zeros((nR*nx, ny, nR*nz), dtype=complex)
        else:
            self.printtxt('No wavefunction output !')
            return
            
        for iS in range(self.nS_start, self.nS_end):

            k, n, m = self.Sindex_S3[iS]
            ibzkpt1 = kd.bz2ibz_k[k]
            ibzkpt2 = kd.bz2ibz_k[kq_k[k]]
            print('hole wavefunction', iS, (k,n,m),A_S[iS])
            
            psitold_g = self.get_wavefunction(ibzkpt1, n)
            psit1_g = kd.transform_wave_function(psitold_g, k)

            psitold_g = self.get_wavefunction(ibzkpt2, m)
            psit2_g = kd.transform_wave_function(psitold_g, kq_k[k])

            if re_c is not None:
                # given electron position, plot hole wavefunction
                tmp = A_S[iS] * psit1_g[re_c].conj() * psit2_g
                psith_R += tmp

                k_c = self.kd.bzk_kc[k] + self.q_c
                for i in range(nR):
                    for j in range(nR):
                        R_c = np.array([i-nR2, j-nR2, 0])
                        psith2_R[i*nx:(i+1)*nx, j*ny:(j+1)*ny, 0:nz] += \
                                                tmp * np.exp(1j*2*pi*np.dot(k_c,R_c))
                
            elif rh_c is not None:
                # given hole position, plot electron wavefunction
                tmp = A_S[iS] * psit1_g.conj() * psit2_g[rh_c] * self.expqr_g
                psite_R += tmp

                k_c = self.kd.bzk_kc[k]
                k_v = np.dot(k_c, self.bcell_cv)
                for i in range(nR):
                    for j in range(nR):
                        R_c = np.array([i-nR2, 0, j-nR2])
                        R_v = np.dot(R_c, self.acell_cv)
                        assert np.abs(np.dot(k_v, R_v) - np.dot(k_c, R_c) * 2*pi).sum() < 1e-5
                        psite2_R[i*nx:(i+1)*nx, 0:ny, j*nz:(j+1)*nz] += \
                                                tmp * np.exp(-1j*np.dot(k_v,R_v))
                
            else:
                pass

        if re_c is not None:
            self.Scomm.sum(psith_R)
            self.Scomm.sum(psith2_R)
            if rank == 0:
                write('psit_h.cube',self.calc.atoms, format='cube', data=psith_R)

                atoms = self.calc.atoms
                shift = atoms.cell[0:2].copy()
                atoms.cell[0:2] *= nR2
                atoms.positions += shift * (nR2 - 1)
                
                write('psit_bigcell_h.cube',atoms, format='cube', data=psith2_R)
        elif rh_c is not None:
            self.Scomm.sum(psite_R)
            self.Scomm.sum(psite2_R)
            if rank == 0:
                write('psit_e.cube',self.calc.atoms, format='cube', data=psite_R)

                atoms = self.calc.atoms
#                shift = atoms.cell[0:2].copy()
                atoms.cell[0:2] *= nR2
#                atoms.positions += shift * (nR2 - 1)
                
                write('psit_bigcell_e.cube',atoms, format='cube', data=psite2_R)
                
        else:
            pass

        world.barrier()
            
        return
    

    def load(self, filename):

        data = pickle.load(open(filename))
        self.w_S  = data['w_S']
        self.v_SS = data['v_SS']

        self.printtxt('Read succesfully !')
        

    def save(self, filename):
        """Dump essential data"""

        data = {'w_S'  : self.w_S,
                'v_SS' : self.v_SS}
        
        if rank == 0:
            pickle.dump(data, open(filename, 'w'), -1)

        world.barrier()

    def redistribute_H(self, H_sS):

        g1 = BlacsGrid(world, size, 1)
        g2 = BlacsGrid(world, 1, size)
        N = self.nS
        nndesc1 = g1.new_descriptor(N, N, self.nS_local,  N) 
        nndesc2 = g2.new_descriptor(N, N, N, self.nS_local)
        
        H_Ss = nndesc2.empty(dtype=H_sS.dtype)
        redistributor = Redistributor(world, nndesc1, nndesc2)
        redistributor.redistribute(H_sS, H_Ss)

        return H_Ss

    def scalapack_diagonalize(self, H_sS):

        mb = 32
        N = self.nS
        
        g1 = BlacsGrid(world, size,    1)
        g2 = BlacsGrid(world, size//2, 2)
        nndesc1 = g1.new_descriptor(N, N, self.nS_local,  N) 
        nndesc2 = g2.new_descriptor(N, N, mb, mb)
        
        A_ss = nndesc2.empty(dtype=H_sS.dtype)
        redistributor = Redistributor(world, nndesc1, nndesc2)
        redistributor.redistribute(H_sS, A_ss)
        
        # diagonalize
        v_ss = nndesc2.zeros(dtype=A_ss.dtype)
        w_S = np.zeros(N,dtype=float)
        nndesc2.diagonalize_dc(A_ss, v_ss, w_S, 'L')
        
        # distribute the eigenvectors to master
        v_sS = np.zeros_like(H_sS)
        redistributor = Redistributor(world, nndesc2, nndesc1)
        redistributor.redistribute(v_ss, v_sS)

#        v2_SS = np.zeros((self.nS, self.nS), dtype=complex)
#        world.all_gather(v_sS, v2_SS)
        
        return w_S, v_sS.conj()

    """
    def get_chi(self, w):
        H_SS = self.calculate()
        self.printtxt('Diagonalizing BSE matrix.')
        self.diagonalize(H_SS)
        self.printtxt('Calculating BSE response function.')

        w_S = self.w_S
        if not self.coupling:
            v_SS = self.v_SS.T # v_SS[:,lamda]
        else:
            v_SS = self.v_SS
        rhoG0_S = self.rhoG0_S
        focc_S = self.focc_S

        # get overlap matrix
        if self.coupling:
            tmp = np.dot(v_SS.conj().T, v_SS )
            overlap_SS = np.linalg.inv(tmp)
    
        # get chi
        chi_wGG = np.zeros((len(w), self.npw, self.npw), dtype=complex)
        t0 = time()

        A_S = np.dot(rhoG0_S, v_SS)
        B_S = np.dot(rhoG0_S*focc_S, v_SS)
        if self.coupling:
            C_S = np.dot(B_S.conj(), overlap_SS.T) * A_S
        else:
            C_S = B_S.conj() * A_S

        return chi_wGG
    """
