from __future__ import print_function
import sys
from time import ctime
import numpy as np

from ase.parallel import paropen
from ase.units import Ha
from ase.utils import devnull

from gpaw import GPAW
from gpaw.response.df0 import DF
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.mpi import rank, size, world
from gpaw.response.parallel import parallel_partition, \
     parallel_partition_list, set_communicator
from scipy.special.orthogonal import p_roots

class RPACorrelation:

    def __init__(self,
                 calc,
                 vcut=None,
                 txt=None,
                 tag=None,
                 qsym=True):
        
        self.calc = calc
        self.tag = tag
        
        if txt is None:
            if rank == 0:
                self.txt = sys.stdout
            else:
                sys.stdout = devnull
                self.txt = devnull
        else:
            assert type(txt) is str
            from ase.parallel import paropen
            self.txt = paropen(txt, 'w')

        self.qsym = qsym
        self.vcut = vcut
        self.nspins = calc.wfs.nspins
        self.bz_k_points = calc.wfs.kd.bzk_kc
        self.atoms = calc.get_atoms()
        self.setups = calc.wfs.setups
        self.bz_q_points = calc.wfs.kd.get_bz_q_points(first=True)
        if qsym == False:
            self.ibz_q_points = self.bz_q_points
            self.q_weights = (np.ones(len(self.bz_q_points))
                              / len(self.bz_q_points))
        else:
            op_scc = calc.wfs.kd.symmetry.op_scc
            self.ibz_q_points = calc.wfs.kd.get_ibz_q_points(self.bz_q_points,
                                                             op_scc)[0]
            self.q_weights = calc.wfs.kd.q_weights
        
        self.print_initialization()
        self.initialized = 0
   
    def get_rpa_correlation_energy(self,
                                   kcommsize=None,
                                   dfcommsize=world.size,
                                   directions=None,
                                   skip_gamma=False,
                                   ecut=10,
                                   nbands=None,
                                   gauss_legendre=None,
                                   frequency_cut=None,
                                   frequency_scale=None,
                                   w=None,
                                   restart=None):
            
        self.initialize_calculation(w,
                                    ecut,
                                    nbands,
                                    kcommsize,
                                    gauss_legendre,
                                    frequency_cut,
                                    frequency_scale)

        if dfcommsize == world.size:
            self.dfcomm = world
            E_q = []
            if restart is not None:
                assert type(restart) is str
                try:
                    f = paropen(restart, 'r')
                    lines = f.readlines()
                    for line in lines:
                        E_q.append(eval(line))
                    f.close()
                    print('Correlation energy obtained ' \
                          +'from %s q-points obtained from restart file: ' \
                          % len(E_q), restart, file=self.txt)
                    print(file=self.txt)
                except:
                    IOError
    
            for index, q in zip(range(len(E_q), len(self.ibz_q_points)),
                                self.ibz_q_points[len(E_q):]):
                if abs(np.dot(q, q))**0.5 < 1.e-5:
                    E_q0 = 0.
                    if skip_gamma:
                        print('Not calculating q at the Gamma point', file=self.txt)
                        print(file=self.txt)
                    else:
                        if directions is None:
                            directions = [[0, 1/3.], [1, 1/3.], [2, 1/3.]]
                        for d in directions:                                   
                            E_q0 += self.E_q(q,
                                             index=index,
                                             direction=d[0]) * d[1]
                    E_q.append(E_q0)
                else:
                    E_q.append(self.E_q(q, index=index))
                    
                if restart is not None:
                    f = paropen(restart, 'a')
                    print(E_q[-1], file=f)
                    f.close()
    
            E = np.dot(np.array(self.q_weights), np.array(E_q).real)

        else: # parallelzation over q points
            print('parallelization over q point ! ', file=self.txt)
            # creates q list
            qlist = []
            qweight = []
            id = 0
            for iq, q in enumerate(self.ibz_q_points):
                if abs(np.dot(q, q))**0.5 < 1.e-5:
                    if skip_gamma:
                        continue
                    else:
                        if directions is None:
                            directions = [[0, 1/3.], [1, 1/3.], [2, 1/3.]]
                        for d in directions:
                            qlist.append((id, q, d[0], d[1]))
                            qweight.append(self.q_weights[iq])
                            id += 1
                        continue
                qlist.append((id, q, 0, 1))
                qweight.append(self.q_weights[iq])
                id += 1
            nq = len(qlist)
    
            # distribute q list
            self.dfcomm, qcomm = set_communicator(world,
                                                  world.rank,
                                                  world.size,
                                                  kcommsize=dfcommsize)[:2]
            nq, nq_local, qlist_local = parallel_partition_list(nq,
                                                              qcomm.rank,
                                                              qcomm.size)
    
            E_q = np.zeros(nq)
            
            for iq in qlist_local:
                try:
                    ff = open('E_q_%s_%s.dat' %(self.tag,iq), 'r')
                    E_q[iq] = ff.readline().split()[-2]
                    print('Reading E_q[%s] '%(iq), E_q[iq], file=self.txt) 
                except:
                    E_q[iq] = self.E_q(qlist[iq][1],
                                   index=iq,
                                   direction=qlist[iq][2]) * qlist[iq][3]

                    if self.tag is not None and self.dfcomm.rank == 0:
                        ff = open('E_q_%s_%s.dat' %(self.tag,iq), 'a')
                        print(qlist[iq][1:4], E_q[iq], qweight[iq], file=ff)
                        ff.close()
                    
            qcomm.sum(E_q)

            print('(q, direction, weight), E_q, qweight', file=self.txt)
            for iq in range(nq):
                print(qlist[iq][1:4], E_q[iq], qweight[iq], file=self.txt)
                
            E = np.dot(np.array(qweight), np.array(E_q))

        print('RPA correlation energy:', file=self.txt)
        print('E_c = %s eV' % E, file=self.txt)
        print(file=self.txt)
        print('Calculation completed at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)
        return E

    def get_E_q(self,
                kcommsize=1,
                index=None,
                q=[0., 0., 0.],
                direction=0,
                integrated=True,
                ecut=10,
                nbands=None,
                gauss_legendre=None,
                frequency_cut=None,
                frequency_scale=None,
                w=None):

        self.initialize_calculation(w,
                                    ecut,
                                    nbands,
                                    kcommsize,
                                    gauss_legendre,
                                    frequency_cut,
                                    frequency_scale)
        self.dfcomm = world
        E_q = self.E_q(q,
                       direction=direction,
                       integrated=integrated)
        
        print('Calculation completed at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)

        return E_q


    def E_q(self,
            q,
            index=None,
            direction=0,
            integrated=True):

        if abs(np.dot(q, q))**0.5 < 1.e-5:
            q = [0.,0.,0.]
            q[direction] = 1.e-5
            optical_limit = True
        else:
            optical_limit = False

        dummy = DF(calc=self.calc,
                   eta=0.0,
                   w=self.w * 1j,
                   q=q,
                   ecut=self.ecut,
                   G_plus_q=True,
                   optical_limit=optical_limit,
                   hilbert_trans=False)
        dummy.txt = devnull
        dummy.initialize(simple_version=True)
        npw = dummy.npw
        del dummy

        if self.nbands is None:
            nbands = npw
        else:
            nbands = self.nbands

        if self.txt is sys.stdout:
            txt = 'response.txt'
        else:
            txt='response_'+self.txt.name
        df = DF(calc=self.calc,
                xc=None,
                nbands=nbands,
                eta=0.0,
                q=q,
                txt=txt,
                vcut=self.vcut,
                w=self.w * 1j,
                ecut=self.ecut,
                G_plus_q=True,
                kcommsize=self.kcommsize,
                comm=self.dfcomm,
                optical_limit=optical_limit,
                hilbert_trans=False)

        if index is None:
            print('Calculating KS response function at:', file=self.txt)
        else:
            print('#', index, \
                  '- Calculating KS response function at:', file=self.txt)
        if optical_limit:
            print('q = [0 0 0] -', 'Polarization: ', direction, file=self.txt)
        else:
            print('q = [%1.6f %1.6f %1.6f] -' \
                  % (q[0],q[1],q[2]), '%s planewaves' % npw, file=self.txt)

        e_wGG = df.get_dielectric_matrix(xc='RPA', overwritechi0=True)
        df.chi0_wGG = None
        Nw_local = len(e_wGG)
        local_E_q_w = np.zeros(Nw_local, dtype=complex)
        E_q_w = np.empty(len(self.w), complex)
        for i in range(Nw_local):
            local_E_q_w[i] = (np.log(np.linalg.det(e_wGG[i]))
                              + len(e_wGG[0]) - np.trace(e_wGG[i]))
            #local_E_q_w[i] = (np.sum(np.log(np.linalg.eigvals(e_wGG[i])))
            #                  + len(e_wGG[0]) - np.trace(e_wGG[i]))
        df.wcomm.all_gather(local_E_q_w, E_q_w)
        del df
        
        if self.gauss_legendre is not None:
            E_q = np.sum(E_q_w * self.gauss_weights * self.transform) \
                  / (4*np.pi)
        else:   
            dws = self.w[1:] - self.w[:-1]
            E_q = np.dot((E_q_w[:-1] + E_q_w[1:])/2., dws) / (2.*np.pi)


        print('E_c(q) = %s eV' % E_q.real, file=self.txt)
        print(file=self.txt)

        if integrated:
            return E_q.real
        else:
            return E_q_w.real               

    def initialize_calculation(self,
                               w,
                               ecut,
                               nbands,
                               kcommsize,
                               gauss_legendre,
                               frequency_cut,
                               frequency_scale):
        if kcommsize is None:
            if len(self.calc.wfs.kd.bzk_kc) == 1:
                kcommsize = 1
            else:
                kcommsize = world.size
            
        if w is not None:
            assert (gauss_legendre is None and
                    frequency_cut is None and
                    frequency_scale is None)
        else:
            if gauss_legendre is None:
                gauss_legendre = 16
            self.gauss_points, self.gauss_weights = p_roots(gauss_legendre)
            if frequency_scale is None:
                frequency_scale = 2.0
            if frequency_cut is None:
                frequency_cut = 800.
            ys = 0.5 - 0.5 * self.gauss_points
            ys = ys[::-1]
            w = (-np.log(1-ys))**frequency_scale
            w *= frequency_cut/w[-1]
            alpha = (-np.log(1-ys[-1]))**frequency_scale/frequency_cut
            transform = (-np.log(1-ys))**(frequency_scale-1) \
                        / (1-ys)*frequency_scale/alpha
            self.transform = transform
            
        dummy = DF(calc=self.calc,
                   eta=0.0,
                   w=w * 1j,
                   q=[0.,0.,0.0001],
                   ecut=ecut,
                   optical_limit=True,
                   hilbert_trans=False,
                   kcommsize=kcommsize)
        dummy.txt = devnull
        dummy.initialize(simple_version=True)

        self.npw = dummy.npw
        self.ecut = ecut
        self.w = w
        self.gauss_legendre = gauss_legendre
        self.frequency_cut = frequency_cut
        self.frequency_scale = frequency_scale
        self.kcommsize = kcommsize
        self.nbands = nbands

        print(file=self.txt)
        print('Planewave cutoff              : %s eV' % ecut, file=self.txt)
        print('Number of Planewaves at Gamma : %s' % self.npw, file=self.txt)
        if self.nbands is None:
            print('Response function bands       :'\
                  + ' Equal to number of Planewaves', file=self.txt)
        else:
            print('Response function bands       : %s' \
                  % self.nbands, file=self.txt)
        print('Frequencies', file=self.txt)
        if self.gauss_legendre is not None:
            print('    Gauss-Legendre integration '\
                  + 'with %s frequency points' % len(self.w), file=self.txt)
            print('    Frequency cutoff is '\
                  + '%s eV and scale (B) is %s' % (self.w[-1],
                                                  self.frequency_scale), file=self.txt)
        else:
            print('    %s specified frequency points' \
                  % len(self.w), file=self.txt)
            print('    Frequency cutoff is %s eV' \
                  % self.w[-1], file=self.txt)
        print(file=self.txt)
        print('Parallelization scheme', file=self.txt)
        print('     Total CPUs        : %d' % dummy.comm.size, file=self.txt)
        if dummy.kd.nbzkpts == 1:
            print('     Band parsize      : %d' % dummy.kcomm.size, file=self.txt)
        else:
            print('     Kpoint parsize    : %d' % dummy.kcomm.size, file=self.txt)
        print('     Frequency parsize : %d' % dummy.wScomm.size, file=self.txt)
        print('Memory usage estimate', file=self.txt)
        print('     chi0_wGG(Q)       : %f M / cpu' \
              % (dummy.Nw_local * self.npw**2 * 16. / 1024**2), file=self.txt)
        print(file=self.txt)
        del dummy

    def print_initialization(self):
        print('------------------------------------------------------', file=self.txt)
        print('Non-self-consistent RPA correlation energy', file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print('Started at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('Atoms                          :   %s' \
              % self.atoms.get_chemical_formula(mode="hill"), file=self.txt)
        print('Ground state XC functional     :   %s' \
              % self.calc.hamiltonian.xc.name, file=self.txt)
        print('Valence electrons              :   %s' \
              % self.setups.nvalence, file=self.txt)
        print('Number of Bands                :   %s' \
              % self.calc.wfs.bd.nbands, file=self.txt)
        print('Number of Converged Bands      :   %s' \
              % self.calc.input_parameters['convergence']['bands'], file=self.txt)
        print('Number of Spins                :   %s' \
              % self.nspins, file=self.txt)
        print('Number of k-points             :   %s' \
              % len(self.calc.wfs.kd.bzk_kc), file=self.txt)
        print('Number of q-points             :   %s' \
              % len(self.bz_q_points), file=self.txt)
        print('Number of Irreducible k-points :   %s' \
              % len(self.calc.wfs.kd.ibzk_kc), file=self.txt)
        if self.qsym:
            print('Number of Irreducible q-points :   %s' \
                  % len(self.ibz_q_points), file=self.txt)
        else:
            print('No reduction of q-points', file=self.txt) 
        print(file=self.txt)
        for q, weight in zip(self.ibz_q_points, self.q_weights):
            print('q: [%1.4f %1.4f %1.4f] - weight: %1.3f' \
                  % (q[0],q[1],q[2], weight), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)


    def get_C6_coefficient(self,
                           ecut=100.,
                           nbands=None,
                           kcommsize=None,
                           gauss_legendre=None,
                           frequency_cut=None,
                           frequency_scale=None,
                           direction=2):

        self.initialize_calculation(None,
                                    ecut,
                                    nbands,
                                    kcommsize,
                                    gauss_legendre,
                                    frequency_cut,
                                    frequency_scale)

        d = direction
        d_pro = []
        for i in range(3):
            if i != d:
                d_pro.append(i)
        
        dummy = DF(calc=self.calc,
                   eta=0.0,
                   w=self.w * 1j,
                   ecut=self.ecut,
                   hilbert_trans=False)
        dummy.txt = devnull
        dummy.initialize(simple_version=True)
        npw = dummy.npw
        del dummy

        q = [0.,0.,0.]
        q[d] = 1.e-5

        if self.nbands is None:
            nbands = npw
        else:
            nbands = self.nbands

        if self.txt is sys.stdout:
            txt = 'response.txt'
        else:
            txt='response_'+self.txt.name
        df = DF(calc=self.calc,
                xc=None,
                nbands=nbands,
                eta=0.0,
                q=q,
                txt=txt,
                vcut=self.vcut,
                w=self.w * 1j,
                ecut=self.ecut,
                comm=world,
                optical_limit=True,
                G_plus_q=True,
                kcommsize=self.kcommsize,
                hilbert_trans=False)
        
        print('Calculating RPA response function', file=self.txt)
        print('Polarization: %s' % d, file=self.txt)

        chi_wGG = df.get_chi(xc='RPA')
        chi0_wGG = df.chi0_wGG

        Nw_local = len(chi_wGG)
        local_a0_w = np.zeros(Nw_local, dtype=complex)
        a0_w = np.empty(len(self.w), complex)
        local_a_w = np.zeros(Nw_local, dtype=complex)
        a_w = np.empty(len(self.w), complex)

        Gvec_Gv = np.dot(df.Gvec_Gc + np.array(q), df.bcell_cv)
        gd = self.calc.density.gd
        n_d = gd.get_size_of_global_array()[d]
        d_d = gd.get_grid_spacings()[d]
        r_d = np.array([i*d_d for i in range(n_d)])

        print('Calculating real space integrals', file=self.txt)

        int_G = np.zeros(npw, complex)
        for iG in range(npw):
            if df.Gvec_Gc[iG, d_pro[0]] == 0 and df.Gvec_Gc[iG, d_pro[1]] == 0:
                int_G[iG] = np.sum(r_d * np.exp(1j*Gvec_Gv[iG, d] * r_d))*d_d
        int2_GG = np.outer(int_G, int_G.conj())

        print('Calculating dynamic polarizability', file=self.txt)

        for i in range(Nw_local):
            local_a0_w[i] = np.trace(np.dot(chi0_wGG[i], int2_GG))
            local_a_w[i] = np.trace(np.dot(chi_wGG[i], int2_GG))
        df.wcomm.all_gather(local_a0_w, a0_w)
        df.wcomm.all_gather(local_a_w, a_w)

        A = df.vol / gd.cell_cv[d,d]
        a0_w *= A**2 / df.vol
        a_w *= A**2 / df.vol

        del df
        
        C06 = np.sum(a0_w**2 * self.gauss_weights
                     * self.transform) * 3 / (2*np.pi)
        C6 = np.sum(a_w**2 * self.gauss_weights
                    * self.transform) * 3 / (2*np.pi)

        print('C06 = %s Ha*Bohr**6' % (C06.real / Ha), file=self.txt)
        print('C6 = %s Ha*Bohr**6' % (C6.real / Ha), file=self.txt)
        print(file=self.txt)

        return C6.real / Ha, C06.real / Ha
