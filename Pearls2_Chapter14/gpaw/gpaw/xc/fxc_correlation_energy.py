from __future__ import print_function
import sys
from time import ctime, time
import numpy as np

from ase.parallel import paropen
from ase.units import Hartree, Bohr
from ase.utils import devnull

from gpaw import GPAW
from gpaw.xc import XC
from gpaw.xc.libxc import LibXC
from gpaw.response.df0 import DF
from gpaw.utilities.blas import gemmdot, axpy
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.mpi import rank, size, world, serial_comm
from gpaw.blacs import BlacsGrid, BlacsDescriptor, Redistributor
from gpaw.response.parallel import parallel_partition, set_communicator
from gpaw.fd_operators import Gradient, Laplace
from gpaw.sphere.lebedev import weight_n, R_nv
from gpaw.io.tar import Writer, Reader
from gpaw import debug
from scipy.special import sici
from scipy.special.orthogonal import p_roots

class FXCCorrelation:

    def __init__(self,
                 calc,
                 txt=None,
                 tag='',
                 qsym=True,
                 xc=None,
                 num=0,
                 write=False,
                 method='standard',
                 lambda_points=8,
                 density_cut=None,
                 paw_correction=1,
                 unit_cells=None):
        
        self.calc = calc

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
        
        self.tag = tag
        self.qsym = qsym
        self.num = num
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
        
        if xc is None:
            self.xc = 'RPA'
        else:
            self.xc = xc

        self.lambda_points = lambda_points
        self.method = method
        self.density_cut = density_cut
        if self.density_cut is None:
            self.density_cut = 1.e-6
        assert paw_correction in range(4)
        self.paw_correction = paw_correction
        self.unit_cells = unit_cells
        self.print_initialization()
        self.initialized = 0

   
    def get_fxc_correlation_energy(self,
                                   kcommsize=None,
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
        
        if directions is None:
            directions = [[0, 1/3.], [1, 1/3.], [2, 1/3.]]
        if skip_gamma:
            directions = []
                
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

        if len(E_q) == 0:
            kernel = Kernel(self.calc,
                            self.xc,
                            self.method,
                            self.ibz_q_points,
                            self.w,
                            self.ecut,
                            self.txt,
                            self.tag,
                            paw_correction=self.paw_correction,
                            unit_cells=self.unit_cells,
                            density_cut=self.density_cut)
            kernel.calculate_fhxc(directions=[d[0] for d in directions])
            del kernel
 
        for index, q in zip(range(len(E_q), len(self.ibz_q_points)),
                            self.ibz_q_points[len(E_q):]):
            if abs(np.dot(q, q))**0.5 < 1.e-5:
                E_q0 = 0.
                if skip_gamma:
                    print('Not calculating q at the Gamma point', file=self.txt)
                    print(file=self.txt)
                else:
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

        print('%s correlation energy:' % self.xc, file=self.txt)
        print('E_c = %s eV' % E, file=self.txt)
        print(file=self.txt)
        print('Calculation completed at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)

        return E


    def get_E_q(self,
                kcommsize=None,
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

        ns = self.nspins

        if self.xc[0] == 'r':
            name = self.xc + '_' + self.method + '_' + self.tag
            if optical_limit:
                r = Reader('fhxc_%s_%s_0%s.gpw' % (name, self.ecut, direction))
            else:
                r = Reader('fhxc_%s_%s_%s.gpw' % (name, self.ecut, index))
        else:
            name = self.xc + '_' + self.tag
            if optical_limit:
                r = Reader('fhxc_%s_%s_0%s.gpw' % (name, self.ecut, direction))
            else:
                r = Reader('fhxc_%s_%s_%s.gpw' % (name, self.ecut, index))
                
        npw = r.dimension('sG') / ns

        if index is not None:
            print('#', index, file=self.txt)
        if optical_limit:
            print('q = [0 0 0] -', 'Polarization: ', direction, file=self.txt)
        else:
            print('q = [%1.6f %1.6f %1.6f] -' \
                % (q[0],q[1],q[2]), '%s planewaves' % npw, file=self.txt)
        print('Calculating KS response function', file=self.txt)

        if self.nbands is None:
            nbands = npw
        else:
            nbands = self.nbands

        if self.txt is sys.stdout:
            txt = 'response.txt'
        else:
            txt='response_' + self.txt.name
        df = DF(calc=self.calc,
                xc=None,
                nbands=nbands,
                eta=0.0,
                q=q,
                txt=txt,
                w=self.w * 1j,
                ecut=self.ecut,
                G_plus_q=True,
                density_cut=self.density_cut,
                kcommsize=self.kcommsize,
                comm=world,
                optical_limit=optical_limit,
                hilbert_trans=False)
        
        df.initialize()
        Nw_local = df.Nw_local

        Kc_G = np.zeros(npw, dtype=complex)

        for iG in range(npw):
            qG = np.dot(q + df.Gvec_Gc[iG], df.bcell_cv)
            Kc_G[iG] = 4 * np.pi / np.dot(qG, qG)

        chi0 = np.zeros((Nw_local, ns*npw, ns*npw), dtype=complex)

        df.calculate(seperate_spin=0)
        chi0[:, :npw, :npw] = df.chi0_wGG[:]
        if ns == 2:
            df.ecut *= Hartree
            df.xc = 'RPA'
            df.initialize()
            df.calculate(seperate_spin=1)
            chi0[:, npw:, npw:] = df.chi0_wGG[:]
        del df.chi0_wGG
                
        local_E_q_w = np.zeros(Nw_local, dtype=complex)
        E_q_w = np.empty(len(self.w), complex)

        print('Calculating interacting response function', file=self.txt)
        ls, l_ws = p_roots(self.lambda_points)
        ls = (ls + 1.0) * 0.5
        l_ws *= 0.5
        for i in range(Nw_local):
            chi0_fhxc = np.dot(chi0[i], r.get('fhxc_sGsG'))
            for l, l_w in zip(ls, l_ws):
                chi_l = np.linalg.solve(np.eye(npw*ns, npw*ns)
                                        - l*chi0_fhxc, chi0[i]).real
                for s1 in range(ns):
                    for s2 in range(ns):
                        X_ss = chi_l[s1*npw:(s1+1)*npw, s2*npw:(s2+1)*npw]
                        local_E_q_w[i] -= np.dot(np.diag(X_ss), Kc_G)*l_w
            local_E_q_w[i] += np.dot(np.diag(chi0[i]), np.tile(Kc_G, ns))
        df.wcomm.all_gather(local_E_q_w, E_q_w)

        del df, chi0, chi0_fhxc, chi_l, X_ss, Kc_G
        
        if self.gauss_legendre is not None:
            E_q = np.sum(E_q_w * self.gauss_weights * self.transform) \
                  / (4 * np.pi)
        else:
            dws = self.w[1:] - self.w[:-1]
            E_q = np.dot((E_q_w[:-1] + E_q_w[1:])/2., dws) / (2 * np.pi)

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
                   xc='RPA',
                   eta=0.0,
                   w=w * 1j,
                   q=[0.,0.,0.0001],
                   ecut=ecut,
                   optical_limit=True,
                   hilbert_trans=False,
                   kcommsize=kcommsize)
        dummy.txt = devnull
        dummy.initialize(simple_version=True)
            
        self.ecut = ecut
        self.w = w
        self.gauss_legendre = gauss_legendre
        self.frequency_cut = frequency_cut
        self.frequency_scale = frequency_scale
        self.kcommsize = kcommsize
        self.nbands = nbands

        if debug:
            print(file=self.txt)
            print('DEBUG MODE', file=self.txt)
            print(file=self.txt)
    
        print('Numerical coupling constant integration' \
              + ' with % s Gauss-Legendre points' % self.lambda_points, file=self.txt)
        print(file=self.txt)
        print('Planewave cutoff              : %s eV' % ecut, file=self.txt)
        print('Number of Planewaves at Gamma : %s' % dummy.npw, file=self.txt)
        if self.nbands is None:
            print('Response function bands       :' \
                  + ' Equal to number of Planewaves', file=self.txt)
        elif type(self.nbands) is float:
            print('Response function bands       : %s' \
                  % int(dummy.npw * self.nbands), file=self.txt)
        else:
            print('Response function bands       : %s' \
                  % self.nbands, file=self.txt)
        if self.density_cut is not None:
            print(file=self.txt)
            print('Min value of pseudo density   : %1.2e Bohr^-3'\
                  % np.min(self.calc.density.nt_sG), file=self.txt)
            print('Max value of pseudo density   : %1.2e Bohr^-3'\
                  % np.max(self.calc.density.nt_sG), file=self.txt)
            print('Density cutoff in fxc at      : %1.2e Bohr^-3'\
                  % self.density_cut, file=self.txt)
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
        print('     chi0_wGG(q)         : %6.3f MB / frequency point' \
              % (self.nspins**2 * dummy.npw**2 * 16. / 1024**2), file=self.txt)
        print('     Kernel_GG(q)        : %6.3f MB / CPU' \
              % ((2*self.nspins + 3%self.nspins) * dummy.npw**2 * 16. / 1024**2), file=self.txt)
        if self.method == 'standard':
            n0 = dummy.gd.N_c[0] * dummy.gd.N_c[1] * dummy.gd.N_c[2]
            print('     Kernel_rG(q) (Int.) : %6.3f MB / CPU' \
                      % ((2*self.nspins + 3%self.nspins)
                         * dummy.npw * float(n0)/world.size * 16. / 1024**2), file=self.txt)
            
        print(file=self.txt)
        del dummy


    def print_initialization(self):
        
        print('------------------------------------------------------', file=self.txt)
        print('Non-self-consistent %s correlation energy' \
              % self.xc, file=self.txt)
        if self.xc is not 'RPA':
            if self.xc[0] == 'r':
                if self.method == 'solid':
                    print('Periodic average density', file=self.txt)
                else:
                    print('Non-periodic two-point density', file=self.txt)
            if self.paw_correction == 0:
                print('Using pseudo-density with ALDA PAW corrections', file=self.txt)
            elif self.paw_correction == 1:
                print('Using all-electron density', file=self.txt)
            elif self.paw_correction == 2:
                print('Using pseudo-density with average ALDA PAW corrections', file=self.txt)
            elif self.paw_correction == 3:
                print('Using pseudo-density', file=self.txt)
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
        
        q = [0.,0.,0.]
        q[d] = 1.e-5

        kernel = Kernel(self.calc,
                        self.xc,
                        self.method,
                        self.ibz_q_points,
                        self.w,
                        self.ecut,
                        self.txt,
                        self.tag,
                        paw_correction=self.paw_correction,
                        unit_cells=self.unit_cells,
                        density_cut=self.density_cut)

        kernel.calculate_fhxc(directions=[d])

        del kernel

        if self.xc[0] == 'r':
            name = self.xc + '_' + self.method + '_' + self.tag
            r = Reader('fhxc_%s_%s_0%s.gpw' % (name, self.ecut, direction))
        else:
            name = self.xc + '_' + self.tag
            r = Reader('fhxc_%s_%s_0%s.gpw' % (name, self.ecut, direction))
                
        ns = self.nspins
        npw = r.dimension('sG') / ns

        if self.nbands is None:
            nbands = npw
        else:
            nbands = self.nbands

        if self.txt is sys.stdout:
            txt = 'response.txt'
        else:
            txt='response_' + self.txt.name
        df = DF(calc=self.calc,
                xc=None,
                nbands=nbands,
                eta=0.0,
                q=q,
                txt=txt,
                w=self.w * 1j,
                ecut=self.ecut,
                comm=world,
                optical_limit=True,
                G_plus_q=True,
                kcommsize=self.kcommsize,
                hilbert_trans=False)
        
        print('Calculating %s response function' % self.xc, file=self.txt)
        print('Polarization: %s' % d, file=self.txt)

        df.initialize()
        Nw_local = df.Nw_local
        chi0 = np.zeros((Nw_local, ns*npw, ns*npw), dtype=complex)

        df.calculate(seperate_spin=0)
        chi0[:, :npw, :npw] = df.chi0_wGG[:]
        if ns == 2:
            df.ecut *= Hartree
            df.xc = 'RPA'
            df.initialize()
            df.calculate(seperate_spin=1)
            chi0[:, npw:2*npw, npw:2*npw] = df.chi0_wGG[:]
        del df.chi0_wGG
        
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
            chi0_fhxc = np.dot(chi0[i], r.get('fhxc_sGsG'))
            chi = np.linalg.solve(np.eye(npw*ns, npw*ns)
                                  - chi0_fhxc, chi0[i]).real
            for s1 in range(ns):
                X0 = chi0[i, s1*npw:(s1+1)*npw, s1*npw:(s1+1)*npw]
                local_a0_w[i] += np.trace(np.dot(X0, int2_GG))
                for s2 in range(ns):
                    X_ss = chi[s1*npw:(s1+1)*npw, s2*npw:(s2+1)*npw]
                    local_a_w[i] += np.trace(np.dot(X_ss, int2_GG))
        df.wcomm.all_gather(local_a0_w, a0_w)
        df.wcomm.all_gather(local_a_w, a_w)

        A = df.vol / gd.cell_cv[d,d]
        a0_w *= A**2 / df.vol
        a_w *= A**2 / df.vol

        C06 = np.sum(a0_w**2 * self.gauss_weights
                     * self.transform) * 3 / (2*np.pi)
        C6 = np.sum(a_w**2 * self.gauss_weights
                    * self.transform) * 3 / (2*np.pi)

        print('C06 = %s Ha*Bohr**6' % (C06.real / Hartree), file=self.txt)
        print('C6 = %s Ha*Bohr**6' % (C6.real / Hartree), file=self.txt)
        print(file=self.txt)

        return C6.real / Hartree, C06.real / Hartree


class Kernel:
    
    def __init__(self, calc, xc, method, q_points, w, ecut, txt, tag,
                 paw_correction=1, unit_cells=None, density_cut=None):

        self.calc = calc
        self.xc = xc
        self.method = method
        self.q_points = q_points
        self.ecut = ecut
        self.txt = txt
        self.tag = tag
        self.paw_correction = paw_correction
        self.unit_cells = unit_cells
        self.density_cut = density_cut
        self.Gvec_qGc = {}
        self.npw_q = []

        self.A_x = -(3/4.) * (3/np.pi)**(1/3.)

        if calc.wfs.nspins == 1:
            self.n_g = calc.get_all_electron_density(gridrefinement=1)
        else:
            self.n_g = (calc.get_all_electron_density(gridrefinement=1, spin=0)
                     + calc.get_all_electron_density(gridrefinement=1, spin=1))
        self.n_g *= Bohr**3
        
        if xc == 'rAPBE':
            if calc.wfs.nspins == 1:
                nf_g = calc.get_all_electron_density(gridrefinement=2)
            else:
                nf_g = (calc.get_all_electron_density(gridrefinement=2, spin=0)
                     + calc.get_all_electron_density(gridrefinement=2, spin=1))
            nf_g *= Bohr**3
            gdf = calc.density.gd.refine()
            grad_v = [Gradient(gdf, v, n=1).apply for v in range(3)]
            #laplace = Laplace(gdf, n=1).apply
            gradnf_vg = gdf.empty(3)
            for v in range(3):
                grad_v[v](nf_g, gradnf_vg[v])
            self.gradn_vg = gradnf_vg[:,::2,::2,::2]
            
        for iq, q in enumerate(q_points):
            dummy = DF(calc=self.calc,
                       q=q,
                       hilbert_trans=False,
                       w=1j*w,
                       ecut=ecut,
                       G_plus_q=True)
            dummy.txt = devnull
            dummy.initialize(simple_version=True)
            self.Gvec_qGc[iq] = dummy.Gvec_Gc
            self.npw_q.append(len(dummy.Gvec_Gc))
        self.gd = dummy.gd
        self.nG = dummy.gd.N_c
        self.vol = dummy.vol
        self.bcell_cv = dummy.bcell_cv
        self.acell_cv = dummy.acell_cv
        del dummy

    def calculate_fhxc(self, directions=[0,1,2]):

        print('Calculating %s kernel' % self.xc, file=self.txt)

        if self.xc == 'RPA':
            self.calculate_rpa_kernel(directions)
            
        else:
            if self.xc[0] == 'r':
                if self.method == 'solid':
                    self.calculate_rkernel_solid(directions)
                elif self.method == 'standard':
                    self.calculate_rkernel(directions)
                else:
                    raise 'Method % s not known' % self.method
            else:
                self.calculate_local_kernel()


    def calculate_rpa_kernel(self, directions=[0,1,2]):

        Gvec_qGc = self.Gvec_qGc
        npw_q = self.npw_q
        bcell_cv = self.bcell_cv

        ns = self.calc.wfs.nspins
        
        for iq, q in enumerate(self.q_points):
            if abs(np.dot(q, q))**0.5 < 1.e-5:
                iq0 = iq
                continue
            Gvec_Gc = Gvec_qGc[iq]
            npw = npw_q[iq]
            qG_Gv = np.dot(q + Gvec_Gc, bcell_cv)

            Kc_GG = np.zeros((npw, npw), dtype=float)
            for iG in range(npw):
                Kc_GG[iG,iG] = 4 * np.pi / np.dot(qG_Gv[iG], qG_Gv[iG])

            fhxc_sGsG = np.tile(Kc_GG, (ns, ns))
            
            if rank == 0:
                w = Writer('fhxc_RPA_%s_%s_%s.gpw' % (self.tag, self.ecut, iq))
                w.dimension('sG', ns*npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG)
                w.close()
            world.barrier()

        # Gamma term
        Gvec_Gc = Gvec_qGc[iq0]
        npw = npw_q[iq0]
        for d in directions:
            q = np.array([0., 0., 0.])
            q[d] = 1.e-5
            qG_Gv = np.dot(q + Gvec_Gc, bcell_cv)
            Kc_GG = np.zeros((npw, npw), dtype=float)
            for iG in range(npw):
                Kc_GG[iG,iG] = 4 * np.pi / np.dot(qG_Gv[iG], qG_Gv[iG])
            
            fhxc_sGsG = np.tile(Kc_GG, (ns, ns))

            if rank == 0:
                w = Writer('fhxc_RPA_%s_%s_0%s.gpw' % (self.tag, self.ecut, d))
                w.dimension('sG', ns*npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG)
                w.close()
            world.barrier()
        
        print(file=self.txt)


    def calculate_rkernel(self, directions=[0,1,2]):

        Gvec_qGc = self.Gvec_qGc
        npw_q = self.npw_q
        gd = self.gd
        nG = self.nG
        vol = self.vol
        bcell_cv = self.bcell_cv
        acell_cv = self.acell_cv
        
        ns = self.calc.wfs.nspins
        n_g = self.n_g

        fx_g = ns * self.get_fxc_g(n_g)
        qc_g = (-4 * np.pi * ns / fx_g)**0.5
        flocal_g = qc_g**3 * fx_g / (6 * np.pi**2)
        #flocal_g = 4 * n_g * fx_g # LDA
        Vlocal_g = 2 * qc_g / np.pi
        #Vlocal_g = 4 * (3 * n_g / np.pi)**(1./3.)
        
        nG0 = np.prod(nG)
        r_vg = gd.get_grid_point_coordinates()
        r_vgx = r_vg[0].flatten()
        r_vgy = r_vg[1].flatten()
        r_vgz = r_vg[2].flatten()

        # Unit cells
        R = []
        R_weight = []
        if self.unit_cells is None:
            N_R = self.calc.wfs.kd.N_c
        else:
            N_R = self.unit_cells
        N_R0 = N_R[0]*N_R[1]*N_R[2]
        for i in range(-N_R[0]+1, N_R[0]):
            for j in range(-N_R[1]+1, N_R[1]):
                for h in range(-N_R[2]+1, N_R[2]):
                    R.append(i*acell_cv[0] + j*acell_cv[1] + h*acell_cv[2])
                    R_weight.append((N_R[0]-abs(i))*
                                    (N_R[1]-abs(j))*
                                    (N_R[2]-abs(h)) / float(N_R0))
        if N_R0 > 1:
            dv = self.calc.density.gd.dv
            gc = (3.*dv/4./np.pi)**(1./3.)
            Vlocal_g -= 2*np.pi * gc**2 /dv
            print('    Lattice point sampling: ' \
                  + '(%s x %s x %s)^2 ' % (N_R[0], N_R[1], N_R[2]) \
                  + ' Reduced to %s lattice points' % len(R), file=self.txt)
        
        l_g_size = -(-nG0 // world.size)
        l_g_range = range(world.rank * l_g_size,
                          min((world.rank+1) * l_g_size, nG0))

        fhxc_qsGr = {}
        for iq in range(len(npw_q)):
            fhxc_qsGr[iq] = np.zeros((ns, npw_q[iq], len(l_g_range)),
                                     dtype=complex)

        inv_error = np.seterr()
        np.seterr(invalid='ignore')
        np.seterr(divide='ignore')

        t0 = time()
        # Loop over Lattice points
        for i, R_i in enumerate(R):
            # Loop over r'.
            # f_rr, V_rr and V_off are functions of r (dim. as r_vg[0])
            if i == 1:
                print('      Finished 1 cell in %s seconds - estimated %s seconds left' % (int(time() - t0), int(len(R)*(time()-t0))), file=self.txt)
            if len(R) > 5:
                if (i+1) % (len(R)/5+1) == 0:
                    print('      Finished %s cells in %s seconds - estimated %s seconds left' % (i, int(time() - t0), int((len(R)-i)*(time() - t0)/i)), file=self.txt)
            for g in l_g_range:
                r_x = r_vgx[g] + R_i[0]
                r_y = r_vgy[g] + R_i[1]
                r_z = r_vgz[g] + R_i[2]

                # |r-r'-R_i|
                rr = ((r_vg[0]-r_x)**2 +
                      (r_vg[1]-r_y)**2 +
                      (r_vg[2]-r_z)**2)**0.5

                n_av = (n_g + n_g.flatten()[g]) / 2.
                fx_g = ns * self.get_fxc_g(n_av, index=g)
                qc_g = (-4 * np.pi * ns / fx_g)**0.5
                x = qc_g * rr
                #kf_g = (3 * np.pi**2 * n_av)**(1./3.)
                #x = 2 * kf_g * rr  # LDA
                f_rr = fx_g * (np.sin(x) - x*np.cos(x)) \
                              / (2 * np.pi**2 * rr**3)
                if N_R0 > 1:
                    V_rr = (sici(x)[0] * 2 / np.pi - 1) / rr
                else:
                    V_rr = (sici(x)[0] * 2 / np.pi) / rr
                    
                # Terms with r = r'
                if (np.abs(R_i) < 0.001).all():
                    tmp_flat = f_rr.flatten()
                    tmp_flat[g] = flocal_g.flatten()[g]
                    f_rr = tmp_flat.reshape(nG)
                    tmp_flat = V_rr.flatten()
                    tmp_flat[g] = Vlocal_g.flatten()[g]
                    V_rr = tmp_flat.reshape(nG)
                    del tmp_flat

                f_rr[np.where(n_av < self.density_cut)] = 0.0
                V_rr[np.where(n_av < self.density_cut)] = 0.0

                f_rr *= R_weight[i]
                V_rr *= R_weight[i]

                # r-r'-R_i
                r_r = np.array([r_vg[0]-r_x, r_vg[1]-r_y, r_vg[2]-r_z])

                # Fourier transform of r
                for iq, q in enumerate(self.q_points):
                    q_v = np.dot(q, bcell_cv)
                    e_q = np.exp(-1j * gemmdot(q_v, r_r, beta=0.0))
                    tmp_fhxc = np.fft.fftn((f_rr+V_rr)*e_q) * vol / nG0
                    if ns == 2:
                        tmp_V_off = np.fft.fftn(V_rr*e_q) * vol / nG0
                    for iG in range(npw_q[iq]):
                        assert (nG / 2 - np.abs(Gvec_qGc[iq][iG]) > 0).all()
                        f_i = Gvec_qGc[iq][iG] % nG
                        fhxc_qsGr[iq][0, iG, g-l_g_range[0]] += \
                                      tmp_fhxc[f_i[0], f_i[1], f_i[2]]
                        if ns == 2:
                            fhxc_qsGr[iq][1, iG, g-l_g_range[0]] += \
                                          tmp_V_off[f_i[0], f_i[1], f_i[2]]
        world.barrier()
        
        np.seterr(**inv_error)

        for iq, q in enumerate(self.q_points):
            npw = npw_q[iq]
            Gvec_Gc = Gvec_qGc[iq]
            fhxc_sGsG = np.zeros((ns*npw, ns*npw), complex)
            l_pw_size = -(-npw // world.size)
            l_pw_range = range(world.rank * l_pw_size,
                               min((world.rank+1) * l_pw_size, npw))

            if world.size > 1 :
                bg1 = BlacsGrid(world, 1, world.size)
                bg2 = BlacsGrid(world, world.size, 1)
                bd1 = bg1.new_descriptor(npw, nG0, npw, - (-nG0 / world.size))
                bd2 = bg2.new_descriptor(npw, nG0, -(-npw / world.size), nG0)

                fhxc_Glr = np.zeros((len(l_pw_range), nG0), dtype=complex)
                if ns == 2:
                    Koff_Glr = np.zeros((len(l_pw_range), nG0), dtype=complex)

                r = Redistributor(bg1.comm, bd1, bd2)
                r.redistribute(fhxc_qsGr[iq][0], fhxc_Glr, npw, nG0)
                if ns == 2:
                    r.redistribute(fhxc_qsGr[iq][1], Koff_Glr, npw, nG0)
            else:
                fhxc_Glr = fhxc_qsGr[iq][0]
                if ns == 2:
                    Koff_Glr = fhxc_qsGr[iq][1]

            # Fourier transform of r'
            for iG in range(len(l_pw_range)):
                tmp_fhxc = np.fft.fftn(fhxc_Glr[iG].reshape(nG)) * vol/nG0
                if ns == 2:
                    tmp_Koff = np.fft.fftn(Koff_Glr[iG].reshape(nG)) * vol/nG0
                for jG in range(npw):
                    assert (nG / 2 - np.abs(Gvec_Gc[jG]) > 0).all()
                    f_i = -Gvec_Gc[jG] % nG
                    fhxc_sGsG[l_pw_range[0] + iG, jG] = \
                                    tmp_fhxc[f_i[0], f_i[1], f_i[2]]
                    if ns == 2:
                        fhxc_sGsG[npw + l_pw_range[0] + iG, jG] += \
                                      tmp_Koff[f_i[0], f_i[1], f_i[2]]

            if ns == 2:
                fhxc_sGsG[:npw, npw:] = fhxc_sGsG[npw:, :npw]
                fhxc_sGsG[npw:, npw:] = fhxc_sGsG[:npw, :npw]

            world.sum(fhxc_sGsG)

            if rank == 0:
                if abs(np.dot(q, q))**0.5 < 1.e-5:
                    for d in directions:
                        q0 = np.array([0., 0., 0.])
                        q0[d] = 1.e-5

                        w = Writer('fhxc_%s_%s_%s_%s_0%s.gpw' % (self.xc,
                                                                 self.method,
                                                                 self.tag,
                                                                 self.ecut,
                                                                 d))
                        w.dimension('sG', ns*npw)
                        w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                        if N_R0 > 1:
                            Gq_Gv = np.dot(Gvec_Gc + q0, bcell_cv)
                            Gq_G = (Gq_Gv[:,0]**2 + Gq_Gv[:,1]**2 + Gq_Gv[:,2]**2)**0.5
                            vq_G = 4 * np.pi / Gq_G**2
                            w.fill(fhxc_sGsG / vol
                                   + np.tile(np.eye(npw) * vq_G, (ns, ns)))
                        else:
                            w.fill(fhxc_sGsG / vol)
                        w.close()
                else:
                    w = Writer('fhxc_%s_%s_%s_%s_%s.gpw' % (self.xc,
                                                            self.method,
                                                            self.tag,
                                                            self.ecut,
                                                            iq))
                    w.dimension('sG', ns*npw)
                    w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                    if N_R0 > 1:
                        Gq_Gv = np.dot(Gvec_Gc + q, bcell_cv)
                        Gq_G = (Gq_Gv[:,0]**2 + Gq_Gv[:,1]**2 + Gq_Gv[:,2]**2)**0.5
                        vq_G = 4 * np.pi / Gq_G**2
                        w.fill(fhxc_sGsG / vol
                               + np.tile(np.eye(npw) * vq_G, (ns, ns)))
                    else:
                        w.fill(fhxc_sGsG / vol)
                    w.close()
            world.barrier()

        print(file=self.txt)

        
    def calculate_rkernel_solid(self, directions=[0,1,2]):

        Gvec_qGc = self.Gvec_qGc
        npw_q = self.npw_q
        gd = self.gd
        nG = self.nG
        vol = self.vol
        bcell_cv = self.bcell_cv

        ns = self.calc.wfs.nspins

        n_g = self.n_g

        fx_g = ns * self.get_fxc_g(n_g)
        fx_g[np.where(n_g < self.density_cut)] = 0.0
        qc_g = (-4 * np.pi * ns / fx_g)**0.5
        #kf_g = (3 * np.pi**2 * n_g)**(1./3.)

        r_vg = gd.get_grid_point_coordinates()

        for iq, q in enumerate(self.q_points):
            if abs(np.dot(q, q))**0.5 < 1.e-5:
                iq0 = iq
                continue
            npw = npw_q[iq]
            Gvec_Gc = Gvec_qGc[iq]
            fhxc_sGsG = np.zeros((ns*npw, ns*npw), dtype=complex)
            l_pw_size = -(-npw // world.size)
            l_pw_range = range(world.rank * l_pw_size,
                               min((world.rank+1) * l_pw_size, npw))

            for iG in l_pw_range:
                for jG in range(npw):
                    dGq_c = (Gvec_Gc[iG] + Gvec_Gc[jG])/ 2. + q
                    if (np.abs(dGq_c) < 1.e-12).all():
                        dGq_c = np.array([0.0, 0.0, 0.0001])
                    dGq_v = np.dot(dGq_c, bcell_cv)
                    fx = fx_g.copy()
                    #fx[np.where(2*kf_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    fx[np.where(qc_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    v_c = 4 * np.pi * np.ones(np.shape(fx_g), float)
                    v_c /= np.dot(dGq_v, dGq_v)
                    #v_c[np.where(2*kf_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    v_c[np.where(qc_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0

                    dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                    dG_v = np.dot(dG_c, bcell_cv)
                    dGr_g = gemmdot(dG_v, r_vg, beta=0.0)

                    fhxc_sGsG[iG, jG] = gd.integrate(np.exp(-1j*dGr_g) \
                                                     * (fx + v_c))
                    if ns == 2:
                        fhxc_sGsG[npw + iG, jG] = \
                                      gd.integrate(np.exp(-1j*dGr_g) * v_c)
            if ns == 2:
                fhxc_sGsG[:npw, npw:] = fhxc_sGsG[npw:, :npw]
                fhxc_sGsG[npw:, npw:] = fhxc_sGsG[:npw, :npw]

            world.sum(fhxc_sGsG)

            if rank == 0:
                w = Writer('fhxc_%s_%s_%s_%s_%s.gpw' % (self.xc,
                                                        self.method,
                                                        self.tag,
                                                        self.ecut,
                                                        iq))
                w.dimension('sG', ns*npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG / vol)
                w.close()
            world.barrier()

        # Gamma term
        Gvec_Gc = Gvec_qGc[iq0]
        npw = npw_q[iq0]
        l_pw_size = -(-npw // world.size)
        l_pw_range = range(world.rank * l_pw_size,
                           min((world.rank+1) * l_pw_size, npw))
        for d in directions:
            q = np.array([0., 0., 0.])
            q[d] = 1.e-5
            qG_Gv = np.dot(q + Gvec_Gc, bcell_cv)
            fhxc_sGsG = np.zeros((ns*npw, ns*npw), dtype=complex)

            for iG in l_pw_range:
                for jG in range(npw):
                    dGq_c = (Gvec_Gc[iG] + Gvec_Gc[jG])/ 2. + q
                    if (np.abs(dGq_c) < 1.e-12).all():
                        dGq_c = np.array([0.0, 0.0, 0.0001])
                    dGq_v = np.dot(dGq_c, bcell_cv)
                    fx = fx_g.copy()
                    #fx[np.where(2*kf_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    fx[np.where(qc_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    v_c = 4 * np.pi * np.ones(np.shape(fx_g), float)
                    v_c /= np.dot(dGq_v, dGq_v)
                    #v_c[np.where(2*kf_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0
                    v_c[np.where(qc_g < np.dot(dGq_v, dGq_v)**0.5)] = 0.0

                    dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                    dG_v = np.dot(dG_c, bcell_cv)
                    dGr_g = gemmdot(dG_v, r_vg, beta=0.0)

                    fhxc_sGsG[iG, jG] = gd.integrate(np.exp(-1j*dGr_g) * (fx + v_c))
                    if ns == 2:
                        fhxc_sGsG[npw + iG, jG] = gd.integrate(np.exp(-1j*dGr_g) * v_c)
            if ns == 2:
                fhxc_sGsG[:npw, npw:] = fhxc_sGsG[npw:, :npw]
                fhxc_sGsG[npw:, npw:] = fhxc_sGsG[:npw, :npw]

            world.sum(fhxc_sGsG)

            if rank == 0:
                w = Writer('fhxc_%s_%s_%s_%s_0%s.gpw' % (self.xc,
                                                         self.method,
                                                         self.tag,
                                                         self.ecut,
                                                         d))
                w.dimension('sG', ns*npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG / vol)
                w.close()
            world.barrier()

        print(file=self.txt)

    def calculate_local_kernel(self):
        # Standard ALDA exchange kernel
        # Use with care. Results are very difficult to converge
        # Sensitive to densitycut
        ns = self.calc.wfs.nspins
        gd = self.gd

        n_g = self.n_g

        A_x = -(3/4.) * (3/np.pi)**(1/3.)
        fxc_sg = ns * (4 / 9.) * A_x * (ns*n_g)**(-2/3.)
        fxc_sg[np.where(n_g < self.density_cut)] = 0.0

        r_vg = gd.get_grid_point_coordinates()

        for iq in range(len(self.q_points)):
            npw = self.npw_q[iq]
            Gvec_Gc = self.Gvec_qGc[iq]
            l_pw_size = -(-npw // world.size)
            l_pw_range = range(world.rank * l_pw_size,
                               min((world.rank+1) * l_pw_size, npw))
            fhxc_sGsG = np.zeros((ns*npw, ns*npw), dtype=complex)
            for s in range(ns):
                for iG in l_pw_range:
                    for jG in range(npw):
                        fxc = fxc_sg[s].copy()
                        dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                        dG_v = np.dot(dG_c, self.bcell_cv)
                        dGr_g = gemmdot(dG_v, r_vg, beta=0.0)
                        fhxc_sGsG[s*npw+iG, s*npw+jG] = \
                                            gd.integrate(np.exp(-1j*dGr_g)*fxc)

            world.sum(fhxc_sGsG)

            Kc_GG = np.zeros((npw, npw), dtype=complex)
            for iG in range(npw):
                qG_Gv = np.dot(Gvec_Gc[iG] + self.q_points[iq], self.bcell_cv)
                Kc_GG[iG, iG] = 4 * np.pi / np.dot(qG_Gv, qG_Gv)
            fhxc_sGsG += np.tile(Kc_GG, (ns, ns)) * self.vol
            
            if self.paw_correction in [0, 2]:
                print('    Calculating PAW correction', file=self.txt)
                f_paw_sGG = self.add_paw_correction(npw,
                                                    Gvec_Gc,
                                                    bcell_cv,
                                                    self.calc.wfs.setups,
                                                    self.calc.density.D_asp,
                                                    self.atoms.positions/Bohr)
                for s in range(ns):
                    fhxc_sGsG[s*npw:(s+1)*npw, s*npw:(s+1)*npw] += f_paw_sGG[s]

            if rank == 0:
                w = Writer('fhxc_%s_%s_%s_%s.gpw' % (self.xc,
                                                     self.tag,
                                                     self.ecut,
                                                     iq))
                w.dimension('sG', ns*npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG / self.vol)
                w.close()
            world.barrier()

        print(file=self.txt)


    def add_paw_correction(self, npw, Gvec_Gc, bcell_cv, setups, D_asp, R_av):
        # By default only used for ALDA
        ns = self.nspins
        A_x = -(3/4.) * (3/np.pi)**(1/3.)
        KxcPAW_sGG = np.zeros((ns, npw, npw), complex)
        dG_GGv = np.zeros((npw, npw, 3))
        for iG in range(npw):
            for jG in range(npw):
                dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                dG_GGv[iG, jG] =  np.dot(dG_c, bcell_cv)

        for a, setup in enumerate(setups):
            rgd = setup.xc_correction.rgd
            ng = len(rgd.r_g)
            myng = -(-ng // world.size)
            n_qg = setup.xc_correction.n_qg
            nt_qg = setup.xc_correction.nt_qg
            nc_g = setup.xc_correction.nc_g
            nct_g = setup.xc_correction.nct_g
            Y_nL = setup.xc_correction.Y_nL
            dv_g = rgd.dv_g

            D_sp = D_asp[a]
            B_pqL = setup.xc_correction.B_pqL
            D_sLq = np.inner(D_sp, B_pqL.T)

            f_sg = rgd.empty(ns)
            ft_sg = rgd.empty(ns)

            n_sLg = np.dot(D_sLq, n_qg)
            nt_sLg = np.dot(D_sLq, nt_qg)

            # Add core density
            n_sLg[:, 0] += (4 * np.pi)**0.5 / ns * nc_g
            nt_sLg[:, 0] += (4 * np.pi)**0.5 / ns * nct_g

            coefatoms_GG = np.exp(-1j * np.inner(dG_GGv, R_av[a]))
            w = weight_n

            if self.paw_correction == 2:
                Y_nL = [Y_nL[0]]
                w = [1.]
                
            for n, Y_L in enumerate(Y_nL):
                f_sg[:] = 0.0
                n_sg = np.dot(Y_L, n_sLg)
                f_sg = ns * (4 / 9.) * A_x * (ns*n_sg)**(-2/3.)
                if self.density_cut is not None:
                    f_sg[np.where(ns*n_sg < self.density_cut)] = 0.0

                ft_sg[:] = 0.0
                nt_sg = np.dot(Y_L, nt_sLg)
                ft_sg = ns * (4 / 9.) * A_x * (ns*nt_sg)**(-2/3.)
                if self.density_cut is not None:
                    ft_sg[np.where(ns*nt_sg < self.density_cut)] = 0.0

                for i in range(world.rank * myng,
                               min((world.rank + 1) * myng, ng)):
                    coef_GG = np.exp(-1j * np.inner(dG_GGv, R_nv[n])
                                     * rgd.r_g[i])
                    for s in range(len(f_sg)):
                        KxcPAW_sGG[s] += w[n] * np.dot(coef_GG,
                                            (f_sg[s,i]-ft_sg[s,i]) * dv_g[i]) \
                                                       * coefatoms_GG

        world.sum(KxcPAW_sGG)

        return KxcPAW_sGG
                

    def get_fxc_g(self, n_g, index=None):

        if self.xc[2:] == 'LDA':
            return self.get_lda_g(n_g)
        else:
            return self.get_pbe_g(n_g, index=index)
    

    def get_lda_g(self, n_g):

        return (4. / 9.) * self.A_x * n_g**(-2./3.)


    def get_pbe_g(self, n_g, index=None):

        if index is None:
            gradn_vg = self.gradn_vg
        else:
            gradn_vg = self.calc.density.gd.empty(3)
            for v in range(3):
                gradn_vg[v] = (self.gradn_vg[v] +
                               self.gradn_vg[v].flatten()[index]) / 2

        kf_g = (3. * np.pi**2 * n_g)**(1./3.)
        s2_g = np.zeros_like(n_g)
        for v in range(3):
            axpy(1.0, gradn_vg[v]**2, s2_g)
        s2_g /= 4 * kf_g**2 * n_g**2

        e_g = self.A_x * n_g**(4./3.)
        v_g = (4. / 3.) * e_g / n_g
        f_g = (1. / 3.) * v_g / n_g

        kappa = 0.804
        mu = 0.2195149727645171

        denom_g = (1 + mu * s2_g / kappa)
        F_g = 1. + kappa - kappa / denom_g
        Fn_g = -mu / denom_g**2 * 8 * s2_g / (3 * n_g)
        Fnn_g = -11 * Fn_g / (3 * n_g) - 2 * Fn_g**2 / kappa
            
        fxc_g = f_g * F_g
        fxc_g += 2 * v_g * Fn_g
        fxc_g += e_g * Fnn_g

        #Fgrad_vg = np.zeros_like(gradn_vg)
        #Fngrad_vg = np.zeros_like(gradn_vg)
        #for v in range(3):
        #    axpy(1.0, mu / den_g**2 * gradn_vg[v] / (2 * kf_g**2 * n_g**2),
        #         Fgrad_vg[v])
        #    axpy(-8.0, Fgrad_vg[v] / (3 * n_g), Fngrad_vg[v])
        #    axpy(-2.0, Fgrad_vg[v] * Fn_g / kappa, Fngrad_vg[v])


        #tmp = np.zeros_like(fxc_g)
        #tmp1 = np.zeros_like(fxc_g)

        #for v in range(3):
            #self.grad_v[v](Fgrad_vg[v], tmp)
            #axpy(-2.0, tmp * v_g, fxc_g)
            #for u in range(3):
                #self.grad_v[u](Fgrad_vg[u] * tmp, tmp1)
                #axpy(-4.0/kappa, tmp1 * e_g, fxc_g)
            #self.grad_v[v](Fngrad_vg[v], tmp)
            #axpy(-2.0, tmp * e_g, fxc_g)
        #self.laplace(mu / den_g**2 / (2 * kf_g**2 * n_g**2), tmp)
        #axpy(1.0, tmp * e_g, fxc_g)
        
        return fxc_g


    def get_fxc_libxc_g(self, n_g):
        
        gd = self.calc.density.gd.refine()

        xc = XC('GGA_X_' + self.xc[2:])
        #xc = XC('LDA_X')
        #sigma = np.zeros_like(n_g).flat[:]
        xc.set_grid_descriptor(gd)
        sigma_xg, gradn_svg = xc.calculate_sigma(np.array([n_g]))

        dedsigma_xg = np.zeros_like(sigma_xg)
        e_g = np.zeros_like(n_g)
        v_sg = np.array([np.zeros_like(n_g)])
        
        xc.calculate_gga(e_g, np.array([n_g]), v_sg, sigma_xg, dedsigma_xg)

        sigma = sigma_xg[0].flat[:]
        gradn_vg = gradn_svg[0]
        dedsigma_g = dedsigma_xg[0]
        
        libxc = LibXC('GGA_X_' + self.xc[2:])
        #libxc = LibXC('LDA_X')
        libxc.initialize(1)
        libxc_fxc = libxc.xc.calculate_fxc_spinpaired
        
        fxc_g = np.zeros_like(n_g).flat[:]
        d2edndsigma_g = np.zeros_like(n_g).flat[:]
        d2ed2sigma_g = np.zeros_like(n_g).flat[:]

        libxc_fxc(n_g.flat[:], fxc_g, sigma, d2edndsigma_g, d2ed2sigma_g)
        fxc_g = fxc_g.reshape(np.shape(n_g))
        d2edndsigma_g = d2edndsigma_g.reshape(np.shape(n_g))
        d2ed2sigma_g = d2ed2sigma_g.reshape(np.shape(n_g))

        tmp = np.zeros_like(fxc_g)
        tmp1 = np.zeros_like(fxc_g)
        
        #for v in range(3):
            #self.grad_v[v](d2edndsigma_g * gradn_vg[v], tmp)
            #axpy(-4.0, tmp, fxc_g)
            
        #for u in range(3):
            #for v in range(3):
                #self.grad_v[v](d2ed2sigma_g * gradn_vg[u] * gradn_vg[v], tmp)
                #self.grad_v[u](tmp, tmp1)
                #axpy(4.0, tmp1, fxc_g)

        #self.laplace(dedsigma_g, tmp)
        #axpy(2.0, tmp, fxc_g)
            
        return fxc_g[::2,::2,::2]


    def get_numerical_fxc_sg(self, n_sg):

        gd = self.calc.density.gd.refine()

        delta = 1.e-4
        if self.xc[2:] == 'LDA':
            xc = XC('LDA_X')
            v1xc_sg = np.zeros_like(n_sg)
            v2xc_sg = np.zeros_like(n_sg)
            xc.calculate(gd, (1+delta)*n_sg, v1xc_sg)
            xc.calculate(gd, (1-delta)*n_sg, v2xc_sg)
            fxc_sg = (v1xc_sg - v2xc_sg) / (2 * delta * n_sg)
        else:
            fxc_sg = np.zeros_like(n_sg)
            xc = XC('GGA_X_' + self.xc[2:])
            vxc_sg = np.zeros_like(n_sg)
            xc.calculate(gd, n_sg, vxc_sg)
            for s in range(len(n_sg)):
                for x in range(len(n_sg[0])):
                    for y in range(len(n_sg[0,0])):
                        for z in range(len(n_sg[0,0,0])):
                            v1xc_sg = np.zeros_like(n_sg)
                            n1_sg = n_sg.copy()
                            n1_sg[s,x,y,z] *= (1 + delta)
                            xc.calculate(gd, n1_sg, v1xc_sg)
                            fxc_sg[s,x,y,z] = \
                (v1xc_sg[s,x,y,z] - vxc_sg[s,x,y,z]) / (delta * n_sg[s,x,y,z])
                            
        return fxc_sg[:,::2,::2,::2]
