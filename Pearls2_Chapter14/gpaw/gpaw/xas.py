from __future__ import print_function
import pickle
from math import log, pi, sqrt

import numpy as np
from ase.units import Hartree
from gpaw.utilities.cg import CG
import gpaw.mpi as mpi


class XAS:
    def __init__(self, paw, mode="xas", center=None, spin=0):
        wfs = paw.wfs
        self.orthogonal = wfs.gd.orthogonal
        self.cell_cv = np.array(wfs.gd.cell_cv)
        assert wfs.world.size == 1 #assert not mpi.parallel
        #assert wfs.gd.orthogonal
        
        #
        # to allow spin polarized calclulation
        #
        
        nkpts = len(wfs.kd.ibzk_kc)

        # the following lines are to stop the user to make mistakes
        #if mode is not "xes" and spin == 1:
        #    raise RuntimeError(
        #        "The core hole is always in spin 0: please use spin=0")

        if wfs.nspins == 1:
            if spin is not 0:
                raise RuntimeError(
                    "use spin=0 for a spin paired calculation")            
            nocc = wfs.setups.nvalence // 2
            self.list_kpts = range(nkpts)
        else:
            self.list_kpts=[]

            if spin is not 0 and spin is not 1:
                print("spin",spin)    
                raise RuntimeError(
                    "use either spin=0 or spin=1")            

            #find kpoints with correct spin
            for i, kpt in  enumerate(wfs.kpt_u):
                if kpt.s == spin:
                    self.list_kpts.append(i)
                print(self.list_kpts)
            assert len(self.list_kpts) == nkpts
                        
            #find number of occupied orbitals, if no fermi smearing
            nocc = 0.
            for i in self.list_kpts:
                nocc += sum(wfs.kpt_u[i].f_n)
            nocc = int(nocc + 0.5)
            print("nocc", nocc)
                 

        # look for the center with the corehole
        if center is not None:
            #print "center", center
            setup = wfs.setups[center]
            a = center
            
        else:
            for a, setup in enumerate(wfs.setups):
                if setup.phicorehole_g is not None:  
                    break

        A_ci = setup.A_ci

        # xas, xes or all modes
        if mode == "xas":
            n_start = nocc
            n_end = wfs.bd.nbands  
            n =  wfs.bd.nbands - nocc
        elif mode == "xes":
            n_start = 0
            n_end = nocc  
            n = nocc
        elif mode == "all":
            n_start = 0
            n_end = wfs.bd.nbands 
            n = wfs.bd.nbands
        else:
            raise RuntimeError(
                "wrong keyword for 'mode', use 'xas', 'xes' or 'all'")

        self.n = n
            
        self.eps_n = np.empty(nkpts * n)
        self.sigma_cn = np.empty((3, nkpts * n), complex)
        n1 = 0
        for kpt in wfs.kpt_u:
            if kpt.s != spin:
                continue
            
            n2 = n1 + n
            self.eps_n[n1:n2] = kpt.eps_n[n_start:n_end] * Hartree
            P_ni = kpt.P_ani[a][n_start:n_end]
            a_cn = np.inner(A_ci, P_ni)
            weight = kpt.weight * wfs.nspins / 2
            print("weight", weight)
            print(a_cn.shape, self.sigma_cn.shape)
            self.sigma_cn[:, n1:n2] =  weight **0.5 * a_cn #.real
            n1 = n2

        self.symmetry = wfs.kd.symmetry

    def get_spectra(self, fwhm=0.5, E_in=None, linbroad=None, N=1000, kpoint=None,
                    proj=None,  proj_xyz=True, stick=False):
        """Calculate spectra.

        Parameters:
        
        fwhm:
          the full width half maximum in eV for gaussian broadening
        linbroad:
          a list of three numbers, the first fwhm2, the second the value
          where the linear increase starts and the third the value where
          the broadening has reached fwhm2. example [0.5, 540, 550]
        E_in:
          a list of energy values where the spectrum is to be computed
          if None the orbital energies will be used to compute the energy 
          range
        N:
          the number of bins in the broadened spectrum. If E_in is given N
          has no effect
        kpoint:
          select a specific k-point to calculate spectrum for
        proj:
          a list of vectors to project the transition dipole on. Default
          is None then only x,y,z components are calculated.  a_stick and
          a_c squares of the transition moments in resp. direction
        proj_xyz:
          if True keep projections in x, y and z. a_stck and a_c will have
          length 3 + len(proj). if False only those projections
          defined by proj keyword, a_stick and a_c will have length len(proj)
        stick:
          if False return broadened spectrum, if True return stick spectrum
          
        Symmtrization has been moved inside get_spectra because we want to
        symmtrice squares of transition dipoles."""
        
        # eps_n = self.eps_n[k_in*self.n: (k_in+1)*self.n -1]
         
        # proj keyword, check normalization of incoming vectors
        if proj_xyz:
            proj_3 = np.array([[1,0,0],[0,1,0],[0,0,1]], float)
        else:
            proj_3 = np.array([],float)

        if proj is not None:
            assert self.orthogonal # does not work for nonorthogonal unit cells
            proj_2 = np.array(proj,float)
            if len(proj_2.shape) == 1:
                proj_2 = np.array([proj],float)

            for i,p in enumerate(proj_2):
                if sum(p ** 2) ** 0.5 != 1.0:
                    print("proj_2 %s not normalized" %i)
                    proj_2[i] /=  sum(p ** 2) ** 0.5
            
            proj_tmp = np.zeros((proj_3.shape[0] + proj_2.shape[0], 3), float)
                
            for i, p in enumerate(proj_3):
                proj_tmp[i,:] = proj_3[i,:]

            for i, p in enumerate(proj_2):
                proj_tmp[proj_3.shape[0] + i,:] = proj_2[i,:]
        
            proj_3 = proj_tmp.copy()
        
        # now symmetrize
        sigma2_cn = np.zeros((proj_3.shape[0], self.sigma_cn.shape[1]),float) 
        
        if self.symmetry is not None:
            for i,p in enumerate(proj_3):
                for op_cc in self.symmetry.op_scc:
                    op_vv = np.dot(np.linalg.inv(self.cell_cv),
                                   np.dot(op_cc, self.cell_cv))
                    s_tmp = np.dot(p, np.dot(op_vv, self.sigma_cn))
                    sigma2_cn[i,:] += (s_tmp * np.conjugate(s_tmp) ).real 
            sigma2_cn /= len(self.symmetry.op_scc)
            
        else:
            for i,p in enumerate(proj_3):
                s_tmp = np.dot(p, self.sigma_cn)
                sigma2_cn[i,:] += (s_tmp * np.conjugate(s_tmp) ).real 
                
        eps_n = self.eps_n[:]

        if kpoint is not None:
            eps_start = kpoint*self.n
            eps_end = (kpoint+1)*self.n
        else: 
            eps_start = 0
            eps_end = len(self.eps_n)
            
       
        # return stick spectrum if stick=True
        if stick:
            e_stick = eps_n[eps_start:eps_end]
            a_stick = sigma2_cn[:,eps_start:eps_end]

            return e_stick, a_stick

        # else return broadened spectrum
        else:
            if E_in is not None:
                e = np.array(E_in)
            else:
                emin = min(eps_n) - 2 * fwhm
                emax = max(eps_n) + 2 * fwhm
                e = emin + np.arange(N + 1) * ((emax - emin) / N)
            
            a_c = np.zeros((len(sigma2_cn), len(e)))
            
            if linbroad is None:
                #constant broadening fwhm
                alpha = 4 * log(2) / fwhm**2
            
                for n, eps in enumerate(eps_n[eps_start:eps_end]):
                    x = -alpha * (e - eps)**2
                    x = np.clip(x, -100.0, 100.0)
                    a_c += np.outer(sigma2_cn[:, n + eps_start],
                                        (alpha / pi)**0.5 * np.exp(x))
            else:

                # constant broadening fwhm until linbroad[1] and a
                # constant broadening over linbroad[2] with fwhm2=
                # linbroad[0]
                fwhm2 = linbroad[0]
                lin_e1 = linbroad[1]
                lin_e2 = linbroad[2]
                print("fwhm", fwhm, fwhm2, lin_e1, lin_e2)
                for n, eps in enumerate(eps_n):
                    if eps < lin_e1:
                        alpha = 4*log(2) / fwhm**2
                    elif eps <=  lin_e2:
                        fwhm_lin = (fwhm + (eps - lin_e1) *
                                    (fwhm2 - fwhm) / (lin_e2 - lin_e1))
                        alpha = 4*log(2) / fwhm_lin**2
                    elif eps >= lin_e2:
                        alpha =  4*log(2) / fwhm2**2
                        
                    x = -alpha * (e - eps)**2
                    x = np.clip(x, -100.0, 100.0)
                    a_c += np.outer(sigma2_cn[:, n],
                                     (alpha / pi)**0.5 * np.exp(x))
                
            return  e, a_c


class RecursionMethod:
    """This class implements the Haydock recursion method. """

    def __init__(self, paw=None, filename=None,
                 tol=1e-10, maxiter=100, proj=None,
                 proj_xyz=True):

        if paw is not None:
            wfs = paw.wfs
            assert wfs.gd.orthogonal

            self.wfs = wfs
            self.hamiltonian = paw.hamiltonian
            self.nkpts = len(wfs.kd.ibzk_kc) * wfs.nspins
            self.nmykpts = len(wfs.kpt_u)

            self.k1 = wfs.kd.comm.rank * self.nmykpts
            self.k2 = self.k1 + self.nmykpts
            
            print("k1", self.k1, "k2",self.k2)

            # put spin and weight index in the columns corresponding
            # to this processors k-points
            self.spin_k = np.zeros(self.nkpts,int)
            self.weight_k = np.zeros(self.nkpts)

            for n, i in enumerate(range(self.k1, self.k2)):
                self.spin_k[i] = wfs.kpt_u[n].s
                self.weight_k[i] = wfs.kpt_u[n].weight
                     
            self.op_scc = None
            if wfs.kd.symmetry is not None:
                self.op_scc = wfs.kd.symmetry.op_scc
        else:
            self.k1 = 0
            self.k2 = None
            self.wfs = None
            wfs = None
            
        self.tol = tol
        self.maxiter = maxiter

     
        if filename is not None:
            self.read(filename)
            if wfs is not None:
                self.allocate_tmp_arrays()
        else:
            self.initialize_start_vector(proj=proj,proj_xyz=proj_xyz)

    def read(self, filename):
        data = pickle.load(open(filename))
        self.nkpts = data['nkpts']
        if 'swaps' in data:
            # This is an old file:
            self.op_scc = np.array([np.identity(3, int)[list(swap)]
                                    for swap in data['swaps']])
        else:
            self.op_scc = data['symmetry operations']
        self.weight_k = data['weight_k']
        self.spin_k = data['spin_k']
        self.dim = data['dim']
        k1, k2 = self.k1, self.k2
        if k2 is None:
            k2 = self.nkpts
        a_kci, b_kci = data['ab']
        self.a_uci = a_kci[k1:k2].copy()
        self.b_uci = b_kci[k1:k2].copy()

        if self.wfs is not None and 'arrays' in data:
            print('reading arrays')
            w_kcG, wold_kcG, y_kcG = data['arrays']
            i = [slice(k1, k2), slice(0, self.dim)] + self.wfs.gd.get_slice()
            self.w_ucG = w_kcG[i].copy()
            self.wold_ucG = wold_kcG[i].copy()
            self.y_ucG = y_kcG[i].copy()
            
    def write(self, filename, mode=''):
        assert self.wfs is not None
        kpt_comm = self.wfs.kd.comm
        gd = self.wfs.gd

        if gd.comm.rank == 0:
            if kpt_comm.rank == 0:
                nmyu, dim, ni = self.a_uci.shape
                a_kci = np.empty((kpt_comm.size, nmyu, dim, ni),
                                  self.wfs.dtype)
                b_kci = np.empty((kpt_comm.size, nmyu, dim, ni),
                                  self.wfs.dtype)

                kpt_comm.gather(self.a_uci, 0, a_kci)
                kpt_comm.gather(self.b_uci, 0, b_kci)
                kpt_comm.sum(self.spin_k, 0)
                kpt_comm.sum(self.weight_k, 0)
                
                a_kci.shape = (self.nkpts, dim, ni)
                b_kci.shape = (self.nkpts, dim, ni)
                data = {'ab': (a_kci, b_kci),
                        'nkpts': self.nkpts,
                        'symmetry operations': self.op_scc,
                        'weight_k':self.weight_k,'spin_k':self.spin_k,
                        'dim':dim}
            else:
                kpt_comm.gather(self.a_uci, 0)
                kpt_comm.gather(self.b_uci, 0)
                kpt_comm.sum(self.spin_k, 0)
                kpt_comm.sum(self.weight_k, 0)
            
        if mode == 'all':
            w0_ucG = gd.collect(self.w_ucG)
            wold0_ucG = gd.collect(self.wold_ucG)
            y0_ucG = gd.collect(self.y_ucG)
            if gd.comm.rank == 0:
                if kpt_comm.rank == 0:
                    w_kcG = gd.empty((self.nkpts, dim), self.wfs.dtype,
                                     global_array=True)
                    wold_kcG = gd.empty((self.nkpts, dim), self.wfs.dtype,
                                        global_array=True)
                    y_kcG = gd.empty((self.nkpts, dim), self.wfs.dtype,
                                     global_array=True)
                    kpt_comm.gather(w0_ucG, 0, w_kcG)
                    kpt_comm.gather(wold0_ucG, 0, wold_kcG)
                    kpt_comm.gather(y0_ucG, 0, y_kcG)
                    data['arrays'] = (w_kcG, wold_kcG, y_kcG)
                else:
                    kpt_comm.gather(w0_ucG, 0)
                    kpt_comm.gather(wold0_ucG, 0)
                    kpt_comm.gather(y0_ucG, 0)

        if self.wfs.world.rank == 0:
            pickle.dump(data, open(filename, 'w'))

    def allocate_tmp_arrays(self):
        
        self.tmp1_cG = self.wfs.gd.zeros(self.dim, self.wfs.dtype)
        self.tmp2_cG = self.wfs.gd.zeros(self.dim, self.wfs.dtype)
        self.z_cG = self.wfs.gd.zeros(self.dim, self.wfs.dtype)
        
    def initialize_start_vector(self, proj=None, proj_xyz=True):
        # proj is one list of vectors [[e1_x,e1_y,e1_z],[e2_x,e2_y,e2_z]]
        #( or [ex,ey,ez] if only one projection )
        # that the spectrum will be projected on 
        # default is to only calculate the averaged spectrum
        # if proj_xyz is True, keep projection in x,y,z, if False
        # only calculate the projections in proj
        
        # Create initial wave function:
        nmykpts = self.nmykpts
        
        for a, setup in enumerate(self.wfs.setups):
            if setup.phicorehole_g is not None:
                break
        A_ci = setup.A_ci

        #
        # proj keyword
        #

        #check normalization of incoming vectors
        if proj is not None:
            proj_2 = np.array(proj,float)
            if len(proj_2.shape) == 1:
                proj_2 = np.array([proj],float)
            
            for i,p in enumerate(proj_2):
                if sum(p ** 2) ** 0.5 != 1.0:
                    print("proj_2 %s not normalized" %i)
                    proj_2[i] /=  sum(p ** 2) ** 0.5

            proj_tmp = []
            for p in proj_2:
               proj_tmp.append(np.dot(p, A_ci))
            proj_tmp = np.array(proj_tmp, float)   

            # if proj_xyz is True, append projections to A_ci
            if proj_xyz:
                A_ci_tmp = np.zeros((3 + proj_2.shape[0], A_ci.shape[1]))
                A_ci_tmp[0:3,:] = A_ci 
                A_ci_tmp[3:,:]= proj_tmp

            # otherwise, replace A_ci by projections
            else:
                A_ci_tmp = np.zeros((proj_2.shape[0], A_ci.shape[1]))
                A_ci_tmp = proj_tmp
            A_ci = A_ci_tmp

        self.dim = len(A_ci)

        self.allocate_tmp_arrays()

        self.w_ucG = self.wfs.gd.zeros((nmykpts, self.dim), self.wfs.dtype)
        self.wold_ucG = self.wfs.gd.zeros((nmykpts, self.dim), self.wfs.dtype)
        self.y_ucG = self.wfs.gd.zeros((nmykpts, self.dim), self.wfs.dtype)
            
        self.a_uci = np.zeros((nmykpts, self.dim, 0), self.wfs.dtype)
        self.b_uci = np.zeros((nmykpts, self.dim, 0), self.wfs.dtype)

        A_aci = self.wfs.pt.dict(3, zero=True)
        if a in A_aci:
            A_aci[a] = A_ci.astype(self.wfs.dtype)
        for u in range(nmykpts):
            self.wfs.pt.add(self.w_ucG[u], A_aci, u)

    def run(self, nsteps, inverse_overlap="exact"):

        if inverse_overlap == "exact":
            self.solver = self.solve
        elif inverse_overlap == "approximate":
            self.solver = self.solve2
        elif inverse_overlap == "noinverse":
            self.solver = self.solve3
        else:
            raise RuntimeError("""Error, inverse_solver must be either 'exact',
            'approximate' or 'noinverse' """)
            

        ni = self.a_uci.shape[2]
        a_uci = np.empty((self.nmykpts, self.dim, ni + nsteps), self.wfs.dtype)
        b_uci = np.empty((self.nmykpts, self.dim, ni + nsteps), self.wfs.dtype)
        a_uci[:, :, :ni]  = self.a_uci
        b_uci[:, :, :ni]  = self.b_uci
        self.a_uci = a_uci
        self.b_uci = b_uci

        for u in range(self.nmykpts):
            for i in range(nsteps):
                self.step(u, ni + i)

    def step(self, u, i):
        print(u, i)
        integrate = self.wfs.gd.integrate
        w_cG = self.w_ucG[u]
        y_cG = self.y_ucG[u]
        wold_cG = self.wold_ucG[u]
        z_cG = self.z_cG
        
        self.solver(w_cG, self.z_cG, u)
        I_c = np.reshape(integrate(np.conjugate(z_cG) * w_cG)**-0.5,
                          (self.dim, 1, 1, 1))
        z_cG *= I_c
        w_cG *= I_c
        
        if i != 0:
            b_c =  1.0 / I_c 
        else:
            b_c = np.reshape(np.zeros(self.dim), (self.dim, 1, 1, 1))
    
        self.hamiltonian.apply(z_cG, y_cG, self.wfs, self.wfs.kpt_u[u])
        a_c = np.reshape(integrate(np.conjugate(z_cG) * y_cG), (self.dim, 1, 1, 1))
        wnew_cG = (y_cG - a_c * w_cG - b_c * wold_cG)
        wold_cG[:] = w_cG
        w_cG[:] = wnew_cG
        self.a_uci[u, :, i] = a_c[:, 0, 0, 0]
        self.b_uci[u, :, i] = b_c[:, 0, 0, 0]


    def continued_fraction(self, e, k, c, i, imax):
        a_i = self.a_uci[k, c]
        b_i = self.b_uci[k, c]
        if i == imax - 2:
            return self.terminator(a_i[i], b_i[i], e)
        return 1.0 / (a_i[i] - e -
                      b_i[i + 1]**2 *
                      self.continued_fraction(e, k, c, i + 1, imax))

    def get_spectra(self, eps_s, delta=0.1, imax=None, kpoint=None, fwhm=None,
                    linbroad=None, spin=0):
        assert not mpi.parallel

        # the following lines are to stop the user to make mistakes
        #if spin == 1:
        #    raise RuntimeError(
        #        "The core hole is always in spin 0: please use spin=0")

        n = len(eps_s)
                
        sigma_cn = np.zeros((self.dim, n))
        if imax is None:
            imax = self.a_uci.shape[2]
        eps_n = (eps_s + delta * 1.0j) / Hartree
                
        # if a certain k-point is chosen
        if kpoint is not None:
             for c in range(self.dim):
                sigma_cn[c] += self.continued_fraction(eps_n, kpoint, c,
                                                       0, imax).imag
        else:
            for k in range(self.nkpts):
                print('kpoint', k, 'spin_k', self.spin_k[k], spin, 'weight', self.weight_k[k]) 
                if self.spin_k[k] == spin:
                    weight = self.weight_k[k]
                    for c in range(self.dim):
                        sigma_cn[c] += weight*self.continued_fraction(eps_n, k, c,
                                                                      0, imax).imag
                        
        if self.op_scc is not None:
            sigma0_cn = sigma_cn
            sigma_cn = np.zeros((self.dim, n))
            for op_cc in self.op_scc:
                sigma_cn += np.dot(op_cc**2, sigma0_cn)
            sigma_cn /= len(self.op_scc)

        # gaussian broadening 
        if fwhm is not None:
            sigma_tmp = np.zeros(sigma_cn.shape)

            #constant broadening fwhm
            if linbroad is None:
                alpha = 4 * log(2) / fwhm**2
                for n, eps in enumerate(eps_s):
                    x = -alpha * (eps_s - eps)**2
                    x = np.clip(x, -100.0, 100.0)
                    sigma_tmp += np.outer(sigma_cn[:,n],
                                        (alpha / pi)**0.5 * np.exp(x))

            else:
                # constant broadening fwhm until linbroad[1] and a
                # constant broadening over linbroad[2] with fwhm2=
                # linbroad[0]
                fwhm2 = linbroad[0]
                lin_e1 = linbroad[1]
                lin_e2 = linbroad[2]
                for n, eps in enumerate(eps_s):
                    if eps < lin_e1:
                        alpha = 4*log(2) / fwhm**2
                    elif eps <=  lin_e2:
                        fwhm_lin = (fwhm + (eps - lin_e1) *
                                (fwhm2 - fwhm) / (lin_e2 - lin_e1))
                        alpha = 4*log(2) / fwhm_lin**2
                    elif eps >= lin_e2:
                        alpha =  4*log(2) / fwhm2**2

                    x = -alpha * (eps_s - eps)**2
                    x = np.clip(x, -100.0, 100.0)
                    sigma_tmp += np.outer(sigma_cn[:, n],
                                        (alpha / pi)**0.5 * np.exp(x))
            sigma_cn = sigma_tmp
                    

        return sigma_cn
    
    def solve(self, w_cG, z_cG, u):
        # exact inverse overlap
        self.wfs.overlap.apply_inverse(w_cG, self.tmp1_cG, self.wfs,
                                       self.wfs.kpt_u[u])
        self.u = u
        CG(self, z_cG, self.tmp1_cG,
           tolerance=self.tol, maxiter=self.maxiter)

    def solve2(self, w_cG, z_cG, u):
        # approximate inverse overlap
        self.wfs.overlap.apply_inverse(w_cG, z_cG, self.wfs, self.wfs.kpt_u[u])

        self.u = u

    def solve3(self, w_cG, z_cG, u):
        # no inverse overlap
        z_cG[:] =  w_cG
        self.u = u

    

    def sum(self, a):
        self.wfs.gd.comm.sum(a)
        return a
    
    def __call__(self, in_cG, out_cG):
        """Function that is called by CG. It returns S~-1Sx_in in x_out
        """

        kpt = self.wfs.kpt_u[self.u]
        self.wfs.overlap.apply(in_cG, self.tmp2_cG, self.wfs, kpt)
        self.wfs.overlap.apply_inverse(self.tmp2_cG, out_cG, self.wfs, kpt)

    def terminator(self, a, b, e):
        """ Analytic formula to terminate the continued fraction from
        [R Haydock, V Heine, and M J Kelly, J Phys. C: Solid State Physics, Vol 8, (1975), 2591-2605]
        """

        return 0.5 * (e - a - ((e - a)**2 - 4 * b**2)**0.5 / b**2)

    def duplicate_coefficients(self, nsteps, ntimes):
        n1 = self.a_uci.shape[0]
        n2 = self.a_uci.shape[1]
        ni = self.a_uci.shape[2]
        type_code = self.a_uci.dtype.name #typecode()
        a_uci = np.empty((n1 , n2, ni + nsteps*ntimes), type_code)
        b_uci = np.empty((n1, n2, ni + nsteps*ntimes), type_code)
        a_uci[:, :, :ni]  = self.a_uci
        b_uci[:, :, :ni]  = self.b_uci
        
        ni1 = ni
        ni2 = ni +nsteps
        for i in range(ntimes):
            a_uci[:, :, ni1: ni2] = a_uci[:, :, ni-nsteps:ni]
            b_uci[:, :, ni1: ni2] = b_uci[:, :, ni-nsteps:ni]
            ni1 += nsteps
            ni2 += nsteps
        self.a_uci = a_uci
        self.b_uci = b_uci


def write_spectrum(a,b, filename):
    f=open(filename, 'w')         
    print(f, a.shape, b.shape)        
    
    for i in range(a.shape[0]):
        print("%g" % a[i], b[0,i] +b[1,i] +b[2,i], end=' ', file=f)
        for b2 in b:
            print("%g" % b2[i], end=' ', file=f)
        print(file=f)
    f.close()      


                                                                                                                 
