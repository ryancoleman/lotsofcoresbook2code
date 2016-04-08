from __future__ import print_function
import numpy as np
from math import sqrt, pi
import pickle
from ase.units import Hartree, Bohr
from gpaw.mpi import rank
from gpaw.response.chi import CHI


class DF(CHI):
    """This class defines dielectric function related physical quantities."""

    def __init__(self,
                 calc=None,
                 nbands=None,
                 w=None,
                 q=None,
                 eshift=None,
                 ecut=10.,
                 density_cut=None,
                 G_plus_q=False,
                 eta=0.2,
                 rpad=None,
                 vcut=None,
                 ftol=1e-7,
                 txt=None,
                 xc='ALDA',
                 print_xc_scf=False,
                 hilbert_trans=True,
                 time_ordered=False,
                 optical_limit=False,
                 comm=None,
                 kcommsize=None):

        CHI.__init__(self, calc=calc, nbands=nbands, w=w, q=q, eshift=eshift,
                     ecut=ecut, density_cut=density_cut,
                     G_plus_q=G_plus_q, eta=eta, rpad=rpad, vcut=vcut,
                     ftol=ftol, txt=txt, xc=xc, hilbert_trans=hilbert_trans,
                     time_ordered=time_ordered, optical_limit=optical_limit,
                     comm=comm, kcommsize=kcommsize)

        self.df_flag = False
        self.print_bootstrap = print_xc_scf
        self.df1_w = None  # NLF RPA
        self.df2_w = None  # LF RPA
        self.df3_w = None  # NLF ALDA
        self.df4_w = None  # LF ALDA

    def get_dielectric_matrix(self,
                              xc='RPA',
                              overwritechi0=False,
                              symmetric=True,
                              chi0_wGG=None,
                              calc=None,
                              vcut=None,
                              dir=None):
        if self.chi0_wGG is None and chi0_wGG is None:
            self.initialize()
            self.calculate()
        elif self.chi0_wGG is None and chi0_wGG is not None:
            #Read from file and reinitialize 
            self.xc = xc

            from gpaw.response.parallel import par_read 
            self.chi0_wGG = par_read(chi0_wGG, 'chi0_wGG')
            self.nvalbands = self.nbands
            #self.parallel_init() # parallelization not yet implemented
            self.Nw_local = self.Nw  # parallelization not yet implemented
            if self.calc is None:
                from gpaw import GPAW
                self.calc = GPAW(calc,txt=None)
            if self.xc == 'ALDA' or self.xc == 'ALDA_X':
                from gpaw.response.kernel import calculate_Kxc
                from gpaw.grid_descriptor import GridDescriptor
                from gpaw.mpi import world, rank, size, serial_comm
                    
                self.pbc = self.calc.atoms.pbc
                self.gd = GridDescriptor(self.calc.wfs.gd.N_c*self.rpad, self.acell_cv,
                                         pbc_c=True, comm=serial_comm)

                R_av = self.calc.atoms.positions / Bohr
                nt_sg = self.calc.density.nt_sG
                    
                if (self.rpad > 1).any() or (self.pbc - True).any():
                    nt_sG = self.gd.zeros(self.nspins)
                    #nt_sG = np.zeros([self.nspins, self.nG[0], self.nG[1], self.nG[2]])
                    for s in range(self.nspins):
                        nt_G = self.pad(nt_sg[s])
                        nt_sG[s] = nt_G
                else:
                    nt_sG = nt_sg
                        
                self.Kxc_sGG = calculate_Kxc(self.gd, 
                                             nt_sG,
                                             self.npw, self.Gvec_Gc,
                                             self.gd.N_c, self.vol,
                                             self.bcell_cv, R_av,
                                             self.calc.wfs.setups,
                                             self.calc.density.D_asp,
                                             functional=self.xc,
                                             density_cut=self.density_cut)

        if overwritechi0:
            dm_wGG = self.chi0_wGG
        else:
            dm_wGG = np.zeros_like(self.chi0_wGG)

        if dir is None:
            q_c = self.q_c
        else:
            q_c = np.diag((1,1,1))[dir] * self.qopt
            self.chi0_wGG[:,0,:] = self.chi00G_wGv[:,:,dir]
            self.chi0_wGG[:,:,0] = self.chi0G0_wGv[:,:,dir]
        
        from gpaw.response.kernel import calculate_Kc, CoulombKernel
        kernel = CoulombKernel(vcut=self.vcut,
                               pbc=self.calc.atoms.pbc,
                               cell=self.acell_cv)
        self.Kc_GG = kernel.calculate_Kc(q_c,
                                         self.Gvec_Gc,
                                         self.bcell_cv,
                                         symmetric=symmetric)
        #self.Kc_GG = calculate_Kc(q_c,
        #                          self.Gvec_Gc,
        #                          self.acell_cv,
        #                          self.bcell_cv,
        #                          self.pbc,
        #                          self.vcut,
        #                          symmetric=symmetric)

        tmp_GG = np.eye(self.npw, self.npw)

        if xc == 'RPA':
            self.printtxt('Use RPA.')
            for iw in range(self.Nw_local):
                dm_wGG[iw] = tmp_GG - self.Kc_GG * self.chi0_wGG[iw]
                
        elif xc == 'ALDA':
            self.printtxt('Use ALDA kernel.')
            # E_LDA = 1 - v_c chi0 (1-fxc chi0)^-1
            # http://prb.aps.org/pdf/PRB/v33/i10/p7017_1 eq. 4
            A_wGG = self.chi0_wGG.copy()
            for iw in range(self.Nw_local):
                A_wGG[iw] = np.dot(self.chi0_wGG[iw], np.linalg.inv(tmp_GG - np.dot(self.Kxc_sGG[0], self.chi0_wGG[iw])))
    
            for iw in range(self.Nw_local):
                dm_wGG[iw] = tmp_GG - self.Kc_GG * A_wGG[iw]                

        return dm_wGG


    def get_inverse_dielectric_matrix(self, xc='RPA'):

        dm_wGG = self.get_dielectric_matrix(xc=xc)
        dminv_wGG = np.zeros_like(dm_wGG)
        for iw in range(self.Nw_local):
            dminv_wGG[iw] = np.linalg.inv(dm_wGG[iw])
        return dminv_wGG


    def get_chi(self, xc='RPA'):
        """Solve Dyson's equation."""

        if self.chi0_wGG is None:
            self.initialize()
            self.calculate()
        else:
            pass # read from file and re-initializing .... need to be implemented

        kernel_GG = np.zeros((self.npw, self.npw), dtype=complex)
        chi_wGG = np.zeros_like(self.chi0_wGG)

        # Coulomb kernel
        for iG in range(self.npw):
            qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
            kernel_GG[iG,iG] = 4 * pi / np.dot(qG, qG)
            
        if xc == 'ALDA':
            kernel_GG += self.Kxc_sGG[0]

        for iw in range(self.Nw_local):
            tmp_GG = np.eye(self.npw, self.npw) - np.dot(self.chi0_wGG[iw], kernel_GG)
            chi_wGG[iw] = np.dot(np.linalg.inv(tmp_GG) , self.chi0_wGG[iw])

        return chi_wGG
    

    def get_dielectric_function(self, xc='RPA', dir=None):
        """Calculate the dielectric function. Returns df1_w and df2_w.

        Parameters:

        df1_w: ndarray
            Dielectric function without local field correction.
        df2_w: ndarray
            Dielectric function with local field correction.
        """

        if not self.optical_limit:
            assert dir is None
            
        if self.df_flag is False:
            dm_wGG = self.get_dielectric_matrix(xc=xc, dir=dir)

            Nw_local = dm_wGG.shape[0]
            dfNLF_w = np.zeros(Nw_local, dtype = complex)
            dfLFC_w = np.zeros(Nw_local, dtype = complex)
            df1_w = np.zeros(self.Nw, dtype = complex)
            df2_w = np.zeros(self.Nw, dtype = complex)

            for iw in range(Nw_local):
                tmp_GG = dm_wGG[iw]
                dfLFC_w[iw] = 1. / np.linalg.inv(tmp_GG)[0, 0]
                dfNLF_w[iw] = tmp_GG[0, 0]

            self.wcomm.all_gather(dfNLF_w, df1_w)
            self.wcomm.all_gather(dfLFC_w, df2_w)

            if xc == 'RPA':
                self.df1_w = df1_w
                self.df2_w = df2_w
            elif xc=='ALDA' or xc=='ALDA_X':
                self.df3_w = df1_w
                self.df4_w = df2_w                

        if xc == 'RPA':
            return self.df1_w, self.df2_w
        elif xc == 'ALDA' or xc=='ALDA_X':
            return self.df3_w, self.df4_w


    def get_surface_response_function(self, z0=0., filename='surf_EELS'):
        """Calculate surface response function."""

        if self.chi0_wGG is None:
            self.initialize()
            self.calculate()


        g_w2 = np.zeros((self.Nw,2), dtype=complex)
        assert self.acell_cv[0, 2] == 0. and self.acell_cv[1, 2] == 0.

        Nz = self.gd.N_c[2] # number of points in z direction
        tmp = np.zeros(Nz, dtype=int)
        nGz = 0         # number of G_z 
        for i in range(self.npw):
            if self.Gvec_Gc[i, 0] == 0 and self.Gvec_Gc[i, 1] == 0:
                tmp[nGz] = self.Gvec_Gc[i, 2]
                nGz += 1
        assert (np.abs(self.Gvec_Gc[:nGz, :2]) < 1e-10).all()

        for id, xc in enumerate(['RPA', 'ALDA']):
            chi_wGG = self.get_chi(xc=xc)
    
            # The first nGz are all Gx=0 and Gy=0 component
            chi_wgg_LFC = chi_wGG[:, :nGz, :nGz]
            del chi_wGG
            chi_wzz_LFC = np.zeros((self.Nw_local, Nz, Nz), dtype=complex)        
    
            # Fourier transform of chi_wgg to chi_wzz
            Gz_g = tmp[:nGz] * self.bcell_cv[2,2]
            z_z = np.linspace(0, self.acell_cv[2,2]-self.gd.h_cv[2,2], Nz)
            phase1_zg = np.exp(1j  * np.outer(z_z, Gz_g))
            phase2_gz = np.exp(-1j * np.outer(Gz_g, z_z))
    
            for iw in range(self.Nw_local):
                chi_wzz_LFC[iw] = np.dot(np.dot(phase1_zg, chi_wgg_LFC[iw]), phase2_gz)
            chi_wzz_LFC /= self.acell_cv[2,2]        
    
            # Get surface response function
    
            z_z -= z0 / Bohr
            q_v = np.dot(self.q_c, self.bcell_cv)
            qq = sqrt(np.inner(q_v, q_v))
            phase1_1z = np.array([np.exp(qq*z_z)])
            phase2_z1 = np.exp(qq*z_z)
    
            tmp_w = np.zeros(self.Nw_local, dtype=complex)        
            for iw in range(self.Nw_local):
                tmp_w[iw] = np.dot(np.dot(phase1_1z, chi_wzz_LFC[iw]), phase2_z1)[0]            
    
            tmp_w *= -2 * pi / qq * self.gd.h_cv[2,2]**2        
            g_w = np.zeros(self.Nw, dtype=complex)
            self.wcomm.all_gather(tmp_w, g_w)
            g_w2[:, id] = g_w
    
        if rank == 0:
            f = open(filename,'w')
            for iw in range(self.Nw):
                energy = iw * self.dw * Hartree
                print(energy, np.imag(g_w2[iw, 0]), np.imag(g_w2[iw, 1]), file=f)
            f.close()

        # Wait for I/O to finish
        self.comm.barrier()


    def check_sum_rule(self, df1_w=None, df2_w=None):
        """Check f-sum rule."""

        if df1_w is None:
            df1_w = self.df1_w
            df2_w = self.df2_w

        N1 = N2 = 0
        for iw in range(self.Nw):
            w = iw * self.dw
            N1 += np.imag(df1_w[iw]) * w
            N2 += np.imag(df2_w[iw]) * w
        N1 *= self.dw * self.vol / (2 * pi**2)
        N2 *= self.dw * self.vol / (2 * pi**2)

        self.printtxt('')
        self.printtxt('Sum rule for ABS:')
        nv = self.nvalence
        self.printtxt('Without local field: N1 = %f, %f  %% error' %(N1, (N1 - nv) / nv * 100) )
        self.printtxt('Include local field: N2 = %f, %f  %% error' %(N2, (N2 - nv) / nv * 100) )

        N1 = N2 = 0
        for iw in range(self.Nw):
            w = iw * self.dw
            N1 -= np.imag(1/df1_w[iw]) * w
            N2 -= np.imag(1/df2_w[iw]) * w
        N1 *= self.dw * self.vol / (2 * pi**2)
        N2 *= self.dw * self.vol / (2 * pi**2)
                
        self.printtxt('')
        self.printtxt('Sum rule for EELS:')
        nv = self.nvalence
        self.printtxt('Without local field: N1 = %f, %f  %% error' %(N1, (N1 - nv) / nv * 100) )
        self.printtxt('Include local field: N2 = %f, %f  %% error' %(N2, (N2 - nv) / nv * 100) )


    def get_macroscopic_dielectric_constant(self, xc='RPA'):
        """Calculate macroscopic dielectric constant. Returns eM1 and eM2

        Macroscopic dielectric constant is defined as the real part of dielectric function at w=0.
        
        Parameters:

        eM1: float
            Dielectric constant without local field correction. (RPA, ALDA)
        eM2: float
            Dielectric constant with local field correction.

        """

        assert self.optical_limit
        self.printtxt('')
        self.printtxt('%s Macroscopic Dielectric Constant:' % xc)
        dirstr = ['x', 'y', 'z']

        for dir in range(3):
        
            eM = np.zeros(2)
            df1, df2 = self.get_dielectric_function(xc=xc, dir=dir)
            eps0 = np.real(df1[0])
            eps = np.real(df2[0])
            self.printtxt('  %s direction' %(dirstr[dir]))
            self.printtxt('    Without local field: %f' % eps0 )
            self.printtxt('    Include local field: %f' % eps )        
            
        return eps0, eps


    def get_absorption_spectrum(self, filename='Absorption.dat'):
        """Calculate optical absorption spectrum. By default, generate a file 'Absorption.dat'.

        Optical absorption spectrum is obtained from the imaginary part of dielectric function.
        """

        assert self.optical_limit
        for dir in range(3):
            df1, df2 = self.get_dielectric_function(xc='RPA', dir=dir)
            if self.xc == 'ALDA':
                df3, df4 = self.get_dielectric_function(xc='ALDA', dir=dir)
            if self.xc == 'ALDA_X':
                df3, df4 = self.get_dielectric_function(xc='ALDA_X', dir=dir)
    
            Nw = df1.shape[0]
    
            if self.xc == 'Bootstrap':
                # bootstrap doesnt support all direction spectra yet
                from gpaw.response.fxc import Bootstrap
                Kc_GG = np.zeros((self.npw, self.npw))
                q_c = np.diag((1, 1, 1))[dir] * self.qopt

                for iG in range(self.npw):
                    qG = np.dot(q_c + self.Gvec_Gc[iG], self.bcell_cv)
                    Kc_GG[iG,iG] = 4 * pi / np.dot(qG, qG)
    
                from gpaw.mpi import world
                assert self.wcomm.size == world.size
                df3 = Bootstrap(self.chi0_wGG, Nw, Kc_GG, self.printtxt, self.print_bootstrap, self.wcomm)
    
            if rank == 0:
                f = open('%s.%s' % (filename, 'xyz'[dir]), 'w')
                #f = open(filename+'.%s'%(dirstr[dir]),'w') # ????
                for iw in range(Nw):
                    energy = iw * self.dw * Hartree
                    if self.xc == 'RPA':
                        print(energy, np.real(df1[iw]), np.imag(df1[iw]), \
                              np.real(df2[iw]), np.imag(df2[iw]), file=f)
                    elif self.xc == 'ALDA':
                        print(energy, np.real(df1[iw]), np.imag(df1[iw]), \
                          np.real(df2[iw]), np.imag(df2[iw]), \
                          np.real(df3[iw]), np.imag(df3[iw]), \
                          np.real(df4[iw]), np.imag(df4[iw]), file=f)
                    elif self.xc == 'Bootstrap':
                        print(energy, np.real(df1[iw]), np.imag(df1[iw]), \
                          np.real(df2[iw]), np.imag(df2[iw]), \
                          np.real(df3[iw]), np.imag(df3[iw]), file=f)
                f.close()
    
            # Wait for I/O to finish
            self.comm.barrier()


    def get_EELS_spectrum(self, filename='EELS.dat'):
        """Calculate EELS spectrum. By default, generate a file 'EELS.dat'.

        EELS spectrum is obtained from the imaginary part of the inverse of dielectric function.
        """

        # calculate RPA dielectric function
        df1, df2 = self.get_dielectric_function(xc='RPA')
        if self.xc == 'ALDA':
            df3, df4 = self.get_dielectric_function(xc='ALDA')
        Nw = df1.shape[0]

        if rank == 0:
            f = open(filename,'w')
            for iw in range(self.Nw):
                energy = iw * self.dw * Hartree
                if self.xc == 'RPA':
                    print(energy, -np.imag(1./df1[iw]), -np.imag(1./df2[iw]), file=f)
                elif self.xc == 'ALDA':
                    print(energy, -np.imag(1./df1[iw]), -np.imag(1./df2[iw]), \
                       -np.imag(1./df3[iw]), -np.imag(1./df4[iw]), file=f)
            f.close()

        # Wait for I/O to finish
        self.comm.barrier()


    def get_jdos(self, f_skn, e_skn, kd, kq, dw, Nw, sigma):
        """Calculate Joint density of states"""

        JDOS_w = np.zeros(Nw)
        nbands = f_skn[0].shape[1]

        for k in range(kd.nbzkpts):
            print(k)
            ibzkpt1 = kd.bz2ibz_k[k]
            ibzkpt2 = kd.bz2ibz_k[kq[k]]
            for n in range(nbands):
                for m in range(nbands):
                    focc = f_skn[0][ibzkpt1, n] - f_skn[0][ibzkpt2, m]
                    w0 = e_skn[0][ibzkpt2, m] - e_skn[0][ibzkpt1, n]
                    if focc > 0 and w0 >= 0:
                        w0_id = int(w0 / dw)
                        if w0_id + 1 < Nw:
                            alpha = (w0_id + 1 - w0/dw) / dw
                            JDOS_w[w0_id] += focc * alpha
                            alpha = (w0/dw-w0_id) / dw
                            JDOS_w[w0_id+1] += focc * alpha
                            
        w = np.arange(Nw) * dw * Hartree

        return w, JDOS_w


    def calculate_induced_density(self, q, w):
        """ Evaluate induced density for a certain q and w.

        Parameters:

        q: ndarray
            Momentum tranfer at reduced coordinate.
        w: scalar
            Energy (eV).
        """

        if type(w) is int:
            iw = w
            w = self.wlist[iw] / Hartree
        elif type(w) is float:
            w /= Hartree
            iw = int(np.round(w / self.dw))
        else:
            raise ValueError('Frequency not correct !')

        self.printtxt('Calculating Induced density at q, w (iw)')
        self.printtxt('(%f, %f, %f), %f(%d)' %(q[0], q[1], q[2], w*Hartree, iw))

        # delta_G0
        delta_G = np.zeros(self.npw)
        delta_G[0] = 1.

        # coef is (q+G)**2 / 4pi
        coef_G = np.zeros(self.npw)
        for iG in range(self.npw):
            qG = np.dot(q + self.Gvec_Gc[iG], self.bcell_cv)
            coef_G[iG] = np.dot(qG, qG)
        coef_G /= 4 * pi

        # obtain chi_G0(q,w)
        dm_wGG = self.get_RPA_dielectric_matrix()
        tmp_GG = dm_wGG[iw]
        del dm_wGG
        chi_G = (np.linalg.inv(tmp_GG)[:, 0] - delta_G) * coef_G

        gd = self.gd
        r = gd.get_grid_point_coordinates()

        # calculate dn(r,q,w)
        drho_R = gd.zeros(dtype=complex)
        for iG in range(self.npw):
            qG = np.dot(q + self.Gvec_Gc[iG], self.bcell_cv)
            qGr_R = np.inner(qG, r.T).T
            drho_R += chi_G[iG] * np.exp(1j * qGr_R)

        # phase = sum exp(iq.R_i)
        return drho_R


    def get_induced_density_z(self, q, w):
        """Get induced density on z axis (summation over xy-plane). """

        drho_R = self.calculate_induced_density(q, w)

        drho_z = np.zeros(self.gd.N_c[2],dtype=complex)
#        dxdy = np.cross(self.h_c[0], self.h_c[1])

        for iz in range(self.gd.N_c[2]):
            drho_z[iz] = drho_R[:,:,iz].sum()

        return drho_z

    def get_eigenmodes(self,filename = None, chi0 = None, calc = None, dm = None, 
                       xc = 'RPA', sum = None, vcut = None, checkphase = False, 
                       return_full = False):
        """
        Calculate the plasmonic eigenmodes as eigenvectors of the dielectric matrix.  

        Parameters:

        filename:  pckl file
                   output from response calculation.
         
        chi0:  gpw file
               chi0_wGG from response calculation.

        calc:  gpaw calculator instance
               ground state calculator used in response calculation.
               Wavefunctions only needed if chi0 is calculated from scratch

        dm:  gpw file
             dielectric matrix from response calculation

        xc:  str 'RPA'or 'ALDA' XC- Kernel
        
        sum:  str
              '2D': sum in the x and y directions
              '1D': To be implemented

        vcut:  str '0D','1D' or '2D'
               Cut the Coulomb potential 

        checkphase:   Bool
                      if True, the eigenfunctions id rotated in the complex
                      plane, to be made as real as posible

        return_full:  Bool
                      if True, the eigenvectors in reciprocal space is also
                      returned. 
           
        """
        self.read(filename)
        self.pbc = [1,1,1]
        #self.calc.atoms.pbc = [1,1,1]
        npw = self.npw
        self.w_w = np.linspace(0, self.dw * (self.Nw - 1)*Hartree, self.Nw)
        self.vcut = vcut
        dm_wGG = self.get_dielectric_matrix(xc=xc,
                                            symmetric=False,
                                            chi0_wGG=chi0,
                                            calc=calc,
                                            vcut=vcut)
    
        q = self.q_c
        

        # get grid on which the eigenmodes are calculated
        #gd = self.calc.wfs.gd
        #r = gd.get_grid_point_coordinates()
        #rrr = r*Bohr 
        from gpaw.utilities.gpts import get_number_of_grid_points
        from gpaw.grid_descriptor import GridDescriptor
        grid_size = [1,1,1]
        h=0.2
        cell_cv = self.acell_cv*np.diag(grid_size)
        mode = 'fd'
        realspace = True
        h /= Bohr
        N_c = get_number_of_grid_points(cell_cv, h, mode, realspace)
        gd = GridDescriptor(N_c, cell_cv, self.pbc) 
        #gd = self.calc.wfs.gd
        r = gd.get_grid_point_coordinates()
        rrr = r*Bohr
        
        eig_0 = np.array([], dtype = complex)
        eig_left = np.array([], dtype = complex)
        eig_right = np.array([], dtype = complex)
        vec_modes = np.zeros([1, self.npw], dtype = complex)
        vec_modes_dual = np.zeros([1, self.npw], dtype = complex)
        vec_modes_density = np.zeros([1, self.npw], dtype = complex)
        vec_modes_norm = np.zeros([1, self.npw], dtype = complex)
        eig_all = np.zeros([1, self.npw], dtype = complex)
        eig_dummy = np.zeros([1, self.npw], dtype = complex)
        v_dummy = np.zeros([1, self.npw], dtype = complex)
        vec_dummy = np.zeros([1, self.npw], dtype = complex)
        vec_dummy2 = np.zeros([1, self.npw], dtype = complex)
        w_0 = np.array([]) 
        w_left = np.array([])
        w_right = np.array([])
     
        if sum == '2D':
            v_ind = np.zeros([1, r.shape[-1]], dtype = complex)
            n_ind = np.zeros([1, r.shape[-1]], dtype = complex)
        elif sum == '1D':            
            self.printtxt('1D sum not implemented')
            return 
        else:
            v_ind = np.zeros([1, r.shape[1], r.shape[2], r.shape[3]], dtype = complex)
            n_ind = np.zeros([1, r.shape[1], r.shape[2], r.shape[3]], dtype = complex)

        eps_GG_plus = dm_wGG[0]
        eig_plus, vec_plus = np.linalg.eig(eps_GG_plus)  # find eigenvalues and eigenvectors
        vec_plus_dual = np.linalg.inv(vec_plus)
        
        # loop over frequencies, where the eigenvalues for the 2D matrix in G,G' are found.
        for i in np.array(range(self.Nw-1))+1:
            eps_GG = eps_GG_plus
            eig, vec = eig_plus,vec_plus
            vec_dual = vec_plus_dual
            eps_GG_plus = dm_wGG[i] # epsilon_GG'(omega + d-omega)
            eig_plus, vec_plus = np.linalg.eig(eps_GG_plus)
            vec_plus_dual = np.linalg.inv(vec_plus)
            eig_dummy[0,:] = eig
            eig_all = np.append(eig_all, eig_dummy, axis=0) # append all eigenvalues to array         
            # loop to check find the eigenvalues that crosses zero from negative to positive values:
            for k in range(self.npw):
                for m in range(self.npw):
                    if eig[k]< 0 and 0 < eig_plus[m]:
                        # check it's the same mode - Overlap between eigenvectors should be large:
                        if abs(np.inner(vec[:,k], vec_plus_dual[m,:])) > 0.95:                             
                            self.printtxt('crossing found at w = %1.1f eV'%self.w_w[i-1])
                            eig_left = np.append(eig_left, eig[k])   
                            eig_right = np.append(eig_right, eig_plus[m])

                            vec_dummy[0, :] = vec[:,k]
                            vec_modes = np.append(vec_modes, vec_dummy, axis = 0)
                            vec_dummy[0, :] = vec_dual[k, :].T
                            vec_modes_dual = np.append(vec_modes_dual, vec_dummy, axis = 0)
                                                   
                            w1 = self.w_w[i-1]
                            w2 = self.w_w[i]
                            a = np.real((eig_plus[m]-eig[k]) / (w2-w1))
                            w0 = np.real(-eig[k]) / a + w1
                            eig0 = a*(w0-w1)+eig[k]

                            w_0 = np.append(w_0,w0)
                            w_left = np.append(w_left, w1)
                            eig_0 = np.append(eig_0,eig0)                           

                            n_dummy = np.zeros([1, r.shape[1], r.shape[2],
                                                r.shape[3]], dtype = complex)
                            v_dummy = np.zeros([1, r.shape[1], r.shape[2],
                                                r.shape[3]], dtype = complex)    

                            vec_n = np.zeros([self.npw])
                            
                            for iG in range(self.npw):  # Fourier transform
                                qG = np.dot((q + self.Gvec_Gc[iG]), self.bcell_cv)
                                coef_G = np.dot(qG, qG) / (4 * pi)
                                qGr_R = np.inner(qG, r.T).T
                                v_dummy += vec[iG, k] * np.exp(1j * qGr_R) 
                                n_dummy += vec[iG, k] * np.exp(1j * qGr_R) * coef_G
                                                        
                            if checkphase: # rotate eigenvectors in complex plane 
                                integral = np.zeros([81])
                                phases = np.linspace(0,2,81)
                                for ip in range(81):
                                    v_int = v_dummy * np.exp(1j * pi * phases[ip])
                                    integral[ip] = abs(np.imag(v_int)).sum()                                     
                                phase = phases[np.argsort(integral)][0]
                                
                                v_dummy *= np.exp(1j * pi * phase)
                                n_dummy *= np.exp(1j * pi * phase)
                                

                            if sum == '2D':
                                i_xyz = 3 
                                v_dummy_z = np.zeros([1,v_dummy.shape[i_xyz]],
                                                     dtype = complex)
                                n_dummy_z = np.zeros([1,v_dummy.shape[i_xyz]],
                                                     dtype = complex)
                                v_dummy_z[0,:] = np.sum(np.sum(v_dummy, axis = 1),
                                                        axis = 1)[0,:]
                                n_dummy_z[0,:] = np.sum(np.sum(n_dummy, axis = 1),
                                                        axis = 1)[0,:]

                                v_ind = np.append(v_ind, v_dummy_z, axis=0)
                                n_ind = np.append(n_ind, n_dummy_z, axis=0)
                                                    
                            elif sum == '1D':
                                self.printtxt('1D sum not implemented')
                            else :
                                v_ind = np.append(v_ind, v_dummy, axis=0)
                                n_ind = np.append(n_ind, n_dummy, axis=0)
            
        """                        
        returns: grid points, frequency grid, all eigenvalues, mode energies, left point energies,
                 mode eigenvalues, eigenvalues of left and right-side points,
                 (mode eigenvectors, mode dual eigenvectors,)
                 induced potential in real space, induced density in real space
        """
        
        if return_full:
            return rrr, self.w_w, eig_all[1:], w_0, eig_0, w_left, eig_left, \
               eig_right, vec_modes[1:], vec_modes_dual[1:], v_ind[1:], n_ind[1:]
            
        else:
            return rrr, self.w_w, eig_all[1:], w_0, eig_0, w_left, eig_left, \
               eig_right, v_ind[1:], n_ind[1:]
                        
        
    def project_chi_to_LCAO_pair_orbital(self, orb_MG):

        nLCAO = orb_MG.shape[0]
        N = np.zeros((self.Nw, nLCAO, nLCAO), dtype=complex)

        kcoulinv_GG = np.zeros((self.npw, self.npw))
        for iG in range(self.npw):
            qG = np.dot(self.q_c + self.Gvec_Gc[iG], self.bcell_cv)
            kcoulinv_GG[iG, iG] = np.dot(qG, qG)

        kcoulinv_GG /= 4.*pi

        dm_wGG = self.get_RPA_dielectric_matrix()

        for mu in range(nLCAO):
            for nu in range(nLCAO):
                pairorb_R = orb_MG[mu] * orb_MG[nu]
                if not (pairorb_R * pairorb_R.conj() < 1e-10).all():
                    tmp_G = np.fft.fftn(pairorb_R) * self.vol / self.nG0

                    pairorb_G = np.zeros(self.npw, dtype=complex)
                    for iG in range(self.npw):
                        index = self.Gindex[iG]
                        pairorb_G[iG] = tmp_G[index[0], index[1], index[2]]

                    for iw in range(self.Nw):
                        chi_GG = (dm_wGG[iw] - np.eye(self.npw)) * kcoulinv_GG
                        N[iw, mu, nu] = (np.outer(pairorb_G.conj(), pairorb_G) * chi_GG).sum()
#                        N[iw, mu, nu] = np.inner(pairorb_G.conj(),np.inner(pairorb_G, chi_GG))

        return N


    def write(self, filename, all=False):
        """Dump essential data"""

        data = {'nbands': self.nbands,
                'acell': self.acell_cv, #* Bohr,
                'bcell': self.bcell_cv, #/ Bohr,
                'h_cv' : self.gd.h_cv,   #* Bohr,
                'nG'   : self.gd.N_c,
                'nG0'  : self.nG0,
                'vol'  : self.vol,   #* Bohr**3,
                'BZvol': self.BZvol, #/ Bohr**3,
                'nkpt' : self.kd.nbzkpts,
                'ecut' : self.ecut,  #* Hartree,
                'npw'  : self.npw,
                'eta'  : self.eta,   #* Hartree,
                'ftol' : self.ftol,
                'Nw'   : self.Nw,
                'NwS'  : self.NwS,
                'dw'   : self.dw,    # * Hartree,
                'q_red': self.q_c,
                'q_car': self.qq_v,    # / Bohr,
                'qmod' : np.dot(self.qq_v, self.qq_v), # / Bohr
                'vcut' : self.vcut,
                'pbc'  : self.pbc,
                'nvalence'     : self.nvalence,                
                'hilbert_trans' : self.hilbert_trans,
                'optical_limit' : self.optical_limit,
                'e_skn'         : self.e_skn,          # * Hartree,
                'f_skn'         : self.f_skn,
                'bzk_kc'       : self.kd.bzk_kc,
                'ibzk_kc'      : self.kd.ibzk_kc,
                'kq_k'         : self.kq_k,
                'Gvec_Gc'      : self.Gvec_Gc,
                'dfNLFRPA_w'   : self.df1_w,
                'dfLFCRPA_w'   : self.df2_w,
                'dfNLFALDA_w'  : self.df3_w,
                'dfLFCALDA_w'  : self.df4_w,
                'df_flag'      : True}

        if all:
            from gpaw.response.parallel import par_write
            par_write('chi0' + filename,'chi0_wGG',self.wcomm,self.chi0_wGG)
        
        if rank == 0:
            pickle.dump(data, open(filename, 'w'), -1)

        self.comm.barrier()


    def read(self, filename):
        """Read data from pickle file"""

        data = pickle.load(open(filename))
        
        self.nbands = data['nbands']
        self.acell_cv = data['acell']
        self.bcell_cv = data['bcell']
        self.nG0   = data['nG0']
        self.vol   = data['vol']
        self.BZvol = data['BZvol']
        self.ecut  = data['ecut']
        self.npw   = data['npw']
        self.eta   = data['eta']
        self.ftol  = data['ftol']
        self.Nw    = data['Nw']
        self.NwS   = data['NwS']
        self.dw    = data['dw']
        self.q_c   = data['q_red']
        self.qq_v  = data['q_car']
        self.qmod  = data['qmod']
        #self.vcut  = data['vcut']
        #self.pbc = data['pbc']
        
        self.hilbert_trans = data['hilbert_trans']
        self.optical_limit = data['optical_limit']
      
        self.e_skn  = data['e_skn']
        self.f_skn  = data['f_skn']
                
        self.nvalence= data['nvalence']
        self.kq_k    = data['kq_k']
        self.Gvec_Gc  = data['Gvec_Gc']
        self.df1_w   = data['dfNLFRPA_w']
        self.df2_w   = data['dfLFCRPA_w']
        self.df3_w   = data['dfNLFALDA_w']
        self.df4_w   = data['dfLFCALDA_w']
        self.df_flag = data['df_flag']
        
        self.printtxt('Read succesfully !')
