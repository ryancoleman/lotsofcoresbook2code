import numpy as np
from ase.units import Bohr
from ase.transport.tools import dagger

from gpaw.transport.selfenergy import LeadSelfEnergy
from gpaw.transport.greenfunction import GreenFunction
#from ase.transport.selfenergy import LeadSelfEnergy
#from ase.transport.greenfunction import GreenFunction

from gpaw.fd_operators import Laplace
from gpaw.grid_descriptor import GridDescriptor
from gpaw.lfc import NewLocalizedFunctionsCollection as LFC
import gpaw.mpi as mpi
from gpaw.utilities.tools import tri2full

''' differ from version 4 by selected scaning range in z direction by values in Angstrom '''

class LocalizedFunctions:

    def __init__(self, gd, f_iG, corner_c, index=None, vt_G=None):
        self.gd = gd
        #assert not gd.is_non_orthogonal()
        self.size_c = np.array(f_iG.shape[1:4])
        self.f_iG = f_iG
        self.corner_c = corner_c
        self.index = index
        self.vt_G = vt_G
   
    def periodic(self):
        trans_c = [np.array([0, 0, 0]), np.array([0, 1, 0]), np.array([-1,  0,  0]), np.array([ 0, -1,  0]), 
                   np.array([1, 1, 0]), np.array([ 1, -1,  0]), np.array([-1,  1,  0]), np.array([-1, -1,  0]), 
                   np.array([ 1,  0,  0])]
        trans_c[:] *= np.array([self.gd.N_c[0],self.gd.N_c[1],0])
        
        self.periodic_list = trans_c+self.corner_c
 
        list = []
        for corner in self.periodic_list:
            list.append(LocalizedFunctions(self.gd,self.f_iG,
                                        corner_c=corner,
                                        index=self.index,
                                        vt_G=self.vt_G))
        return list



    def __len__(self):
        return len(self.f_iG)

    def apply_t(self):
        """Apply kinetic energy operator and return new object."""
        p = 2  # padding
        newsize_c = self.size_c + 2 * p
        gd = GridDescriptor(N_c=newsize_c + 1,
                            cell_cv=self.gd.h_c * (newsize_c + 1),
                            pbc_c=False,
                            comm=mpi.serial_comm)
        T = Laplace(gd, scale =1/2., n=p)
        f_ig = np.zeros((len(self.f_iG),) + tuple(newsize_c))
        f_ig[:, p:-p, p:-p, p:-p] = self.f_iG
        Tf_iG = np.empty_like(f_ig)
        T.apply(f_ig, Tf_iG)
        return LocalizedFunctions(self.gd, Tf_iG, self.corner_c - p,
                                  self.index)
        
    def overlap(self, other):
        start_c = np.maximum(self.corner_c, other.corner_c)
        stop_c = np.minimum(self.corner_c + self.size_c,
                            other.corner_c + other.size_c)
        if (start_c < stop_c).all():
            astart_c = start_c - self.corner_c
            astop_c = stop_c - self.corner_c
            a_iG = self.f_iG[:,
                astart_c[0]:astop_c[0],
                astart_c[1]:astop_c[1],
                astart_c[2]:astop_c[2]].copy().reshape((len(self.f_iG), -1))
            bstart_c = start_c - other.corner_c
            bstop_c = stop_c - other.corner_c
            b_iG = other.f_iG[:,
                bstart_c[0]:bstop_c[0],
                bstart_c[1]:bstop_c[1],
                bstart_c[2]:bstop_c[2]].copy().reshape((len(other.f_iG), -1))
            if self.vt_G is not None:
                x = self.vt_G[start_c[0]:stop_c[0],
                              start_c[1]:stop_c[1],
                              start_c[2]:stop_c[2]].copy().reshape((-1,))
                #print 'start_c:',start_c
                #print 'stop_c:',stop_c
                #print 'shape of x', np.shape(x)
                #print 'shape of a_iG', np.shape(a_iG)
                a_iG *= x    
     
            return self.gd.dv * np.inner(a_iG, b_iG)
        else:
            return None
    

    def __or__(self, other):
        if isinstance(other, LocalizedFunctions):
            return self.overlap(other)

        # other is a potential:
        vt_G = other
        return LocalizedFunctions(self.gd, self.f_iG, self.corner_c,
                                  self.index, vt_G)

class WannierFunction(LocalizedFunctions):
    def __init__(self, gd, wanf_G, corner_c, index=None):
        LocalizedFunctions.__init__(self, gd, wanf_G[np.newaxis, :, :, :],
                                    corner_c, index)

class AtomCenteredFunctions(LocalizedFunctions):
    def __init__(self, gd, spline_j, spos_c, index=None):
        rcut = max([spline.get_cutoff() for spline in spline_j])
        corner_c = np.ceil(spos_c * gd.N_c - rcut / gd.h_c).astype(int)
        size_c = np.ceil(spos_c * gd.N_c + rcut / gd.h_c).astype(int) - corner_c
        smallgd = GridDescriptor(N_c=size_c + 1,
                                 cell_cv=gd.h_c * (size_c + 1),
                                 pbc_c=False,
                                 comm=mpi.serial_comm)
        lfc = LFC(smallgd, [spline_j])
        lfc.set_positions((spos_c[np.newaxis, :] * gd.N_c - corner_c + 1) /
                          smallgd.N_c)
        ni = lfc.Mmax
        f_iG = smallgd.zeros(ni)
        lfc.add(f_iG, {0: np.eye(ni)})
        LocalizedFunctions.__init__(self, gd, f_iG, corner_c, 
                                     index=index)


class STM:
    def __init__(self, srf, tip, srf_pl, tip_pl, eta1=1e-4, eta2=1e-4):
        self.tip = tip                     # tip calcualtor
        self.tip_pl = tip_pl               # tip periodic layers calculator
        self.srf = srf                     # surface calculator
        self.srf_pl = srf_pl               # surface periodic layers calculator
        self.eta1 = eta1
        self.eta2 = eta2 
        self.tgd = tip.gd                  
        self.sgd = srf.gd 
        
        #assert tgd.h_cv == sgd.h_cv
        #assert not (tgd.h_c - sgd.h_c).any()

    def initialize(self, srf_atoms, tip_atoms, srf_coupling_atoms, tip_coupling_atoms, 
                   srf_pl_atoms, tip_pl_atoms, energies, bias=0):
        """
        srf_atoms: the atom index of onsite surface atoms 
        tip_atoms: the atom index of onsite tip atoms
        energies: list of energies, for which the transmission function should be evaluated
        bias: will be used to precalculate the greenfunctions of surface and tip
        """
        self.bias = bias
        self.energies = energies
        self.tip_atoms = tip_atoms
        self.srf_atoms = srf_atoms
        self.srf_coupling_atoms = srf_coupling_atoms
        self.tip_coupling_atoms = tip_coupling_atoms
        #align potentials for onsite surface, onsite tip and tip principle layers according to surface_pl
        srf_alignv = self.align_v(calc1 = self.srf_pl, index1=srf_pl_atoms[0], calc2=self.srf, index2=srf_atoms[0])         
        tip_alignv = srf_alignv + self.tip.get_fermi_level() - self.srf.get_fermi_level() #- bias
        tip_pl_alignv = self.align_v(calc1=self.tip, index1=tip_atoms[-1], calc2=self.tip_pl, index2=tip_pl_atoms[-1])
        tip_pl_alignv += tip_alignv
        self.srf_alignv = srf_alignv
        self.tip_alignv = tip_alignv

        #hamiltoian and overlap matrix of onsite tip and srf
        self.h1, self.s1 = Lead_HS(calc = self.tip, atom_list = tip_atoms, shift = srf_alignv-bias).get_HS()
        print('shape of onsite tip', np.shape(self.h1), np.shape(self.s1))
        self.h2, self.s2 = Lead_HS(calc = self.srf, atom_list = srf_atoms, shift = tip_alignv).get_HS()        
        print('shape of onsite srf', np.shape(self.h2), np.shape(self.s2))

        #hamiltoian and overlap matrix of two tip and srf principle layers
        self.h10, self.s10 = Lead_HS(calc = self.tip_pl, atom_list = tip_pl_atoms, shift = tip_pl_alignv).get_HS() #2 and only 2 principle layers
        print('shape of pl tip', np.shape(self.h10), np.shape(self.s10)) 
        self.h20, self.s20 = Lead_HS(calc = self.srf_pl, atom_list = srf_pl_atoms, shift = -bias).get_HS() #2 and only 2 principle layers
        print('shape of pl srf', np.shape(self.h20), np.shape(self.s20))

        nbf1, nbf2 = len(self.h1), len(self.h2)       #No. of basis functions of onsite tip and surface atoms
        pl1, pl2 = len(self.h10)/2, len(self.h20)/2   #No. of basis functions per principle layer of tip and srf
        nenergies = len(energies)        

        #periodic part of the tip
        hs1_dii = self.h10[:pl1, :pl1], self.s10[:pl1, :pl1]
        hs1_dij = self.h10[:pl1, pl1:2*pl1], self.s10[:pl1, pl1:2*pl1]
        #coupling betwen per. and non. per part of the tip
        h1_im = np.zeros((pl1, nbf1), complex)
        s1_im = np.zeros((pl1, nbf1), complex)
        h1_im[:pl1, :pl1], s1_im[:pl1, :pl1] = hs1_dij
        hs1_dim = [h1_im, s1_im]

        #periodic part the surface 
        hs2_dii = self.h20[:pl2, :pl2], self.s20[:pl2, :pl2]
        hs2_dij = self.h20[pl2:2*pl2, :pl2], self.s20[pl2:2*pl2, :pl2]
        #coupling betwen per. and non. per part of the surface
        h2_im = np.zeros((pl2, nbf2), complex)
        s2_im = np.zeros((pl2, nbf2), complex)
        h2_im[-pl2:, -pl2:], s2_im[-pl2:, -pl2:] = hs2_dij
        hs2_dim = [h2_im, s2_im]
        #tip and surface greenfunction 
        self.selfenergy1 = LeadSelfEnergy(hs1_dii, hs1_dij, hs1_dim)
        self.selfenergy2 = LeadSelfEnergy(hs2_dii, hs2_dij, hs2_dim)
        self.greenfunction1 = GreenFunction(self.h1, self.s1,[self.selfenergy1])
        self.greenfunction2 = GreenFunction(self.h2, self.s2,[self.selfenergy2])
        #print 'shape of greenfunction tip:', np.shape(self.greenfunction1)
        #print 'shape of greenfunction srf:', np.shape(self.greenfunction2)
        #shift the bands due to the bias
        self.selfenergy1.set_bias(0.0)
        self.selfenergy2.set_bias(bias)

        # extend the surface potentials and basis functions
        self.setup_stm()


        #tip and surface greenfunction matrices.
        nbf1_small = self.ni #XXX Change this for efficiency in the future
        nbf2_small = self.nj #XXX -||-
        coupling_list1 = range(nbf1-nbf1_small,nbf1)# XXX -||-
        coupling_list2 = range(0,nbf2_small)# XXX -||-
        self.gft1_emm = np.zeros((nenergies, nbf1_small, nbf1_small), complex)
        self.gft2_emm = np.zeros((nenergies, nbf2_small, nbf2_small), complex)

        for e, energy in enumerate(self.energies):
            gft1_mm = self.greenfunction1(energy)[coupling_list1]
            gft1_mm = np.take(gft1_mm, coupling_list1, axis=1)

            gft2_mm = self.greenfunction2(energy)[coupling_list2]
            gft2_mm = np.take(gft2_mm, coupling_list2, axis=1)

            self.gft1_emm[e] = gft1_mm
            self.gft2_emm[e] = gft2_mm

    def scan(self,height=None,zmin=None,zmax=None):
        if height is not None:
            self.current = np.zeros((self.srf.gd.N_c[0],self.srf.gd.N_c[1]))
            h = 0
            h = np.round(height/Bohr/self.srf.gd.h_c[2]).astype(int)
            for px in range(self.srf.gd.N_c[0]):
                for py in range(self.srf.gd.N_c[1]):
                    pos = np.array([px,py,h])
                    V_12 = self.get_Vts(pos)
                    I = self.get_current(V_12)
                    self.current[px,py] = I
                    print('I = ', I)
            return self.current
        if height is None:
            zmin_c = np.round(zmin/Bohr/self.srf.gd.h_c[2]).astype(int)
            zmax_c = np.round(zmax/Bohr/self.srf.gd.h_c[2]).astype(int)
            nz = zmax_c - zmin_c
            self.current = np.zeros((self.srf.gd.N_c[0],self.srf.gd.N_c[1],nz))
            for px in range(self.srf.gd.N_c[0]):
                for py in range(self.srf.gd.N_c[1]):
                    for pz in range(zmin_c,zmax_c):
                        # get coupleing between tip and surface Vts for a given tip position on original srf grids 
                        pos = np.array([px,py,pz])
                        V_12 = self.get_Vts(pos)
                        # get current values
                        I = self.get_current(V_12)
                        self.current[px,py,pz-zmin_c] = I
                        print(I)
            return self.current



    def get_current(self, v_12):
        # get transmission
        T_e = self.get_transmission(v_12)
        I = -np.trapz(x = self.energies, y = T_e)         
        return I



    def get_transmission(self, v_12):
        nenergies = len(self.energies)
        T_e = np.empty(nenergies,float)
        v_21 = dagger(v_12)
        for e, energy in enumerate(self.energies):
            gft1 = self.gft1_emm[e]
            gft2 = self.gft2_emm[e]            
            a1 = (gft1 - dagger(gft1))
            a2 = (gft2 - dagger(gft2))
            v12_a2 = np.dot(v_12, a2)
            v21_a1 = np.dot(v_21, a1)
            T = -np.trace(np.dot(v12_a2, v21_a1))
            T_e[e] = T
            self.T_e = T_e
        return T_e



    def get_Vts(self, tip_pos):
        # tip_pos is in terms of grid points, showing the tip apex position respect to the surface grid
        print('tip_pos is', tip_pos)
        tip_apex_c = self.tip_apex_c[0] 
        srf_pos_av = self.srf.atoms.get_positions() / Bohr
        srf_zmax = srf_pos_av[:,2].max()
        srf_zmax_c = np.round(srf_zmax / self.srf.gd.h_c[2]).astype(int)
        srf_ox = self.srf_extension[0] + tip_pos[0]
        srf_oy = self.srf_extension[1] + tip_pos[1]
        srf_oz = srf_zmax_c + tip_pos[2]
        shift = np.array([srf_ox, srf_oy, srf_oz])
        shift -= tip_apex_c
        for tbox in self.tip_coupling_functions:
            tbox.corner_c += shift
        for tbox in self.tip_coupling_kinetics:
            tbox.corner_c += shift
        # get the combined effective potential on the extended surface grid
        """
        This part needs to be updated for a separate xc calculation
        """
        tip_ox = self.ox
        tip_oy = self.oy
        tip_oz = self.oz
        srf_sx = srf_ox - self.nmx
        srf_ex = srf_ox + self.npx
        srf_sy = srf_oy - self.nmy
        srf_ey = srf_oy + self.npy
        srf_sz = srf_oz - self.nmz
        srf_ez = srf_oz + self.npz
        tip_sx = tip_ox - self.nmx
        tip_ex = tip_ox + self.npx
        tip_sy = tip_oy - self.nmy
        tip_ey = tip_oy + self.npy
        tip_sz = tip_oz - self.nmz
        tip_ez = tip_oz + self.npz
        V_tip = self.tip.get_effective_potential() - self.tip_alignv 
        V_couple = self.extvt_G.copy()
        V_couple[srf_sx:srf_ex, srf_sy:srf_ey, srf_sz:srf_ez] = V_tip[tip_sx:tip_ex, tip_sy:tip_ey, tip_sz:tip_ez]

        #get V_ts 
        V_ts = np.zeros((self.ni,self.nj))
        for ns in range(len(self.srf_coupling_functions)):
            sf = self.srf_coupling_functions[ns]
            sk = self.srf_coupling_kinetics[ns]
            j1 = sf[0].index
            j2 = j1 + len(sf[0])
            for t in self.tip_coupling_functions:
                i1 = t.index
                i2 = i1 + len(t)
                test_overlap1 = []
                test_overlap2 = []
                for x in range(9):
                    overlap = (t | sf[x])                    
                    if overlap is not None:
                        test_overlap1.append(overlap.copy())
                        test_overlap2.append((np.abs(overlap)).sum())
                    else:
                        test_overlap1.append(overlap)
                        test_overlap2.append(overlap)
                k = np.argmax(test_overlap2)
                if test_overlap1[k] is not None:
                    V_ts[i1:i2, j1:j2] += (t|sk[k])
                    V_ts[i1:i2, j1:j2] += (t|V_couple|sf[k])

        # return the tbox back :) 
        for tbox in self.tip_coupling_functions:
            tbox.corner_c -= shift
        for tbox in self.tip_coupling_kinetics:
            tbox.corner_c -= shift
       
        #print 'np.shape(V_ts)',np.shape(V_ts)
        return V_ts       

        
    def setup_stm(self):
        tip_pos_av = self.tip.atoms.get_positions() / Bohr
        srf_pos_av = self.srf.atoms.get_positions() / Bohr
        
        # get selected tip functions on small grids and apply kinetics
        self.tip_coupling_functions = []
        self.tip_coupling_kinetics = []
        i = 0
        for a in self.tip_coupling_atoms:
            setup = self.tip.wfs.setups[a]
            spos_c = tip_pos_av[a] / self.tip.gd.cell_c            #gives the scaled position of tip atoms in tip cell
            for phit in setup.phit_j:
                f = AtomCenteredFunctions(self.tip.gd, [phit], spos_c, i)
                f_kinetic = f.apply_t()
                self.tip_coupling_functions.append(f)
                self.tip_coupling_kinetics.append(f_kinetic)
                i += len(f.f_iG)
        self.ni = i
        #print self.tip_coupling_functions[1].size_c
        #print self.tip_coupling_functions[1].corner_c
        #print np.shape(self.tip_coupling_functions[1].f_iG)

        # get tip cutoff and surface extension
        tip_corner_c = []
        temp_min = self.tip_coupling_kinetics[0].corner_c
        temp_max = self.tip_coupling_kinetics[0].corner_c
        for f in self.tip_coupling_kinetics:
                start_c = np.minimum(f.corner_c, temp_min)
                end_c = np.maximum(f.corner_c+f.size_c, temp_max)
                temp_min = start_c
                temp_max = end_c
        #print "start_c",start_c         # in terms of grid points
        #print "end_c",end_c             # in terms of grid points
        self.tip_zmin = tip_pos_av[:,2].min()
        self.tip_apex_index = np.where(tip_pos_av[:,2] < self.tip_zmin + 0.01)
        self.tip_apex_c = np.round(tip_pos_av[self.tip_apex_index] / self.tgd.h_c).astype(int)
        self.ox, self.oy, self.oz = self.tip_apex_c[0]
        self.npx, self.npy, self.npz = np.round(end_c).astype(int) - self.tip_apex_c[0]
        self.nmx, self.nmy, self.nmz = self.tip_apex_c[0] - np.round(start_c).astype(int)

        # get the extension of srf basis functions and effective potentials
        srf_extension = np.array([self.srf.gd.N_c[0], self.srf.gd.N_c[1], 0])
        self.srf_extension = srf_extension
        print('extension of surface in terms of grids:', srf_extension)
        
        # get surface fucntions on small grid with periodic conditions with entension 
        self.srf_coupling_functions = []
        self.srf_coupling_kinetics = []
        j = 0
        for a in self.srf_coupling_atoms: 
            setup = self.srf.wfs.setups[a]
            spos_c = srf_pos_av[a] / self.srf.gd.cell_c
            for phit in setup.phit_j:
                f = AtomCenteredFunctions(self.srf.gd, [phit], spos_c, j)
                f.corner_c += srf_extension
                f_kinetic = f.apply_t()
                self.srf_coupling_functions.append(f.periodic())
                self.srf_coupling_kinetics.append(f_kinetic.periodic())
                j += len(f.f_iG)
        self.nj = j
        #print np.shape(self.srf_coupling_functions[0][0].f_iG), self.srf_coupling_functions[0][0].corner_c
        #print np.shape(self.srf_coupling_functions[0][8].f_iG), self.srf_coupling_functions[0][8].corner_c
        # Extend the surface grid and surface effective potential
        svt_G = self.srf.get_effective_potential() - self.srf_alignv + self.bias 
        ex,ey,ez = srf_extension
        newsize_c = 2 * srf_extension + self.sgd.N_c
        extvt_G = np.zeros(newsize_c)
        extvt_G[ex:-ex,ey:-ey,:] = svt_G
        extvt_G[:ex,ey:-ey,:]  = svt_G[-ex:,:,:]
        extvt_G[-ex:,ey:-ey,:] = svt_G[:ex,:,:]
        extvt_G[:,:ey,:] = extvt_G[:,-2*ey:-ey,:]
        extvt_G[:,-ey:,:] = extvt_G[:,ey:2*ey,:]
        self.extvt_G = extvt_G  

        extsgd = GridDescriptor(N_c=newsize_c,
                                cell_cv=self.sgd.h_c*(newsize_c),
                                pbc_c=False,
                                comm=mpi.serial_comm)

        ## Transfer the coner_c of surface basis-functions and srf_kinetics to the extended surface cell
        #for sbox in self.srf_coupling_functions:
        #    sbox.corner_c += srf_extension
        #for sbox in self.srf_coupling_kinetics:
        #    sbox.corner_c += srf_extension



    def align_v(self,calc1,index1,calc2,index2):
        pos1 = calc1.atoms.positions[index1]
        pos1_c = np.round(pos1/Bohr/calc1.wfs.gd.h_c).astype(int)
        pos1_v = calc1.get_effective_potential()[pos1_c[0],pos1_c[1],pos1_c[2]]
        pos2 = calc2.atoms.positions[index2]
        pos2_c = np.round(pos2/Bohr/calc2.wfs.gd.h_c).astype(int)
        pos2_v = calc2.get_effective_potential()[pos2_c[0],pos2_c[1],pos2_c[2]]
        return pos2_v - pos1_v


class Lead_HS:
    def __init__(self, calc, atom_list, extension = None, shift = 0.0):
        self.calc = calc
        self.gd = calc.wfs.gd
        self.atom_list = atom_list
        self.shift = shift
        self.extension = extension
        if self.extension is None:
            self.extension = self.gd.N_c     ####################### needs to be modified later #######################
 
    def get_HS(self):
        # get functions on small grid without periodic conditions
        self.basis_functions_fix = []
        pos_av = self.calc.atoms.get_positions() / Bohr

        # get functions on small grid with periodic conditions without entension 
        self.basis_functions = []
        self.basis_kinetics = []
        Hi = 0
        for a in self.atom_list:
            setup = self.calc.wfs.setups[a]
            spos_c = pos_av[a] / self.gd.cell_c
            for phit in setup.phit_j:
                f = AtomCenteredFunctions(self.gd, [phit], spos_c, Hi)
                f.corner_c += self.extension
                f_kinetic = f.apply_t()
                self.basis_functions_fix.append(f)
                self.basis_functions.append(f.periodic())
                self.basis_kinetics.append(f_kinetic.periodic())
                Hi += len(f.f_iG)
        self.nHi = Hi
        #print np.shape(self.srfH_kinetics[0].f_iG)
        # extend the effective potential for surface
        svt_G = petential = self.calc.get_effective_potential() - self.shift
        extension = self.extension
        ex,ey,ez = extension
        newsize_c = 2 * extension + self.gd.N_c
        ext_svt_G = np.zeros(newsize_c)
        ext_svt_G[ex:-ex,ey:-ey,ez:-ez] = svt_G
        ext_svt_G[:ex,ey:-ey,ez:-ez]  = svt_G[-ex:,:,:]
        ext_svt_G[-ex:,ey:-ey,ez:-ez] = svt_G[:ex,:,:]
        ext_svt_G[:,:ey,ez:-ez] = ext_svt_G[:,-2*ey:-ey,ez:-ez]
        ext_svt_G[:,-ey:,ez:-ez] = ext_svt_G[:,ey:2*ey,ez:-ez]
        ext_svt_G[:,:,:ez] = ext_svt_G[:,:,-2*ez:-ez]
        ext_svt_G[:,:,-ez:] = ext_svt_G[:,:,ez:2*ez]

        # get H 
        self.H = np.zeros((self.nHi, self.nHi))
        self.S = np.zeros((self.nHi, self.nHi))
        for ns in range(len(self.basis_functions)):
            sf = self.basis_functions[ns]
            sk = self.basis_kinetics[ns]
            j1 = sf[0].index
            j2 = j1 + len(sf[0])

            for s in self.basis_functions_fix:
                i1 = s.index
                i2 = i1 + len(s)
                test_overlap1 = []
                test_overlap2 = []

                for x in range(9):
                    overlap = (s | sf[x])
                    if overlap is not None:
                        test_overlap1.append(overlap.copy())
                        test_overlap2.append((np.abs(overlap)).sum())
                    else:
                        test_overlap1.append(overlap)
                        test_overlap2.append(overlap)
                k = np.argmax(test_overlap2)
                if test_overlap1[k] is not None:
                    self.H[i1:i2, j1:j2] += (s|sk[k])
                    self.H[i1:i2, j1:j2] += (s| ext_svt_G |sf[k])
                    self.S[i1:i2, j1:j2] += (s|sf[k])

        return self.H, self.S

