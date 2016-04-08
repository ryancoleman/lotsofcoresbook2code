import numpy as np
from ase.units import Hartree
from gpaw.transport.tools import aa1d, interpolate_array, \
                          collect_atomic_matrices, distribute_atomic_matrices
from gpaw.transport.io import Transport_IO
from gpaw.grid_descriptor import GridDescriptor
from gpaw import parsize_domain
from gpaw.domain import decompose_domain
''' 
     ---------------------------------------
      side    |                     |  side
          o  o| o  o  o  o  o  o  o |o  o
       -      |                     |     +
          o  o| o  o  o  o  o  o  o |o  o
              |                     |
       ---------------------------------------
         Left |                     |Right
         Lead | <-- Scattering -->  |Lead
              |       Region        |

 class Surrounding is used to deal with the projection close to the
 boundary and combine the potential or density information of the
 scattering region and leads.
'''

class Side:
    #Describe the electrode boundary
    def __init__(self, type, direction, kpt_comm, domain_comm, h):
        # Direction: '-' for the left electrode, '+' for the right electrode

        self.type = type
        self.direction = direction
        self.kpt_comm = kpt_comm
        self.domain_comm = domain_comm
        self.tio = Transport_IO(kpt_comm, domain_comm)
        self.h_cz = h

    def abstract_boundary(self):
        # Abtract the effective potential, hartree potential, and average density
        #out from the electrode calculation.
        map = {'-': '0', '+': '1'}
        data = self.tio.read_data(filename='Lead_' + 
                                  map[self.direction], option='Lead')
        nn = data['fine_N_c'][2]
        ns = data['nspins']
        N_c = data['N_c']
        cell_cv = data['cell_cv']
        pbc_c = data['pbc_c']
        #parsize_c = data['parsize_c']
        if type(parsize_domain) is int:
            parsize_c = None
            assert parsize_domain == self.domain_comm.size
        else:
            parsize_c = parsize_domain
        if parsize_c is None:
            parsize_c = decompose_domain(N_c, self.domain_comm.size)
        parsize_c = np.array(parsize_c)
        d1 = N_c[0] // 2
        d2 = N_c[1] // 2
       
        vHt_g = data['vHt_g']
        vt_sg = data['vt_sg']
        nt_sg = data['nt_sg']
        rhot_g = data['rhot_g']
        vt_sG = data['vt_sG']
        nt_sG = data['nt_sG']
        self.D_asp = data['D_asp']
        self.dH_asp = data['dH_asp']
        gd = GridDescriptor(N_c, cell_cv, pbc_c, self.domain_comm, parsize_c)
        finegd = gd.refine()

        self.boundary_vHt_g = None
        self.boundary_vHt_g1 = None
        self.boundary_vt_sg_line = None
        self.boundary_nt_sg = None
        self.boundary_rhot_g_line = None
        self.boundary_vt_sG = None
        self.boundary_nt_sG = None

        if self.tio.domain_comm.rank == 0:
            self.boundary_vHt_g = self.slice(nn, vHt_g)
            self.boundary_nt_sg = self.slice(nn, nt_sg)
            if self.direction == '-':
                other_direction= '+'
            else:
                other_direction= '-'
            h = self.h_cz / 2.0
            b_vHt_g0 = self.boundary_vHt_g.copy()
            b_vHt_g1 = self.boundary_vHt_g.copy()

            self.boundary_vHt_g = interpolate_array(b_vHt_g0,
                                                    finegd, h, 
                                                    self.direction)
            self.boundary_vHt_g1 = interpolate_array(b_vHt_g1,
                                                     finegd, h, 
                                                     other_direction)
            vt_sg = interpolate_array(vt_sg, finegd, 
                                      h, self.direction)
            self.boundary_vt_sg_line =  aa1d(vt_sg)
            self.boundary_nt_sg = interpolate_array(self.boundary_nt_sg,
                                                    finegd, h, self.direction)
            rhot_g = interpolate_array(rhot_g, finegd, h, self.direction)
            self.boundary_rhot_g_line = aa1d(rhot_g)

            nn /= 2
            h *= 2
            self.boundary_vt_sG = self.slice(nn, vt_sG)
            self.boundary_nt_sG = self.slice(nn, nt_sG)
            self.boundary_vt_sG = interpolate_array(self.boundary_vt_sG,
                                                    gd, h, self.direction)
            self.boundary_nt_sG = interpolate_array(self.boundary_nt_sG,
                                                    gd, h, self.direction)

    def slice(self, nn, in_array):
        if self.type == 'LR':
            seq1 = np.arange(nn)
            seq2 = np.arange(nn)
            di = len(in_array.shape) - 1
            if self.direction == '-':
                slice_array = np.take(in_array, seq1, axis=di)
            else:
                slice_array = np.take(in_array, seq2, axis=di)
        return slice_array


class Surrounding:
    #The potential and density enviroment of the scattering region
    def __init__(self, tp, type='LR'):
        self.type = type
        self.lead_num = tp.lead_num
        self.initialize(tp)

    def initialize(self, tp):
        if self.type == 'LR':
            self.sides = {}
            self.bias_index = {}
            self.side_basis_index = {}
            self.directions = ['-', '+']
            for i in range(self.lead_num):
                direction = self.directions[i]
                kpt_comm = tp.wfs.kd.comm
                gd_comm = tp.gd.comm
                side = Side('LR', direction, kpt_comm, 
                            gd_comm, tp.gd.h_cv[2,2])
                self.sides[direction] = side
                self.bias_index[direction] = tp.bias[i]
                self.side_basis_index[direction] = tp.lead_index[i]
            self.nn = tp.bnc[:]
            self.nn = np.array(self.nn)
            self.operator = \
                        tp.extended_calc.hamiltonian.poisson.operators[0]
        elif self.type == 'all':
            raise NotImplementError()
        self.calculate_sides()
        self.initialized = True

    def reset_bias(self, tp):
        self.bias = tp.bias
        for i in range(self.lead_num):
            direction = self.directions[i]
            self.bias_index[direction] = self.bias[i]
        self.combine(tp)

    def calculate_sides(self):
        if self.type == 'LR':
            for name, in self.sides:
                self.sides[name].abstract_boundary()
        if self.type == 'all':
            raise NotImplementError('type all not yet')

    def get_extra_density(self, tp, vHt_g):
        #help to solve poisson euqation for different left and right electrode
        if self.type == 'LR':
            rhot_g = tp.finegd1.zeros()
            self.operator.apply(vHt_g, rhot_g)
            nn = self.nn * 2
            self.extra_rhot_g = self.uncapsule(tp, nn, rhot_g, tp.finegd1,
                                               tp.finegd)

    def capsule(self, tp, nn, loc_in_array, in_cap_array, gd, gd0):
        ns = tp.nspins
        cap_array = gd.collect(in_cap_array)
        in_array = gd0.collect(loc_in_array)

        if gd.comm.rank == 0:
            if len(loc_in_array.shape) == 4:
                local_cap_array = gd.zeros(ns)
                cap_array[:, :, :, nn[0]:-nn[1]] = in_array
            else:
                local_cap_array = gd.zeros()
                cap_array[:, :, nn[0]:-nn[1]] = in_array
        else:
            if len(loc_in_array.shape) == 4:
                local_cap_array = gd.zeros(ns)
            else:
                local_cap_array = gd.zeros()
        gd.distribute(cap_array, local_cap_array)
        return local_cap_array

    def uncapsule(self, tp, nn, in_array, gd, gd0):
        nn1 = nn[0]
        nn2 = nn[1]
        ns = tp.nspins
        di = 2
        if len(in_array.shape) == 4:
            di += 1
            local_uncap_array = gd0.zeros(ns)
        else:
            local_uncap_array = gd0.zeros()
        global_in_array = gd.collect(in_array)
        if gd.comm.rank == 0:
            seq = np.arange(nn1, global_in_array.shape[di] - nn2)
            uncap_array = np.take(global_in_array, seq, axis=di)
        else:
            uncap_array = None
        gd0.distribute(uncap_array, local_uncap_array)
        return local_uncap_array

    def combine(self, tp):
        if self.type == 'LR':
            nn = self.nn * 2
            ham = tp.extended_calc.hamiltonian
            if ham.vt_sg is None:
                ham.vt_sg = ham.finegd.empty(ham.nspins)
                ham.vHt_g = ham.finegd.zeros()
                ham.vt_sG = ham.gd.zeros(ham.nspins)
                ham.poisson.initialize()
                if not tp.fixed:
                    tp.inner_poisson.initialize()

            bias_shift0 = self.bias_index['-'] / Hartree
            bias_shift1 = self.bias_index['+'] / Hartree
            b_vHt_g0 = self.sides['-'].boundary_vHt_g
            b_vHt_g1 = self.sides['+'].boundary_vHt_g
            b_vHt_g01 = self.sides['-'].boundary_vHt_g1
            b_vHt_g11 = self.sides['+'].boundary_vHt_g1
            if tp.fixed:
                if tp.gd.comm.rank == 0:
                    tp.inner_poisson.initialize(b_vHt_g0 + bias_shift0,
                                                    b_vHt_g1 + bias_shift1)
                else:
                    tp.inner_poisson.initialize(None, None)

            if tp.gd.comm.rank == 0:
                vHt_g = ham.finegd.zeros(global_array=True)
                extra_vHt_g = ham.finegd.zeros(global_array=True)
                nt_sg = ham.finegd.zeros(tp.nspins, global_array=True)
                nt_sG = ham.gd.zeros(tp.nspins, global_array=True)

                vHt_g[:, :, :nn[0]] = b_vHt_g0 + bias_shift0
                vHt_g[:, :, -nn[1]:] = b_vHt_g1 + bias_shift1

                mnn = np.min(nn)
                extra_vHt_g[:, :, nn[0]-mnn:nn[0]] = bias_shift0 + \
                                 b_vHt_g0[:,:,-mnn:] - b_vHt_g11[:,:,-mnn:]
                if nn[1] != mnn:
                    extra_vHt_g[:, :, -nn[1]:-nn[1]+mnn] = bias_shift1 + \
                                 b_vHt_g1[:,:,:mnn] - b_vHt_g01[:,:,:mnn]
                else:
                    extra_vHt_g[:, :, -nn[1]:] = bias_shift1 + \
                                 b_vHt_g1[:,:,:mnn] - b_vHt_g01[:,:,:mnn]

                nt_sg[:, :, :, :nn[0]] = self.sides['-'].boundary_nt_sg
                nt_sg[:, :, :, -nn[1]:] = self.sides['+'].boundary_nt_sg

                nn /= 2
                nt_sG[:, :, :, :nn[0]] = self.sides['-'].boundary_nt_sG
                nt_sG[:, :, :, -nn[1]:] = self.sides['+'].boundary_nt_sG
            else:
                nt_sG = None
                nt_sg = None
                vHt_g = None
                extra_vHt_g = None

            loc_extra_vHt_g = ham.finegd.zeros()
            self.nt_sg = ham.finegd.zeros(tp.nspins)
            self.nt_sG = ham.gd.zeros(tp.nspins)

            ham.gd.distribute(nt_sG, self.nt_sG)
            ham.finegd.distribute(nt_sg, self.nt_sg)
            ham.finegd.distribute(vHt_g, ham.vHt_g)
            ham.finegd.distribute(extra_vHt_g, loc_extra_vHt_g)
            self.get_extra_density(tp, loc_extra_vHt_g)

    def combine_vHt_g(self, tp, vHt_g):
        nn = self.nn * 2
        extended_vHt_g = tp.extended_calc.hamiltonian.vHt_g
        tp.extended_calc.hamiltonian.vHt_g = self.capsule(tp, nn, vHt_g,
                                                               extended_vHt_g,
                                                               tp.finegd1,
                                                               tp.finegd)

    def combine_nt_sG(self, tp, nt_sG):
        nn = self.nn
        self.nt_sG = self.capsule(tp, nn, nt_sG, self.nt_sG, tp.gd1,
                                  tp.gd)
        return self.nt_sG

    def combine_dH_asp(self, tp, dH_asp):
        ham = tp.extended_calc.hamiltonian
        all_dH_asp = dH_asp[:]
        for i in range(self.lead_num):
            direction = self.directions[i]
            side = self.sides[direction]
            for dH_ap in side.dH_asp:
                all_dH_asp.append(dH_ap)
        ham.dH_asp = {}
        for a, D_sp in tp.extended_D_asp.items():
            ham.dH_asp[a] = np.zeros_like(D_sp)
        distribute_atomic_matrices(all_dH_asp, ham.dH_asp, ham.setups)

    def refresh_vt_sG(self, tp):
        nn = self.nn
        gd = tp.extended_calc.gd
        bias_shift0 = self.bias_index['-'] / Hartree
        bias_shift1 = self.bias_index['+'] / Hartree
        vt_sG = gd.collect(tp.extended_calc.hamiltonian.vt_sG)
        if gd.comm.rank == 0:
            vt_sG[:, :, :, :nn[0]] = self.sides['-'].boundary_vt_sG + \
                                                                   bias_shift0
            vt_sG[:, :, :, -nn[1]:] = self.sides['+'].boundary_vt_sG + \
                                                                   bias_shift1
        gd.distribute(vt_sG, tp.extended_calc.hamiltonian.vt_sG)

