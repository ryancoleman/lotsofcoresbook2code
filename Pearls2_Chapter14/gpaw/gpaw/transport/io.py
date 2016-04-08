#This module is used to store and read some temperary data.
import cPickle
import os
from gpaw.mpi import world
from gpaw.transport.tools import collect_atomic_matrices, gather_ndarray_dict
import numpy as np

class Transport_IO:
    def __init__(self, kpt_comm, domain_comm):
        self.kpt_comm = kpt_comm
        self.domain_comm = domain_comm
        assert self.kpt_comm.size * self.domain_comm.size == world.size
        self.dir_name = 'temperary_data'
        if world.rank == 0:
            if not os.access(self.dir_name, os.F_OK):
                os.mkdir(self.dir_name)
        world.barrier()
        self.filenames = self.default_temperary_filenames()

    def default_temperary_filenames(self):
        filenames = {}
        for name in ['Lead']:
            filenames[name] = name
        for i in range(self.kpt_comm.size):
            for j in range(self.domain_comm.size):
                name = 'KC_' + str(i) + '_DC_' + str(j) +'AD'
                #Local analysis data on kpt_comm i and domain_comm j
                filenames[name] = name
        return filenames

    def read_data(self, bias_step=0, filename=None, option='Analysis'):
        if option == 'Lead':
            if filename is None:
                filename = self.filenames[option]
            fd = file(self.dir_name + '/' + filename, 'r')
            data = cPickle.load(fd)
            fd.close()
        elif option == 'Analysis':
            name = 'KC_' + str(self.kpt_comm.rank) + '_DC_' + \
                str(self.domain_comm.rank) +'AD'
            fd = file(self.dir_name + '/' + self.filenames[name] + '_bias_' + str(bias_step), 'r')
            data = cPickle.load(fd)
            fd.close()
        else:
            raise NotImplementError
        return data
           
    def save_data(self, obj, bias_step=0, filename=None, option='Analysis'):
        #       option ------map------ obj
        #       Analysis            Transport.Analysor
        #        Lead                Lead_Calc
        data = self.collect_data(obj, option)
        if option == 'Lead':
            if world.rank == 0:
                if filename is None:
                    filename = self.filenames[option]
                fd = file(self.dir_name + '/' + filename, 'wb')
                cPickle.dump(data, fd, 2)
                fd.close()
        elif option == 'Analysis':
            name = 'KC_' + str(self.kpt_comm.rank) + '_DC_' + \
                str(self.domain_comm.rank) +'AD'
            fd = file(self.dir_name + '/' + self.filenames[name] + '_bias_' + str(bias_step), 'wb')
            cPickle.dump(data, fd, 2)
            fd.close()
        else:
            raise NotImplementError()

    def collect_data(self, obj, option):
        if option == 'Lead':
            data = self.collect_lead_data(obj)
        elif option == 'Analysis':
            data = self.collect_analysis_data(obj)
        else:
            raise NotImplementError
        return data

    def collect_lead_data(self, obj):
        data = {}
        data['bzk_kc'] = obj.wfs.kd.bzk_kc
        #data['bzk_qc'] = obj.wfs.bzk_qc
        data['ibzk_kc'] = obj.wfs.kd.ibzk_kc
        data['ibzk_qc'] = obj.wfs.kd.ibzk_qc
        #data['setups'] = obj.wfs.setups
        data['cell_cv'] = obj.gd.cell_cv
        data['parsize_c'] = obj.gd.parsize_c
        data['pbc_c'] = obj.gd.pbc_c
        data['N_c'] = obj.gd.N_c
        data['nao'] = obj.wfs.setups.nao
        data['fine_N_c'] = obj.finegd.N_c
        data['nspins'] = obj.wfs.nspins
        data['fermi'] = obj.get_fermi_level()

        den, ham = obj.density, obj.hamiltonian
        gd = obj.gd
        vt_sG = gd.collect(ham.vt_sG)
        nt_sG = gd.collect(den.nt_sG)
        data['vt_sG'] = vt_sG
        data['nt_sG'] = nt_sG

        finegd = obj.finegd
        vt_sg = finegd.collect(ham.vt_sg)
        nt_sg = finegd.collect(den.nt_sg)
        vHt_g = finegd.collect(ham.vHt_g)
        rhot_g = finegd.collect(den.rhot_g)
        data['vt_sg'] = vt_sg
        data['nt_sg'] = nt_sg
        data['vHt_g'] = vHt_g
        data['rhot_g'] = rhot_g

        D_asp = collect_atomic_matrices(den.D_asp, den.setups,
                                        den.nspins, den.gd.comm,
                                        den.rank_a)
        dH_asp = collect_atomic_matrices(ham.dH_asp, ham.setups,
                                         ham.nspins, ham.gd.comm,
                                         ham.rank_a)
        data['D_asp'] = D_asp
        data['dH_asp'] = dH_asp
        return data

    def collect_analysis_data(self, obj):
        return obj.data
        
    def arrange_analysis_data(self, n_bias_step, n_ion_step, analysis_mode):           
        data = self.read_data(bias_step=n_bias_step, option='Analysis')
        parsize = data['domain_parsize']
        parpos = data['domain_parpos']
        domain_rank = data['domain_rank']
        kpt_rank = data['kpt_rank']
        kpt_size = data['kpt_size']
        transmission_dimension = data['transmission_dimension']
        dos_dimension = data['dos_dimension']
        assert kpt_rank == self.kpt_comm.rank
        assert domain_rank == self.domain_comm.rank

        #collect transmission, dos, current
        global_data = {}
        flag = 'K_' + str(kpt_rank) + 'D_' + str(domain_rank) + '_'
        global_data[flag + 'parpos'] = data['domain_parpos']
        global_data[flag + 'tc'] = data['tc']
        global_data[flag + 'dos'] = data['dos']

        global_data = gather_ndarray_dict(global_data, self.kpt_comm)
        
        global_data[flag + 'vt'] = data['vt']
        global_data[flag + 'vtx'] = data['vtx']
        global_data[flag + 'vty'] = data['vty']
        global_data[flag + 'nt'] = data['nt']
        global_data[flag + 'ntx'] = data['ntx']
        global_data[flag + 'nty'] = data['nty']

        global_data = gather_ndarray_dict(global_data, self.domain_comm)

        global_data['lead_fermi'] = data['lead_fermi']
        global_data['bias'] = data['bias']
        global_data['gate'] = data['gate']
        global_data['charge'] = data['charge']
        global_data['magmom'] = data['magmom']
        global_data['local_magmom'] = data['local_magmom']

        global_data['newform'] = True
        global_data['domain_parsize'] = data['domain_parsize']
        global_data['kpt_size'] = data['kpt_size']
        global_data['transmission_dimension'] = data['transmission_dimension']
        global_data['dos_dimension'] = data['dos_dimension']

        world.barrier()
        if world.rank == 0:
            if analysis_mode:
                filename = '/abias_step_' + str(n_bias_step)
            else:
                filename = '/bias_step_' + str(n_bias_step)
            fd = file('analysis_data/ionic_step_' + str(n_ion_step)
                      + filename, 'wb')
            cPickle.dump(global_data, fd, 2)
            fd.close()
            for root, dirs, files in os.walk('temperary_data'):
                for name in files:
                    if 'AD' in name:
                        os.remove(os.path.join(root, name))
