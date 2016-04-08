from ase import Atoms, Atom
from ase.parallel import parprint, paropen
from ase.units import Hartree, Bohr
from gpaw import GPAW, debug, dry_run, PoissonSolver
from gpaw.mixer import Mixer, MixerSum, MixerDif, BroydenMixer, BroydenMixerSum
from gpaw.poisson import FixedBoundaryPoissonSolver

from gpaw import restart as restart_gpaw
from gpaw.grid_descriptor import GridDescriptor
from gpaw.transformers import Transformer
from gpaw.mpi import world
from gpaw import mpi
from gpaw.utilities.memory import maxrss

from gpaw.transport.tools import tri2full, dot, \
          get_atom_indices, substract_pk, get_lcao_density_matrix, \
          get_pk_hsd, get_matrix_index, aa1d, aa2d, collect_atomic_matrices,\
          distribute_atomic_matrices, fermidistribution

from gpaw.transport.sparse_matrix import Tp_Sparse_HSD, Banded_Sparse_HSD, \
                                  CP_Sparse_HSD, Se_Sparse_Matrix

from gpaw.transport.contour import Contour, EnergyNode
from gpaw.transport.surrounding import Surrounding
from gpaw.transport.selfenergy import LeadSelfEnergy
from gpaw.transport.analysor import Transport_Analysor, Transport_Plotter
from gpaw.transport.io import Transport_IO
from gpaw.utilities import pack2,unpack,unpack2

import gpaw
import numpy as np
import cPickle
import time

class DBG(object):
    """A simple debugging and logging class"""
    def __init__(self, lfiled=None, enable=False):
        super(DBG, self).__init__()
        self.enable = enable
        if lfiled is None:
            self.logfile = paropen('transport_dbg.log', 'w')
            #self.logfile = open('transport_dbg.log', 'a')
        else:
            self.logfile = lfiled
        self.start_time = time.time()
    
    def log(self, dbgstr):
        """Write a debug string to the given logfile"""
        if not self.enable:
            return
        t_now = time.time() - self.start_time
        #self.logfile.write('N{0}> {1} : {2}\n'\
        #    .format(mpi.rank, t_now, dbgstr))
        self.logfile.write('N0> {0} : {1}\n'\
            .format(t_now, dbgstr))
        self.logfile.flush()
        

class FixedBC_GridDescriptor(GridDescriptor):
    def get_boxes(self, spos_c, rcut, cut=True):
        # The effect of this grid descriptor class is to always use cut=True
        return GridDescriptor.get_boxes(self, spos_c, rcut, True)


class Lead_Calc(GPAW):
    def dry_run(self):
        pass

    def initialize(self, *args, **kwargs):
        GPAW.initialize(self, *args, **kwargs)
        self.gd = self.wfs.gd
        self.finegd = self.density.finegd
        
class Transport(GPAW):
    def __init__(self, **transport_kwargs):
        self.dbg = DBG()
        self.log = self.dbg.log
        self.grid_descriptor_class = FixedBC_GridDescriptor
        self.set_transport_kwargs(**transport_kwargs)
        GPAW.__init__(self, **self.gpw_kwargs)
            
    def initialize(self, *args, **kwargs):
        GPAW.initialize(self, *args, **kwargs)
        self.gd = self.wfs.gd
        self.finegd = self.density.finegd
            
    def set_transport_kwargs(self, **transport_kwargs):
        """ Illustration of keywords::
          o  o  o  o  o  o  o  o  o  o  o  o  o
          0  1  2  3  4  5  6  7  8  9  10 11 12
          .            | mol_atoms |
          | pl_atoms1 |              | pl_atoms2 |
        """
        kw = transport_kwargs  
        p =  self.set_default_transport_parameters()
        self.gpw_kwargs = kw.copy()
        for key in kw:
            if key in [ 'use_lead',
                        'identical_leads',
                        'pl_atoms',
                        'pl_cells',
                        'pl_kpts',
                        'leads',
                        'use_buffer',
                        'buffer_atoms',
                        'edge_atoms',
                        'bias',
                        'lead_restart',
                        'special_datas',
                        'neutral_steps',
                        'plot_eta',
                        'plot_energy_range',
                        'plot_energy_point_num',
                        'vaccs',
                        'lead_guess',
                        'neutral',
                        'buffer_guess',
                        'fill_guess',
                        'cell_atoms',
                        'sf',
                        'lead_sf',
                        'eta',
                        'lead_atoms',
                        'nleadlayers',
                        'mol_atoms',
                        'la_index',
                        'total_charge',
                        'alpha',
                        'beta_guess',
                        'theta',
                        'LR_leads',
                        'gate',
                        'gate_mode',
                        'gate_atoms',
                        'gate_fun',
                        'recal_path',
                        'min_energy',
                        'fix_contour',
                        'hubbard_parameters',
                        'use_qzk_boundary',
                        'n_bias_step',
                        'n_ion_step',
                        'scat_restart',
                        'save_file',
                        'restart_file',
                        'save_lead_hamiltonian',
                        'non_sc',
                        'fixed_boundary',
                        'guess_steps',
                        'foot_print',
                        'data_file', 
                        'extended_atoms',
                        'restart_lead_hamiltonian',
                        'analysis_data_list',
                        'save_bias_data',
                        'perturbation_steps',
                        'perturbation_magmom',
                        'perturbation_charge',
                        'perturbation_atoms',
                        'analysis_mode',
                        'normalize_density',
                        'se_data_path',                     
                        'neintmethod',
                        'neintstep',
                        'eqinttol',
                        'extra_density',
                        'lead_calculators',
                        'enable_dbglog']:
                
                del self.gpw_kwargs[key]
                p[key] = kw[key]
            if key in ['spinpol']:
                p['spinpol'] = kw['spinpol']
            if key in ['verbose']:
                p['verbose'] = kw['verbose']
        
        self.transport_parameters = p
        self.use_lead = p['use_lead']
        self.identical_leads = p['identical_leads']
        self.leads = p['leads']
        self.pl_atoms = p['pl_atoms']
        self.lead_num = len(self.pl_atoms)
        self.bias = p['bias']
        self.lead_calculators = p['lead_calculators']
        self.enable_dbglog = p['enable_dbglog']
        self.dbg.enable =  self.enable_dbglog
        #parprint(">>> lead_calculators: {0}".format(self.lead_calculators))
        #parprint(">>> gpw_kwargs:")
        #for key in self.gpw_kwargs:
        #    parprint(">>> key {0} value {1}".format(key, self.gpw_kwargs[key]))
        
        if self.use_lead:
            self.pl_cells = p['pl_cells']
            self.pl_kpts = p['pl_kpts']
            self.lead_restart = p['lead_restart']
            self.use_buffer = p['use_buffer']
            self.buffer_atoms = p['buffer_atoms']
            self.edge_atoms = p['edge_atoms']

            assert self.lead_num == len(self.pl_cells)
            #assert self.lead_num == len(self.buffer_atoms)
            #assert self.lead_num == len(self.edge_atoms[0])
            #assert self.lead_num == len(self.bias)
            
            self.lead_atoms = p['lead_atoms']
            if self.lead_atoms is None:
                self.lead_atoms = self.pl_atoms
            self.nleadlayers = p['nleadlayers']
            self.mol_atoms = p['mol_atoms']
            self.la_index = p['la_index']
            
        self.LR_leads = p['LR_leads']
        self.gate = p['gate']
        self.gate_mode = p['gate_mode']
        self.gate_fun = p['gate_fun']
        self.gate_atoms = p['gate_atoms']
        self.recal_path = p['recal_path']
        self.plot_eta = p['plot_eta']
        self.eta = p['eta']
        self.plot_energy_range = p['plot_energy_range']
        self.plot_energy_point_num = p['plot_energy_point_num']
        self.alpha = p['alpha']
        self.beta_guess = p['beta_guess']
        self.theta = p['theta']
        self.min_energy = p['min_energy']
        self.use_qzk_boundary = p['use_qzk_boundary']
        self.scat_restart = p['scat_restart']
        self.guess_steps = p['guess_steps']
        self.foot_print = p['foot_print']
        self.save_file = p['save_file']
        self.restart_file = p['restart_file']
        self.restart_lead_hamiltonian = p['restart_lead_hamiltonian']
        self.neintmethod = p['neintmethod']
        self.neintstep = p['neintstep']
        self.n_bias_step = p['n_bias_step']
        self.n_ion_step = p['n_ion_step']
        self.fixed = p['fixed_boundary']
        self.fix_contour = p['fix_contour']
        self.vaccs = p['vaccs']
        self.lead_guess = p['lead_guess']
        self.fill_guess = p['fill_guess']
        self.cell_atoms = p['cell_atoms']
        self.sf = p['sf']
        self.lead_sf = p['lead_sf']
        self.buffer_guess = p['buffer_guess']
        self.neutral = p['neutral']
        self.neutral_steps = p['neutral_steps']
        self.total_charge = p['total_charge']
        self.non_sc = p['non_sc']
        self.data_file = p['data_file']
        self.extended_atoms = p['extended_atoms']
        self.special_datas = p['special_datas']
        self.analysis_data_list = p['analysis_data_list']
        self.save_bias_data = p['save_bias_data']
        self.save_lead_hamiltonian = p['save_lead_hamiltonian']
        self.se_data_path = p['se_data_path']
        self.analysis_mode = p['analysis_mode']
        self.normalize_density = p['normalize_density']
        self.extra_density = p['extra_density']
        self.hubbard_parameters = p['hubbard_parameters']
        self.eqinttol = p['eqinttol']
        self.spinpol = p['spinpol']
        self.perturbation_steps = p['perturbation_steps']
        self.perturbation_charge = p['perturbation_charge']
        self.perturbation_magmom = p['perturbation_magmom']
        self.perturbation_atoms = p['perturbation_atoms']
        
        self.verbose = p['verbose']
        self.d = p['d']
       
        if self.scat_restart and self.restart_file is None:
            self.restart_file = 'bias_data1'
        
        self.master = (world.rank==0)
        
        bias = self.bias
        
        if self.LR_leads and self.lead_num != 2:
            raise RuntimeError('wrong way to use keyword LR_leads')
        
        self.initialized_transport = False
        self.analysis_parameters = []
        self.optimize = False
        kpts = kw['kpts']
        if np.product(kpts) == kpts[self.d]:
            self.gpw_kwargs['symmetry'] = 'off'
        else:
            self.gpw_kwargs['symmetry'] = {'point_group': False}
        self.scat_ntk = 1
        if kpts[2] != 1:
            if self.non_sc:
                self.scat_ntk = kpts[2]
            else:
                self.scat_ntk = 1
        self.gpw_kwargs['kpts'] = kpts[:2] + (1,)
        # ! THa: Hack: 
        # ! Also save the the parameters
        # !'plot_energy_range' and for later analysis 'plot_energy_point_num'
        if self.master:
            ee = \
             np.linspace(self.transport_parameters['plot_energy_range'][0], \
             self.transport_parameters['plot_energy_range'][1], \
             num=self.transport_parameters['plot_energy_point_num'])
            np.save('plot_energy_range.npy', ee)
        

    def set_analysis_parameters(self, **analysis_kwargs):
        self.analysis_parameters = analysis_kwargs
        for key in self.analysis_parameters:
            if key not in ['energies', 'lead_pairs', 'dos_project_atoms',
                       'project_molecular_levels', 'isolate_atoms', 'project_equal_atoms',
                        'dos_project_orbital',
                        'trans_project_orbital', 'eig_trans_channel_energies',
                        'eig_trans_channel_num', 'dos_realspace_energies']:
                raise ValueError('no keyword %s for analysis' % key)    

    def set_default_transport_parameters(self):
        self.log('set_default_transport_parameters()')
        p = {}
        p['use_lead'] = True
        p['identical_leads'] = False
        p['pl_atoms'] = []
        p['pl_cells'] = []
        p['pl_kpts'] = []
        p['use_buffer'] = False
        p['buffer_atoms'] = None
        p['edge_atoms'] = None
        p['mol_atoms'] = None
        p['leads'] = None
        p['bias'] = [0, 0]
        p['lead_calculators'] = None
        p['d'] = 2
        p['lead_restart'] = False
        p['restart_lead_hamiltonian'] = False
        p['lead_atoms'] = None
        p['extended_atoms'] = None
        p['nleadlayers'] = [1, 1]
        p['la_index'] = None
        p['data_file'] = None
        p['analysis_data_list'] = ['tc', 'force']
        p['special_datas'] = []
        p['save_bias_data'] = True
        p['analysis_mode'] = False
        p['normalize_density'] = True
        p['extra_density'] = False
        p['se_data_path'] = None
        p['neintmethod'] = 0
        p['neintstep'] = 0.02
        p['n_bias_step'] = 0
        p['n_ion_step'] = 0
        p['eqinttol'] = 1e-4
        p['plot_eta'] = 0.0001
        p['eta'] = 0.01
        p['plot_energy_range'] = [-5.,5.]
        p['plot_energy_point_num'] = 201
        p['alpha'] = 0.0
        p['beta_guess'] = 0.1
        p['theta'] = 1.
        p['vaccs'] = None
        p['LR_leads'] = True
        p['lead_guess'] = False
        p['buffer_guess'] = False
        p['fill_guess'] = False
        p['cell_atoms'] = None
        p['sf'] = None
        p['lead_sf'] = None
        p['neutral'] = True
        p['neutral_steps'] = None
        p['total_charge'] = 0
        p['perturbation_steps'] = []
        p['perturbation_magmom'] = None
        p['perturbation_charge'] = 0
        p['perturbation_atoms'] = []
        p['hubbard_parameters'] = None
        
        p['gate'] = 0
        p['gate_mode'] = 'VG'
        p['gate_fun'] = None
        p['gate_atoms'] = None
        p['recal_path'] = False
        p['min_energy'] = -700
        p['guess_steps'] = 30
        p['foot_print'] = True
        p['use_qzk_boundary'] = False
        p['scat_restart'] = False
        p['save_file'] = False
        p['save_lead_hamiltonian'] =False
        p['restart_file'] = None
        p['fixed_boundary'] = True
        p['fix_contour'] = False
        p['non_sc'] = False
        p['spinpol'] = False
        p['verbose'] = False
        p['enable_dbglog'] = False
        return p     

    def set_atoms(self, atoms):
        self.log('set_atoms()')
        self.adjust_atom_positions(atoms)
        self.atoms = atoms.copy()
        if self.edge_atoms is None:
            self.edge_atoms = [[0, len(self.pl_atoms[1]) - 1],
                                [0, len(self.atoms) -1]]
        if self.mol_atoms is None:
            self.mol_atoms = range(len(self.atoms))
            for i in range(self.lead_num):
                for ind in self.pl_atoms[i]:
                    self.mol_atoms.remove(ind)            

    def adjust_atom_positions(self, atoms):
        self.log('adjust_atom_positions')
        # match the scattering region and the lead region
        # to get a correct boundary which is used to
        # solve the Poisson equation
        if self.identical_leads or self.vaccs is None:
            atoms.center()
        else:
            atoms.center(axis=0)
            atoms.center(axis=1)
            lb = np.min(atoms.positions[:,2])
            rb = np.max(atoms.positions[:,2])
            dis = self.vaccs[0] - lb
            atoms.positions[:,2] += dis
            assert abs(np.diag(atoms.cell)[2] - rb - dis -self.vaccs[1]) < 0.005

    def init_leads_from_gpw(self):
        """
        Initialize lead calculators from given, well converged gpaw alculator
        objects stored in self.lead_calculators.
        
        Try to mimic behaviour of calculate_leads()
        
        """
        self.log("init_leads_from_gpw()")
        self.setups_lead = []
        self.lead_fermi = []
        ##p = self.gpw_kwargs.copy()
        ##p['kpts'] = self.pl_kpts
        ##p['parallel'] = self.input_parameters['parallel'].copy()
        ##p['parallel'].update(band=self.wfs.bd.comm.size,
        ##                      # kpt=self.wfs.kd.comm.size
        ##                      domain=self.wfs.gd.parsize_c)

        for i in range(self.lead_num):
            if not (self.identical_leads and i > 0):
                atoms = self.get_lead_atoms(i, init_calc=False)
                # get all the keywords from the saved lead
                # calculator except the 'parallel' options
                self.log("  check")
                _tmp = GPAW(self.lead_calculators[i])
                pl_params = _tmp.input_parameters.copy()
                pl_params['parallel'] = self.input_parameters['parallel'].copy()
                self.log('  self.input_parameters[parallel] {0}'.format(pl_params['parallel']))
                #pl_params['parallel'].update(band=self.wfs.bd.comm.size, domain=self.wfs.gd.parsize_c)
                #parprint('>>> self.input_parameters[parallel] updated {0}'.format(pl_params['parallel']))
                _lc = Lead_Calc(self.lead_calculators[i], **pl_params)
                #atoms = _lc.get_atoms()
                _lc.initialize(atoms)
                # initialize *.gpw
                self.log("  calling _lc.set_positions(atoms)")
                _lc.set_positions(atoms)
                _lc.wfs.eigensolver.iterate(_lc.hamiltonian, _lc.wfs)
                _lc.scf.converged = True
                atoms.set_calculator(_lc)
                atoms.get_potential_energy()
                
                fermi = atoms.calc.get_fermi_level()
                if i > 0:
                    shift = self.lead_fermi[0] - fermi
                    fermi += shift
                    atoms.calc.hamiltonian.vt_sG += shift / Hartree
                self.lead_fermi.append(fermi)    
                self.setups_lead.append(atoms.calc.wfs.setups)
                self.collect_leads_matrices(atoms.calc, i)
            self.tio.save_data(atoms.calc, option='Lead',
                               filename='Lead_'+ str(i))

    def calculate_leads(self):
        self.log('calculate_leads()')
        self.setups_lead = []
        self.lead_fermi = []
        for i in range(self.lead_num):
            if not (self.identical_leads and i > 0):
                atoms = self.get_lead_atoms(i)
                atoms.get_potential_energy()
                fermi = atoms.calc.get_fermi_level()
                if i > 0:
                    shift = self.lead_fermi[0] - fermi
                    fermi += shift
                    atoms.calc.hamiltonian.vt_sG += shift / Hartree
                self.lead_fermi.append(fermi)    
                self.setups_lead.append(atoms.calc.wfs.setups)
                self.collect_leads_matrices(atoms.calc, i)
            self.tio.save_data(atoms.calc, option='Lead',
                               filename='Lead_'+ str(i))

    def restart_leads(self):
        self.log('restart_leads()')
        self.setups_lead = []
        for i in range(self.lead_num):
            atoms = self.get_lead_atoms(i)
            calc = atoms.calc
            calc.initialize(atoms)
            self.setups_lead.append(atoms.calc.wfs.setups)
            data = self.tio.read_data(filename='Lead_' + str(i),
                                                  option='Lead')
            vt_sG = data['vt_sG']
            dH_asp = data['dH_asp']
            calc.hamiltonian.dH_asp = {}
            for i, dH_p in enumerate(dH_asp):
                calc.hamiltonian.dH_asp[i] = dH_p
            #calc.hamiltonian.vt_sG = calc.gd.empty()
            calc.gd.distribute(vt_sG, calc.hamiltonian.vt_sG) 
            self.collect_leads_matrices(calc, i)
            
    def parameters_test(self):
        self.log('parameters_test()')
        for i in range(self.lead_num):
            atoms = self.get_lead_atoms(i)
            atoms.calc.initialize(atoms)
        self.initialize()    

    def define_leads_related_variables(self):
        self.log('define_leads_related_variables()')
        self.nblead = []  # numer of basis in lead
        self.bnc = []    # boundary N_c, gd.N_c[self.d] in lead
        self.bzk_kc_lead = [] 
        self.ibzk_kc_lead = []
        self.ibzk_qc_lead = []
        self.edge_index = [[None] * self.lead_num, 
                           [None] * self.lead_num] #index of edge atom
        world.barrier()
        for i in range(self.lead_num):
            data = self.tio.read_data(filename='Lead_' + str(i),
                                                  option='Lead')
            self.nblead.append(data['nao'])
            self.bnc.append(data['N_c'][self.d])
            self.bzk_kc_lead.append(data['bzk_kc'])
            self.ibzk_kc_lead.append(data['ibzk_kc'])
            self.ibzk_qc_lead.append(data['ibzk_qc'])
        
    def initialize_transport(self):
        self.log('initialize_transport()')
        # calculate the lead and generate a guess hamiltonian for
        # scattering region
        if dry_run:
            self.parameters_test()
        self.cvg_ham_steps = 0
        self.initialize()
        self.nspins = self.wfs.nspins
        self.npk = len(self.wfs.kd.ibzk_kc)
        self.my_npk = len(self.wfs.kd.ibzk_qc)
        self.my_nspins = len(self.wfs.kpt_u) // self.my_npk
        self.ntklead = self.pl_kpts[self.d]
        self.initialize_lead_matrix()

        self.tio = Transport_IO(self.wfs.kd.comm, self.gd.comm)
        if not self.restart_lead_hamiltonian:
            if self.lead_calculators is None:
                self.calculate_leads()
            else:
                self.init_leads_from_gpw()
        else:
            self.restart_leads()
        self.define_leads_related_variables()
        self.get_basis_indices()
        self.initialize_matrix()

        self.get_extended_atoms()
        calc = self.extended_atoms.calc
        calc.initialize(self.extended_atoms)
        if not self.use_qzk_boundary:
            del calc.density
        self.extended_calc = calc
        self.gd1, self.finegd1 = calc.gd, calc.finegd
  
        bzk_kc = self.wfs.kd.bzk_kc 
        self.gamma = len(bzk_kc) == 1 and not bzk_kc[0].any()
        self.nbmol = self.wfs.setups.nao

        self.get_ks_map()
        self.set_local_spin_index(self.wfs)
        self.set_local_spin_index(self.extended_calc.wfs)

        if self.npk == 1:
            self.lead_kpts = self.bzk_kc_lead[0]
        else:
            self.lead_kpts = self.ibzk_kc_lead[0]                
          
        if self.nbmol <= np.sum(self.nblead):
            self.use_buffer = False
            if self.master:
                self.text('Moleucle is too small, force not to use buffer')
           
        if self.use_buffer: 
            self.buffer = [len(self.buffer_index[i])
                                               for i in range(self.lead_num)]
            self.print_index = self.buffer_index
        else:
            self.buffer = [0] * self.lead_num
            self.print_index = self.lead_index
            
        self.set_buffer()
            
        self.current = 0
        self.linear_mm = None

        self.fermi = self.lead_fermi[0]
        #self.leads_fermi_lineup()
        world.barrier()
        
        self.timer.start('init surround')
        self.surround = Surrounding(self)  
        self.timer.stop('init surround')

        self.get_inner_setups()
        self.extended_D_asp = None        

        if not self.non_sc:
            self.timer.start('surround set_position')
            if not self.fixed:
                self.inner_poisson = PoissonSolver(nn=self.hamiltonian.poisson.nn)
            else:
                self.inner_poisson = FixedBoundaryPoissonSolver(nn=1)
            self.inner_poisson.set_grid_descriptor(self.finegd)
            self.interpolator = Transformer(self.gd1, self.finegd1,
                                            self.input_parameters.stencils[1])

            self.surround.combine(self)

            if self.use_qzk_boundary:
                self.extended_calc.set_positions()
            else:
                self.set_extended_positions()
            self.timer.stop('surround set_position')
       
        if not self.analysis_mode:
            if self.scat_restart:
                self.get_hamiltonian_initial_guess3()
            elif self.buffer_guess:
                self.get_hamiltonian_initial_guess2()
            elif self.fill_guess:
                self.get_hamiltonian_initial_guess4()
            else:
                self.get_hamiltonian_initial_guess()                
       
        self.initialize_gate()
        self.initialized_transport = True
        self.matrix_mode = 'sparse'
        if not hasattr(self, 'plot_option'):
            self.plot_option = None
        if np.abs(self.bias[1] - self.bias[0]) < 0.001:
            self.ground = True
        else:
            self.ground = False
        self.F_av = None

    def set_energies(self, energies):
        self.log('set_energies()')
        p = {}
        p['energies'] = energies
        self.set_analysis_parameters(**p)

    #def leads_fermi_lineup(self):
    #    for i in range(1, self.lead_num):
    #        shift = self.lead_fermi[0] - self.lead_fermi[i]
    #        self.lead_fermi[i] = self.lead_fermi[0]
    #        self.atoms_l[i].calc.hamiltonian.vt_sG += shift / Hartree
    #        self.atoms_l[i].calc.hamiltonian.vHt_g += shift / Hartree
    #        for pk in range(self.my_npk):
    #            for s in range(self.my_nspins):
    #                h_mm = self.lead_hsd[i].H[s][pk].recover()
    #                h_cmm = self.lead_couple_hsd[i].H[s][pk].recover()
    #                s_mm = self.lead_hsd[i].S[pk].recover()
    #                s_cmm = self.lead_couple_hsd[i].S[pk].recover()
    #                h_mm += shift * s_mm
    #                h_cmm += shift * s_cmm
    #                self.lead_hsd[i].reset(s, pk, h_mm, 'H')     
    #                self.lead_couple_hsd[i].reset(s, pk, h_cmm,'H')     
       
    def get_ks_map(self):
        self.log('get_ks_map()')
        # ks_map: s, k, rank
        # my_ks_map: s, q, rank (q is the k index in local processor)
        self.ks_map = np.zeros([self.npk * self.nspins, 3], int)
        self.my_ks_map = np.zeros([self.my_npk * self.my_nspins, 3], int)
        for i, kpt in enumerate(self.wfs.kpt_u):
            base = self.wfs.kd.comm.rank * self.my_npk * self.my_nspins
            self.ks_map[i + base, 0] = kpt.s
            self.ks_map[i + base, 1] = kpt.k
            self.ks_map[i + base, 2] = self.wfs.kd.comm.rank
            self.my_ks_map[i, 0] = kpt.s
            self.my_ks_map[i, 1] = kpt.k
            self.my_ks_map[i, 2] = self.wfs.kd.comm.rank            
        self.wfs.kd.comm.sum(self.ks_map)
    
    def set_local_spin_index(self, wfs):
        self.log('set_local_spin_index()')
        for kpt in wfs.kpt_u:
            if kpt.s == 0:
                kpt.v = 0
            elif self.my_nspins == 2:
                kpt.v = kpt.s
            else:
                kpt.v = 0

    def initialize_gate(self):
        self.log('initialize_gate()')
        if self.gate_mode == 'SN':
            assert self.gate_atoms is not None

        elif self.gate_mode == 'VM':
            setups = self.wfs.setups
            self.gate_basis_index = get_atom_indices(self.gate_atoms, setups)

        elif self.gate_mode == 'AN':
            setups = self.wfs.setups
            self.gate_basis_index = get_atom_indices(self.gate_atoms, setups)

    def get_hamiltonian_initial_guess2(self):
        self.log('get_hamiltonian_initial_guess2()')
        # get a hamiltonian guess for scattering region using buffer layer
        atoms = self.atoms.copy()
        cell = np.diag(atoms.cell)
        cell[2] += 15.0
        atoms.set_cell(cell)
        atoms.center()
        #atoms.pbc[self.d] = True
        kwargs = self.gpw_kwargs.copy()
        kwargs['poissonsolver'] = PoissonSolver(nn=2)
        kpts = kwargs['kpts']
        kpts = kpts[:2] + (1,)
        kwargs['kpts'] = kpts
        if hasattr(self.density.mixer, 'mixers'):
            kwargs['mixer'] = Mixer(self.beta_guess, 5, weight=100.0)
        elif hasattr(self.density.mixer, 'mixer'):
            kwargs['mixer'] = MixerDif(self.beta_guess, 5, weight=100.0)
        elif self.spinpol and hasattr(self.density.mixer, 'step'):
            kwargs['mixer'] = BroydenMixerSum(self.beta_guess, 5)
        elif self.spinpol and not hasattr(self.density.mixer, 'step'):
            kwargs['mixer'] = MixerSum(self.beta_guess, 5, weight=100.0)
        else:
            kwargs['mixer'] = BroydenMixer(self.beta_guess, 5)
        if 'txt' in kwargs and kwargs['txt'] != '-':
            kwargs['txt'] = 'guess_' + kwargs['txt']            
        atoms.set_calculator(gpaw.GPAW(**kwargs))
        calc = atoms.calc
        calc.initialize(atoms)
        calc.set_positions(atoms)
        
        wfs = calc.wfs
        hamiltonian = calc.hamiltonian
        occupations = calc.occupations
        density = calc.density
        scf = calc.scf
        
        for iter in range(self.guess_steps):
            wfs.eigensolver.iterate(hamiltonian, wfs)
            occupations.calculate(wfs)
            energy = hamiltonian.get_energy(occupations)
            scf.energies.append(energy)
            scf.check_convergence(density, wfs.eigensolver)
            density.update(wfs)
            hamiltonian.update(density)
            calc.print_iteration(iter)
        self.copy_mixer_history(calc, 'buffer') 
        self.initialize_hamiltonian_matrix(calc)      
        del calc
        self.boundary_align_up()        
            
    def get_hamiltonian_initial_guess(self):
        self.log('get_hamiltonian_initial_guess()')
        # get hamiltonian guess for scattering region using normal DFT
        atoms = self.atoms.copy()
        #atoms.pbc[self.d] = True
        kwargs = self.gpw_kwargs.copy()
        kwargs['poissonsolver'] = PoissonSolver(nn=2)
        kpts = kwargs['kpts']
        #kpts = kpts[:2] + (1,)
        #kwargs['kpts'] = kpts
        if self.non_sc:
            kwargs['kpts'] = kpts[:2] + (self.scat_ntk,)

        if hasattr(self.density.mixer, 'mixers'):
            kwargs['mixer'] = Mixer(self.beta_guess, 5, weight=100.0)
        elif hasattr(self.density.mixer, 'mixer'):
            kwargs['mixer'] = MixerDif(self.beta_guess, 5, weight=100.0)
        elif self.spinpol and hasattr(self.density.mixer, 'step'):
            kwargs['mixer'] = BroydenMixerSum(self.beta_guess, 5)
        elif self.spinpol and not hasattr(self.density.mixer, 'step'):
            kwargs['mixer'] = MixerSum(self.beta_guess, 5, weight=100.0)
        else:
            kwargs['mixer'] = BroydenMixer(self.beta_guess, 5)

        if 'txt' in kwargs and kwargs['txt'] != '-':
            kwargs['txt'] = 'guess_' + kwargs['txt']            
        atoms.set_calculator(gpaw.GPAW(**kwargs))
        calc = atoms.calc
        calc.initialize(atoms)
        calc.set_positions(atoms)
        
        wfs = calc.wfs
        hamiltonian = calc.hamiltonian
        occupations = calc.occupations
        density = calc.density
        scf = calc.scf
        
        if self.non_sc:
            if not self.scat_restart:
                self.guess_total_energy = atoms.get_potential_energy()
                calc = atoms.calc
                if self.save_file:
                    atoms.calc.write('scat.gpw')
                self.hamiltonian = atoms.calc.hamiltonian
                self.density = atoms.calc.density
                self.extended_calc.hamiltonian = self.hamiltonian
                self.extended_calc.gd = self.gd
                self.extended_calc.finegd = self.finegd
                self.extended_calc.wfs.basis_functions = self.wfs.basis_functions
   
            else:
                calc = self
                #self.recover_kpts(atoms.calc)                
                self.extended_calc.hamiltonian = self.hamiltonian
                self.extended_calc.gd = self.gd
                self.extended_calc.finegd = self.finegd
                self.extended_calc.wfs.basis_functions = self.wfs.basis_functions                
        else:        
            for iter in range(self.guess_steps):
                wfs.eigensolver.iterate(hamiltonian, wfs)
                occupations.calculate(wfs)
                energy = hamiltonian.get_energy(occupations)
                scf.energies.append(energy)
                scf.check_convergence(density, wfs.eigensolver)
                density.update(wfs)
                if self.extra_density:
                    density.rhot_g += self.surround.extra_rhot_g
                hamiltonian.update(density)
                calc.print_iteration(iter)
            self.copy_mixer_history(calc)       
        self.initialize_hamiltonian_matrix(calc)      
        if not (self.non_sc and self.scat_restart):
            del calc
        #atoms.get_potential_energy()

    def get_hamiltonian_initial_guess3(self):
        self.log('get_hamiltonian_initial_guess3()')
        #get hamiltonian_guess from a hamiltonian file
        fd = file(self.restart_file, 'r')
        self.bias, vt_sG, dH_asp = cPickle.load(fd)
        fd.close()
        self.surround.combine_dH_asp(self, dH_asp)
        self.gd1.distribute(vt_sG, self.extended_calc.hamiltonian.vt_sG) 
        h_spkmm, s_pkmm = self.get_hs(self.extended_calc)
        nb = s_pkmm.shape[-1]
        dtype = s_pkmm.dtype
        for q in range(self.my_npk):
            self.hsd.reset(0, q, s_pkmm[q], 'S', True)                
            for s in range(self.my_nspins):
                self.hsd.reset(s, q, h_spkmm[s, q], 'H', True)
                self.hsd.reset(s, q, np.zeros([nb, nb], dtype), 'D', True)

    def get_hamiltonian_initial_guess4(self):
        self.log('get_hamiltonian_initial_guess4()')
        #get hamiltonian_guess from some separate calculation and combine them.
        wfs = self.wfs
        self.gd.comm.broadcast(wfs.S_qMM, 0)
        self.gd.comm.broadcast(wfs.T_qMM, 0)        
        s_pkmm = wfs.S_qMM.copy()
        for S_MM in s_pkmm:
            tri2full(S_MM)
        h_spkmm = self.scale_and_combine_hamiltonian()          
        nb = s_pkmm.shape[-1]
        dtype = s_pkmm.dtype
        for q in range(self.my_npk):
            self.hsd.reset(0, q, s_pkmm[q], 'S', True)                
            for s in range(self.my_nspins):
                self.hsd.reset(s, q, h_spkmm[s, q], 'H', True)
                self.hsd.reset(s, q, np.zeros([nb, nb], dtype), 'D', True)

    def copy_mixer_history(self, calc, guess_type='normal'):
        self.log('copy_mixer_history()')
        if hasattr(self.density.mixer, 'mixers'):  #for Mixer
            mixers = calc.density.mixer.mixers
            mixers0 = self.density.mixer.mixers
        elif hasattr(self.density.mixer, 'mixer'):   #for MixerDif
            mixers = [calc.density.mixer.mixer, calc.density.mixer.mixer_m]
            mixers0 = [self.density.mixer.mixer, self.density.mixer.mixer_m]
        else:  #for BroydenMixer or BroydenMixerSum
            mixers = [calc.density.mixer]
            mixers0 = [self.density.mixer]
        
        if guess_type == 'buffer':
            from gpaw.transport.tools import cut_grids_side, \
                               collect_and_distribute_atomic_matrices
            gd = calc.wfs.gd
            gd0 = self.wfs.gd
            setups = calc.density.setups
            setups0 = self.density.setups
            rank_a = calc.density.rank_a
            keys = self.density.D_asp.keys()
            for mixer, mixer0 in zip(mixers, mixers0):
                for nt_G, R_G, D_ap, dD_ap in zip(mixer.nt_iG,
                                                  mixer.R_iG,
                                              mixer.D_iap,
                                              mixer.dD_iap):
                    nt_G0 = cut_grids_side(nt_G, gd, gd0)
                    R_G0 = cut_grids_side(R_G, gd, gd0)
                    lD_ap = collect_and_distribute_atomic_matrices(D_ap,
                                                        setups, setups0,
                                                       rank_a, gd.comm, keys)
                    ldD_ap = collect_and_distribute_atomic_matrices(dD_ap,
                                                        setups, setups0,
                                                       rank_a, gd.comm, keys)
                    mixer0.nt_iG.append(nt_G0)
                    mixer0.R_iG.append(R_G0)
                    mixer0.D_iap.append(lD_ap)
                    mixer0.dD_iap.append(ldD_ap)
                nt_G0 = cut_grids_side(mixer.nt_iG[-1], gd, gd0)    
                lD_ap = collect_and_distribute_atomic_matrices(mixer.D_iap[-1],
                                                    setups, setups0,
                                                    rank_a, gd.comm, keys)
                mixer0.nt_iG.append(nt_G0)
                mixer0.D_iap.append(lD_ap)
                if hasattr(mixer, 'A_ii'):
                    mixer0.A_ii = mixer.A_ii    
                else:
                    for cG, vG, uG, uD, etaD in zip(mixer.c_G,
                                              mixer.v_G,
                                              mixer.u_G,
                                              mixer.u_D,
                                              mixer.eta_D):
                        vG0 = cut_grids_side(vG, gd, gd0)
                        uG0 = cut_grids_side(uG, gd, gd0)
                        uD0 = collect_and_distribute_atomic_matrices(uD,
                                                            setups, setups0,
                                                       rank_a, gd.comm, keys)
                        etaD0 = collect_and_distribute_atomic_matrices(etaD,
                                                            setups, setups0,
                                                       rank_a, gd.comm, keys)
                        mixer0.c_G.append(cG)
                        mixer0.v_G.append(vG0)
                        mixer0.u_G.append(uG0)
                        mixer0.u_D.append(uD0)
                        mixer0.eta_D.append(etaD0)
                    mixer0.step = mixer.step    
            
        else:
            for mixer, mixer0 in zip(mixers, mixers0):
                mixer0.nt_iG = mixer.nt_iG[:]
                mixer0.R_iG = mixer.R_iG[:]
                mixer0.D_iap = mixer.D_iap[:]
                mixer0.dD_iap = mixer.dD_iap[:]
                if hasattr(mixer, 'A_ii'):
                    mixer0.A_ii = mixer.A_ii
                else:
                    mixer0.c_G = mixer.c_G[:]
                    mixer0.v_G = mixer.v_G[:]
                    mixer0.u_G = mixer.u_G[:]
                    mixer0.u_D = mixer.u_D[:]
                    mixer0.eta_D = mixer.eta_D[:]
                    mixer0.step = mixer.step
       
    def scale_and_combine_hamiltonian(self):
        self.log('scale_and_combine_hamiltonian()')
        #assert self.cell_ham_file is not None
        #fd = file(self.cell_ham_file, 'r')
        #cell_ham_data = cPickle.load(fd)
        #cell_s_pkmm, cell_cs_pkmm, cell_h_spkmm, cell_ch_spkmm, fermi = cell_ham_data[self.wfs.kd.comm.rank]
        #fd.close()

        if self.cell_atoms is not None:
            self.cell_atoms.get_potential_energy()
            cell_h_skmm, cell_s_kmm = self.get_hs(self.cell_atoms.calc)
            cell_d_skmm = get_lcao_density_matrix(self.cell_atoms.calc)

            cell_h_spkmm, cell_s_pkmm, cell_d_spkmm,  \
            cell_ch_spkmm, cell_cs_pkmm, cell_cd_spkmm = get_pk_hsd(self.d, self.ntklead,
                                                    self.cell_atoms.calc.wfs.kd.ibzk_qc,
                                                    cell_h_skmm, cell_s_kmm, cell_d_skmm,
                                                    self.text, self.wfs.dtype,
                                                    direction=0)
            fermi = self.cell_atoms.calc.get_fermi_level()

            efloat = fermi - self.lead_fermi[0]
            cell_h_spkmm -= cell_s_pkmm * efloat
            cell_ch_spkmm -= cell_cs_pkmm * efloat

            nn, nbp = self.sf #nn is the scaling factor, nbp is the fractional part, - means in front, + means at the end
            dtype = cell_h_spkmm.dtype
            nbc = cell_h_spkmm.shape[-1]
            nbp = int(nbp * nbc)
        else:
            nbc = 0
            nbp = 0
            nn = 0
            dtype = self.wfs.dtype

        ns, npk = self.my_nspins, self.my_npk
        nb0, nb1 = self.nblead
        nbp0, nbp1 = self.lead_sf #nbp0, nbp1 are both positive
        nbp0 = int(nbp0 * self.nblead[0])
        nbp1 = int(nbp1 * self.nblead[1])
        assert dtype == self.wfs.dtype
        nb = nbc * nn + abs(nbp) + nb0 + nb1 + abs(nbp0) + abs(nbp1)
        h_spkmm = np.zeros([ns, npk, nb, nb], dtype=dtype)
        for s in range(ns):
            for pk in range(npk):
                h_spkmm[s, pk, :nb0, :nb0] = self.lead_hsd[0].H[s][pk].recover()
                h_spkmm[s, pk, nb0:nb0+nbp0, nb0:nb0+nbp0] = self.lead_hsd[0].H[s][pk].recover()[:nbp0,:nbp0]
                h_spkmm[s, pk, nb0:nb0+nbp0, :nb0] = self.lead_couple_hsd[0].H[s][pk].recover()[:,:nbp0].T.conj()
                h_spkmm[s, pk, :nb0, nb0:nb0+nbp0] = self.lead_couple_hsd[0].H[s][pk].recover()[:,:nbp0]
                if nbc != 0:
                    if nbp < 0:
                        a, b = nb0+nbp0, nb0+nbp0-nbp 
                        h_spkmm[s, pk, a:b, a:b] = cell_h_spkmm[s, pk, nbp:,nbp:]
                        h_spkmm[s, pk, a:b, b:b+nbc] = cell_ch_spkmm[s, pk, nbp:,:]            
                        h_spkmm[s, pk, b:b+nbc, a:b] = cell_ch_spkmm[s, pk, nbp:,:].T.conj()            
                        a, b = nb0+nbp0-nbp, nb0+nbp0-nbp+nbc
                        for i in range(nn):
                            h_spkmm[s, pk, a:b, a:b] = cell_h_spkmm[s, pk]
                            a += nbc
                            b += nbc
                        a, b = nb0+nbp0-nbp, nb0+nbp0-nbp+nbc
                        for i in range(nn-1):
                            h_spkmm[s, pk, a:b, a+nbc:b+nbc] = cell_ch_spkmm[s, pk]
                            h_spkmm[s, pk, a+nbc:b+nbc, a:b] = cell_ch_spkmm[s, pk].T.conj()
                            a += nbc
                            b += nbc
                    else:
                        a, b = nb0+nbp0, nb0+nbp0+nbc
                        for i in range(nn):
                            h_spkmm[s, pk, a:b, a:b] = cell_h_spkmm[s, pk]
                            a += nbc
                            b += nbc
                        a, b = nb0+nbp0, nb0+nbp0+nbc
                        for i in range(nn-1):
                            h_spkmm[s, pk, a:b, a+nbc:b+nbc] = cell_ch_spkmm[s, pk]
                            h_spkmm[s, pk, a+nbc:b+nbc, a:b] = cell_ch_spkmm[s, pk].T.conj()
                            a += nbc
                            b += nbc
                        a, b = nb0 + nbp0 + nn *nbc, nb0 + nbp0 + nn *nbc +nbp
                        h_spkmm[s, pk, a:b, a:b] = cell_h_spkmm[s, pk, :nbp, :nbp]
                        h_spkmm[s, pk, a:b, a-nbc:a] = cell_ch_spkmm[s, pk, :, :nbp].T.conj()
                        h_spkmm[s, pk, a-nbc:a, a:b] = cell_ch_spkmm[s, pk, :, :nbp]
                else:
                    b = nb0 + nbp0
                h_spkmm[s, pk, b+nbp1:b+nbp1+nb1, b+nbp1:b+nbp1+nb1] = self.lead_hsd[1].H[s][pk].recover()
                if nbp1 != 0:    
                    h_spkmm[s, pk, b:b+nbp1, b:b+nbp1] = self.lead_hsd[1].H[s][pk].recover()[-nbp1:, -nbp1:]
                    h_spkmm[s, pk, b:b+nbp1, b+nbp1:b+nbp1+nb1] = self.lead_couple_hsd[1].H[s][pk].recover()[:, -nbp1:].T.conj()
                    h_spkmm[s, pk, b+nbp1:b+nbp1+nb1, b:b+nbp1] = self.lead_couple_hsd[1].H[s][pk].recover()[:, -nbp1:]
        return h_spkmm 

    def initialize_hamiltonian_matrix(self, calc):
        self.log('initialize_hamiltonian_matrix()')
        h_skmm, s_kmm =  self.get_hs(calc)
        d_skmm = get_lcao_density_matrix(calc)
        ntk = self.scat_ntk
        kpts = calc.wfs.kd.ibzk_qc
        h_spkmm = substract_pk(self.d, self.my_npk, ntk, kpts, h_skmm, 'h')
        s_pkmm = substract_pk(self.d, self.my_npk, ntk, kpts, s_kmm)
        d_spkmm = substract_pk(self.d, self.my_npk, ntk, kpts, d_skmm, 'h')
        h00 = self.lead_hsd[0].H[0][0].recover()[0,0]
        if not self.buffer_guess:
            h01 = h_spkmm[0,0,0,0]
        else:
            #nb = self.nblead[0]
            nb =0
            h01 = h_spkmm[0,0,nb,nb]            
        s00 = self.lead_hsd[0].S[0].recover()[0,0]
        e_shift = (h00 - h01) / s00
        h_spkmm += e_shift * s_pkmm
   
        if self.wfs.dtype == float:
            h_spkmm = np.real(h_spkmm).copy()
            s_pkmm = np.real(s_pkmm).copy()
            d_spkmm = np.real(d_spkmm).copy()
       
        for q in range(self.my_npk):
            self.hsd.reset(0, q, s_pkmm[q], 'S', True)
            for s in range(self.my_nspins):
                self.hsd.reset(s, q, h_spkmm[s, q], 'H', True)            
                self.hsd.reset(s, q, d_spkmm[s, q] * ntk, 'D', True)
            
    def fill_guess_with_leads(self, flag=None):
        self.log('fill_guess_with_leads()')
        if self.hsd.S[0].extended:
            n = -2
        else:
            n = -1
        for s in range(self.my_nspins):
            for pk in range(self.my_npk):
                for l in range(self.lead_num):
                    self.hsd.H[s][pk].diag_h[l][n].reset(
                                          self.lead_hsd[l].H[s][pk].recover()
                            + self.bias[l] * self.lead_hsd[l].S[pk].recover())
                    if flag != 'H':
                        self.hsd.D[s][pk].diag_h[l][n].reset(
                                          self.lead_hsd[l].D[s][pk].recover())

    def append_buffer_hsd(self):
        self.log('append_buffer_hsd()')
        tp_mat = self.hsd.S[0]
        if tp_mat.extended:
            ex_index = []
            sum = 0
            for i in range(self.lead_num):
                self.log('  append lead: {0}'.format(i))
                begin = np.min(self.lead_index[i])
                newb = tp_mat.nb + sum
                #ex_index.append(self.lead_index[i] - begin + newb)
                ex_index.append(np.argsort(self.lead_index[i]) + newb)
                sum += self.nblead[i]
                #ex_index = [self.lead_index[0] + tp_mat.nb]
                #ex_index.append(self.lead_index[1] +
                #                       self.nblead[0] + self.nblead[1] +
                #                       tp_mat.nb - self.lead_index[1][-1] - 1)
            self.log('  call hsd.append_lead_as_buffer()')
            self.hsd.append_lead_as_buffer(self.lead_hsd,
                                           self.lead_couple_hsd,
                                           ex_index,
                                           tp=self)
            self.log('  finished hsd.append_lead_as_buffer()')
            
                 
    def get_basis_indices(self):
        self.log('get_basis_indices()')
        setups = self.wfs.setups
        
        edge_index = []
        lead_edge_index = []
        
        for i in range(self.lead_num):
            lead_setups = self.setups_lead[i]
            self.lead_index[i] = get_atom_indices(self.pl_atoms[i], setups)
            edge_index.append(get_atom_indices([self.edge_atoms[1][i]], setups))
            lead_edge_index.append(get_atom_indices([self.edge_atoms[0][i]],
                                                                 lead_setups))            
            self.edge_index[0][i] = lead_edge_index[i][0]
            self.edge_index[1][i] = edge_index[i][0]
        
        edge_index = []
        self.orbital_indices = []
        for n, setup in enumerate(setups):
            for phit in setup.phit_j:
                l = phit.get_angular_momentum_number()
                for i in range(2 * l  + 1):
                    self.orbital_indices.append([n, l])
        self.orbital_indices = np.array(self.orbital_indices)
        self.lead_orbital_indices = []
        for i in range(self.lead_num):
            self.lead_orbital_indices.append([])
            for n, setup in enumerate(self.setups_lead[i]):
                for phit in setup.phit_j:
                    l = phit.get_angular_momentum_number()
                    for j in range(2 * l  + 1):
                        self.lead_orbital_indices[i].append([n, l])                    
            self.lead_orbital_indices[i] = np.array(self.lead_orbital_indices[i])
        
        for i in range(self.lead_num):
            if self.la_index is None:
                n_layer_atoms = len(self.lead_atoms[i]) // self.nleadlayers[i]
                self.lead_layer_index[i][0] = get_atom_indices(self.mol_atoms, setups)
                begin = 0
                for j in range(1, self.nleadlayers[i] + 1):
                    atoms_index = self.lead_atoms[i][begin: begin + n_layer_atoms]
                    self.lead_layer_index[i][j] = get_atom_indices(atoms_index, setups)
                    begin += n_layer_atoms
            else:
                self.lead_layer_index[i][0] = get_atom_indices(self.mol_atoms, setups)                
                for j in range(1, self.nleadlayers[i] + 1):
                    self.lead_layer_index[i][j] = get_atom_indices(self.la_index[i][j - 1], setups)
      
    def initialize_lead_matrix(self):
        self.log('initialize_lead_matrix()')
        if self.use_lead:
            self.lead_hsd = []
            self.lead_couple_hsd = []
            self.lead_index = []
            self.inner_lead_index = []
            self.buffer_index = []
            self.lead_layer_index = []
            self.lead_fermi = np.empty([self.lead_num])

        npk = self.my_npk
        ns = self.my_nspins
        dtype = self.wfs.dtype
       
        for i in range(self.lead_num):
            self.lead_hsd.append(Banded_Sparse_HSD(dtype, ns, npk))
            self.lead_couple_hsd.append(CP_Sparse_HSD(dtype, ns, npk))
            self.lead_index.append([])
            self.inner_lead_index.append([])
            self.buffer_index.append([])
            self.lead_layer_index.append([])
            for j in range(self.nleadlayers[i] + 1):
                self.lead_layer_index[i].append([])
      
        if self.use_lead:
            self.ec = np.zeros([self.lead_num, ns])

    def initialize_matrix(self):
        self.log('initialize_matrix()')
        dtype = self.wfs.dtype
        extended = True
        self.hsd = Tp_Sparse_HSD(dtype, self.my_nspins, self.my_npk,
                                              self.lead_layer_index, extended)              

    def distribute_energy_points(self):
        self.log('distribute_energy_points()')
        self.energy_comm = self.gd.comm
        size, rank = self.energy_comm.size, self.energy_comm.rank
        ns, npk = self.my_nspins, self.my_npk
        self.eq_par_energy_index = []
        self.ne_par_energy_index = []  
 
        for s in range(ns):
            self.eq_par_energy_index.append([])
            self.ne_par_energy_index.append([])
            for k in range(npk):
                self.eq_par_energy_index[s].append([])
                self.ne_par_energy_index[s].append([])
                     
        if not self.ground:
            self.loc_par_energy_index = []
            for s in range(ns):
                self.loc_par_energy_index.append([])
                for k in range(npk):
                    self.loc_par_energy_index[s].append([])            
        
        for s in range(ns):
            for k in range(npk):
                neeq = self.eqpathinfo[s][k].num
                eq_ind = np.array_split(np.arange(neeq), size)
                self.eq_par_energy_index[s][k] = eq_ind[rank]

                if not self.ground:
                    neloc = self.locpathinfo[s][k].num
                    loc_ind = np.array_split(np.arange(neloc), size)
                    self.loc_par_energy_index[s][k] = loc_ind[rank]

                nene = self.nepathinfo[s][k].num
                ne_ind = np.array_split(np.arange(nene), size)
                self.ne_par_energy_index[s][k] = ne_ind[rank]

    def collect_leads_matrices(self, calc, l):
        self.log('collect_leads_matrices()')
        hl_skmm, sl_kmm = self.get_hs(calc)
        dl_skmm = get_lcao_density_matrix(calc)
        lead_direction = l # character l
        hl_spkmm, sl_pkmm, dl_spkmm,  \
        hl_spkcmm, sl_pkcmm, dl_spkcmm = get_pk_hsd(self.d, self.ntklead,
                                                calc.wfs.kd.ibzk_qc,
                                                hl_skmm, sl_kmm, dl_skmm,
                                                self.text, self.wfs.dtype,
                                                direction=lead_direction)
        for pk in range(self.my_npk):
            self.lead_hsd[l].reset(0, pk, sl_pkmm[pk], 'S', init=True)
            self.lead_couple_hsd[l].reset(0, pk, sl_pkcmm[pk], 'S',
                                                              init=True)
            for s in range(self.my_nspins):
                self.lead_hsd[l].reset(s, pk, hl_spkmm[s, pk], 'H', init=True)     
                self.lead_hsd[l].reset(s, pk, dl_spkmm[s, pk], 'D', init=True)
                
                self.lead_couple_hsd[l].reset(s, pk, hl_spkcmm[s, pk],
                                                               'H', init=True)     
                self.lead_couple_hsd[l].reset(s, pk, dl_spkcmm[s, pk],
                                                               'D', init=True)                    

    def recover_lead_info(self, s00, s01, h00, h01, fermi):
        self.log('recover_lead_info()')
        for pk in range(self.npk):
            for l in range(self.lead_num):
                self.lead_hsd[l].reset(0, pk, s00[pk, l], 'S', init=True)
                self.lead_couple_hsd[l].reset(0, pk, s01[pk, l], 'S',
                                                              init=True)
                for s in range(self.nspins):
                    self.lead_hsd[l].reset(s, pk, h00[s, pk, l], 'H', init=True)     
                    self.lead_hsd[l].reset(s, pk, np.zeros_like(h00[s, pk, l]), 'D', init=True)
                
                    self.lead_couple_hsd[l].reset(s, pk, h01[s, pk, l],
                                                               'H', init=True)     
                    self.lead_couple_hsd[l].reset(s, pk, np.zeros_like(h01[s, pk, l]),
                                                               'D', init=True)                    
        for l in range(self.lead_num):
            self.lead_fermi[l] = fermi[l]
             
    def update_scat_hamiltonian(self, atoms):
        self.log('update_scat_hamiltonian()')
        if atoms is None:
            atoms = self.atoms
        else:
            self.atoms = atoms.copy()
            #self.initialize()
            self.get_extended_atoms()
            self.density.reset()
            self.set_extended_positions()
            #del self.wfs
            #self.wfs = self.extended_calc.wfs
                
        #if self.scat_restart:
        #    self.recover_kpts(self)
        if not self.optimize:
            self.append_buffer_hsd()
        if self.lead_guess:
            self.fill_guess_with_leads('H')           
        self.scat_restart = False

    def get_hs(self, calc):
        self.log('get_hs()')
        wfs = calc.wfs
        eigensolver = wfs.eigensolver
        ham = calc.hamiltonian
        self.gd.comm.broadcast(wfs.S_qMM, 0)
        self.gd.comm.broadcast(wfs.T_qMM, 0)        
        S_qMM = wfs.S_qMM.copy()
        for S_MM in S_qMM:
            tri2full(S_MM)
        H_sqMM = np.empty((self.my_nspins,) + S_qMM.shape, wfs.dtype)
        for kpt in wfs.kpt_u:
            H_MM = eigensolver.calculate_hamiltonian_matrix(ham, wfs, kpt)
            tri2full(H_MM)
            H_MM *= Hartree
            if self.my_nspins == 2:
                H_sqMM[kpt.s, kpt.q] = H_MM
            else:
                H_sqMM[0, kpt.q] = H_MM
        return H_sqMM, S_qMM
       
    def get_lead_atoms(self, l, init_calc=True):
        self.log('get_lead_atoms()')
        """Here is a multi-terminal version """
        if self.leads is not None:
            atomsl = self.leads[l]
        else:
            atoms = self.atoms.copy()
            atomsl = atoms[self.pl_atoms[l]]
            atomsl.cell = self.pl_cells[l]
            atomsl.center(axis=2)
            atomsl._pbc[self.d] = True
        if init_calc:
            atomsl.set_calculator(self.get_lead_calc(l))
        return atomsl
    
    def get_lead_calc(self, l):
        self.log('get_lead_calc()')
        p = self.gpw_kwargs.copy()
        if type(p['basis']) is dict and len(p['basis'])==len(self.atoms):
            basis = {}
            for i, a in enumerate(self.pl_atoms[l]):
                basis[i] = p['basis'][a]
            p['basis'] = basis
        if 'setups' in p:
            if type(p['setups']) is dict and len(p['setups'])==len(self.atoms):
                setups = {}
                for i, a in enumerate(self.pl_atoms[l]):
                    setups[i] = p['setups'][a]
                p['setups'] = setups
        p['nbands'] = None
        p['kpts'] = self.pl_kpts
        p['parallel'] = self.input_parameters['parallel'].copy()
        p['parallel'].update(band=self.wfs.bd.comm.size,
                             # kpt=self.wfs.kd.comm.size
                             domain=self.wfs.gd.parsize_c)
        if 'mixer' in p:
            if not self.spinpol:
                p['mixer'] = Mixer(0.1, 5, weight=100.0)
            else:
                p['mixer'] = MixerDif(0.1, 5, weight=100.0)
        p['poissonsolver'] = PoissonSolver(nn=2)
        if 'txt' in p and p['txt'] != '-':
            p['txt'] = 'lead%i_' % (l + 1) + p['txt']
        return Lead_Calc(**p)

    def negf_prepare(self, atoms=None):
        self.log('negf_prepare()')
        if not self.initialized_transport:
            self.initialize_transport()
        if not self.analysis_mode:    
            self.update_scat_hamiltonian(atoms)
        #if self.ground and not self.analysis_mode:
        #    self.boundary_align_up()

    def boundary_align_up(self):
        self.log('boundary_align_up()')
        tol = 0.1
        ind = self.edge_index[0][0]
        level_in_lead = self.lead_hsd[0].H[0][0].recover()[ind, ind]
        ind = self.edge_index[1][0]
        level_in_scat = self.hsd.H[0][0].recover()[ind, ind]
        overlap_on_site = self.hsd.S[0].recover()[ind, ind]
        shift = (level_in_scat - level_in_lead) / overlap_on_site
        if abs(shift) > tol:
            for s in range(self.my_nspins):
                for pk in range(self.my_npk):
                    self.hsd.H[s][pk].reset_from_others(self.hsd.H[s][pk],
                                                    self.hsd.S[pk], 1, -shift)
        self.align_shift = shift
       
    def get_selfconsistent_hamiltonian(self):
        self.log('get_selfconsistent_hamiltonian()')
        self.timer.start('init scf')
        if not self.fix_contour or not self.optimize:
            self.initialize_scf()
        else:
            self.initialize_scf_flags()
        self.timer.stop('init scf')
        
        ##temperary lines
        #self.hamiltonian.S = 0
        #self.hamiltonian.Etot = 0
        ##temp
 
        if self.analysor.n_bias_step in self.perturbation_steps:
            #self.density.mixer.reset()
            self.induce_density_perturbation()
     
        while not self.cvgflag and self.step < self.max_steps:
            self.iterate()
            self.cvgflag = self.d_cvg and self.h_cvg
            self.step +=  1
        
        if self.foot_print:
            self.analysor.save_bias_step(self)
            magmom = self.occupations.magmom
            self.text('Total Magnetic Moment: %f' % magmom)
            self.text('Local Magnetic Moments:')
            try:
                for a, mom in enumerate(self.get_magnetic_moments()):
                    self.text(a, mom)
            except:
                pass
        self.scf.converged = self.cvgflag
   
        if self.save_bias_data:
            vt_sG = self.gd1.collect(self.extended_calc.hamiltonian.vt_sG)
            if not self.use_qzk_boundary:
                density = self.density
                ham = self.hamiltonian
            else:
                density = self.extended_calc.density
                ham = self.extended_calc.hamiltonian                
            dH_asp = collect_atomic_matrices(ham.dH_asp, ham.setups,
                                             ham.nspins, ham.gd.comm,
                                             density.rank_a)
            if self.master:
                fd = file('bias_data' + str(self.analysor.n_bias_step), 'wb')
                cPickle.dump((self.bias, vt_sG, dH_asp), fd, 2)
                fd.close()
                
        self.ground = False
        self.alpha = 0
        self.total_charge = 0
        self.linear_mm = None
        if not self.scf.converged:
            raise RuntimeError('Transport do not converge in %d steps' %
                                                              self.max_steps)
        
    def non_sc_analysis(self):
        self.log('non_sc_analysis()')
        if not hasattr(self, 'contour'):
            self.contour = Contour(0.1,
                               self.lead_fermi, self.bias, comm=self.gd.comm,
                                plot_eta=self.plot_eta, eta=self.eta,
                                plot_energy_range=self.plot_energy_range,
                             plot_energy_point_num=self.plot_energy_point_num)
            
        if not hasattr(self, 'analysor'):
            self.analysor = Transport_Analysor(self, True)
            
        #self.analysor.save_ele_step()            
        self.analysor.save_bias_step(self)
        fd = file('eq_hsd', 'w')
        cPickle.dump(self.hsd, fd, 2)
        fd.close()

    def get_density_matrix(self):
        self.log('get_density_matrix()')
        self.timer.start('DenMM')
        if self.use_qzk_boundary:
            self.fill_lead_with_scat()
            #for i in range(self.lead_num):
                #self.selfenergies[i].set_bias(0)
        
        if self.recal_path:
            #self.initialize_path()
            self.calculate_integral_path()
            
        for s in range(self.my_nspins):
            for k in range(self.my_npk):
                #if self.recal_path:
                    #d_mm = self.get_eqintegral_points(s, k) + \
                    #                          self.get_neintegral_points(s, k)
                #else:
                    #d_mm = self.fock2den(s, k)
                d_mm = self.fock2den(s, k)                    
                d_mm = self.spin_coff * (d_mm + d_mm.T.conj()) / (2 * self.npk)
                if self.gate_mode == 'AN':
                    d_mm_gate_plus = self.gate_filling()
                    ind = get_matrix_index(self.gate_basis_index)
                    d_mm[ind.T, ind] += d_mm_gate_plus
                self.hsd.reset(s, k, d_mm, 'D') 
        self.timer.stop('DenMM')
        
        self.print_boundary_charge()
        if self.master:
            self.text('DenMM', self.timer.timers['DenMM', ], 'second')

    def gate_filling(self):
        self.log('gate_filling()')
        ind = get_matrix_index(self.gate_basis_index)
        sub_overlap = self.hsd.S[0].recover()[ind.T, ind]
        unit_charge = np.trace(sub_overlap)        
        dmm_plus = self.gate / unit_charge * np.eye(len(ind))
        return dmm_plus

    def iterate(self):
        self.log('iterate()')
        if self.master:
            self.text('----------------step %d -------------------'
                                                                % self.step)
        self.h_cvg = self.check_convergence('h')
        self.get_density_matrix()
       
        self.timer.start('HamMM')            
        self.get_hamiltonian_matrix()
        self.timer.stop('HamMM')
        if self.master:
            self.text('HamMM', self.timer.timers['HamMM',], 'seconds')         
       
        self.d_cvg = self.check_convergence('d')
        
        # adaptive mixing ?
        if self.h_cvg:
            self.cvg_ham_steps += 1
        else:
            self.cvg_ham_steps = 0
        
        if self.cvg_ham_steps > 4:
            self.text('Ham cvg since {0} iterations -> increase density mixing'.format(self.cvg_ham_steps))
            b = self.density.mixer.beta
            self.density.mixer.beta = b + (b*0.25)
            self.text('New beta: {0}'.format(self.density.mixer.beta))
        
        self.txt.flush()
        
    def check_convergence(self, var):
        self.log('check_convergence()')
        cvg = False
        if var == 'h':
            diag_ham = np.zeros([self.nbmol], self.wfs.dtype)
            for s in range(self.my_nspins):
                for q in range(self.my_npk):
                    diag_ham += np.diag(self.hsd.H[s][q].recover())
            self.wfs.kd.comm.sum(diag_ham)
            diag_ham /= self.npk
            
            self.diff_h = 1.         
            if self.step > 0:
                self.diff_h = np.max(abs(diag_ham - self.diag_ham_old))
                if self.master:
                    self.text('hamiltonian: diff = %f  tol=%f' % (self.diff_h,
                                                  self.diag_ham_tol))
                if self.diff_h < self.diag_ham_tol:
                    cvg = True
            self.diag_ham_old = np.copy(diag_ham)
        if var == 'd':
            if self.step > 0:
                if not self.use_qzk_boundary:
                    density = self.density
                else:
                    density = self.extended_calc.density
                self.diff_d = density.mixer.get_charge_sloshing()
                tol =  self.scf.max_density_error
 
                if self.master:
                    self.text('density: diff = %f  tol=%f' % (self.diff_d,
                                            tol))
                if self.diff_d < tol * self.theta or (self.neutral_steps is
                                not None and self.step > self.neutral_steps):
                    if (self.use_qzk_boundary or self.fixed) and \
                                  not self.normalize_density and self.neutral:
                        self.neutral = False
                    elif self.diff_d < tol:
                        cvg = True
        return cvg
 
    def initialize_scf(self):
        self.log('initialize_scf()')
        self.contour = Contour(self.occupations.width * Hartree,
                               self.lead_fermi, self.bias, comm=self.gd.comm,
                               plot_eta=self.plot_eta, eta=self.eta,
                               neintstep=self.neintstep,
                               eqinttol=self.eqinttol,
                               min_energy=self.min_energy,
                               plot_energy_range=self.plot_energy_range,
                             plot_energy_point_num=self.plot_energy_point_num)

        self.surround.reset_bias(self)
        #if not self.use_qzk_boundary:
        #    self.surround.reset_bias(self)
        #else:
        #    self.surround.reset_bias(self)
        self.initialize_green_function()
        self.calculate_integral_path()
        self.distribute_energy_points()
    
    
        if self.master:
            self.text('------------------Transport SCF-----------------------') 
            bias_info = 'Bias:'
            for i in range(self.lead_num):
                bias_info += 'lead' + str(i) + ': ' + str(self.bias[i]) + 'V'
            self.text(bias_info)
            self.text('Gate (Gate_Mode): %f V %s' % (self.gate, self.gate_mode))

        if not hasattr(self, 'analysor'):
            self.analysor = Transport_Analysor(self)
        
        if not self.fixed and not self.use_buffer:
            self.get_linear_hartree_potential()
        #------for check convergence------
        #self.ham_vt_old = np.empty(self.hamiltonian.vt_sG.shape)
        self.initialize_scf_flags()

    def initialize_scf_flags(self):
        self.log('initialize_scf_flags()')
        self.ham_vt_diff = None
        self.ham_vt_tol = 1e-2
        #self.diag_ham_tol = 5e-3
        self.diag_ham_tol = self.scf.max_energy_error * Hartree / len(self.atoms) * 5
        self.step = 0
        self.cvgflag = False
        self.spin_coff = 3. - self.nspins
        self.max_steps = 400
        self.h_cvg = False
        self.d_cvg = False
        self.ele_data = {}
        
    def initialize_path(self):
        self.log('initialize_path()')
        self.eqpathinfo = []
        self.nepathinfo = []
        self.locpathinfo = []
        for s in range(self.my_nspins):
            self.eqpathinfo.append([])
            self.nepathinfo.append([])
            if not self.ground:
                self.locpathinfo.append([])                
            if not self.ground:
                self.locpathinfo.append([])
            for k in range(self.my_npk):
                self.eqpathinfo[s].append(EnergyNode('eq', self.lead_num))
                self.nepathinfo[s].append(EnergyNode('ne', self.lead_num))    
                if not self.ground:
                    self.locpathinfo[s].append(EnergyNode('eq',
                                                         self.lead_num))
                    
    def calculate_integral_path(self):
        self.log('calculate_integral_path()')
        self.initialize_path()
        for s in range(self.my_nspins):
            for k in range(self.my_npk):      
                self.find_contour(s, k)  
        ne = self.eqpathinfo[0][0].num + self.nepathinfo[0][0].num
        if not self.ground:
            ne += self.locpathinfo[0][0].num
        self.text('energy point' + str(ne))

    def find_contour(self, s, k):
        self.log('find_contour()')
        self.cntint = -1
        self.fint = []
        self.tgtint = []
        for i in range(self.lead_num):
            self.tgtint.append([])
        self.zint = [0] * 500
                
        self.reset_lead_hs(s, k)        
        self.hsd.s = s
        self.hsd.pk = k        
        self.contour.get_optimized_contour(self)
        nids, energies, weights, ses = self.contour.sort_contour(self.nbmol)
        elist = []
        wlist = []
        flist = []
        
        siglist = []
        for i in range(self.lead_num):
            siglist.append([])        

        elist1 = []
        wlist1 = []
        flist1 = []
        
        siglist1 = []
        for i in range(self.lead_num):
            siglist1.append([])      
        
        kt = self.contour.kt
        max_ef = np.max(self.contour.leadfermi)
        min_ef = np.min(self.contour.leadfermi)
        for nid, energy, weight, se in zip(nids, energies, weights, ses):
            if str(nid)[0] == '6':
                elist1.append(energy)
                wlist1.append(weight)
                for i in range(self.lead_num):
                    siglist1[i].append(se[i])
                flist1.append(fermidistribution(energy - max_ef, kt) - 
                                      fermidistribution(energy - min_ef, kt) )            
            else:
                elist.append(energy)
                wlist.append(weight)
                for i in range(self.lead_num):
                    siglist[i].append(se[i])  
                flist.append(fermidistribution(energy - min_ef, kt))            
                

        comm = self.gd.comm
        ne_poles = np.array_split(np.arange(4), comm.size)
        myne = ne_poles[comm.rank]
        
        eq_poles = np.arange(1, 8, 2) * np.pi * 1.j * kt + min_ef
        if not myne.tolist() == []:
            my_eq_poles = eq_poles[myne]
        else:
            my_eq_poles = np.array([], int)
        
        my_eq_ffp = [-2.j * np.pi * kt] * len(myne)
        my_eq_wp = [1] * len(myne)
        #my_eq_ses = [[], []]
        for i in range(self.lead_num):
            for e in my_eq_poles:
                siglist[i].append(self.selfenergies[i](e))

        
        elist += my_eq_poles.tolist()
        wlist += my_eq_wp
        flist += my_eq_ffp
        #siglist += my_eq_ses

        self.eqpathinfo[s][k].add(elist, wlist, flist, siglist)
        
        if not self.ground:
            comm = self.gd.comm
            ne_poles = np.array_split(np.arange(8), comm.size)
            myne = ne_poles[comm.rank]            
            
            loc_poles1 = np.arange(1, 8, 2) * np.pi * 1.j * kt + min_ef
            loc_poles2 = np.arange(1, 8, 2) * np.pi * 1.j * kt + max_ef
            loc_poles = np.append(loc_poles1, loc_poles2)
            
            if not myne.tolist() == []:
                my_loc_poles = loc_poles[myne]
            else:
                my_loc_poles = np.array([], int)
            my_loc_ffp = [-2.j * np.pi * kt] * len(myne)
            loc_wp = [-1.] * 4 + [1.] * 4
            my_loc_wp = np.array_split(loc_wp, comm.size)[comm.rank]
            
            for i in range(self.lead_num):
                for e in my_loc_poles:
                    siglist1[i].append(self.selfenergies[i](e))
          
            elist1 += my_loc_poles.tolist()
            wlist1 += my_loc_wp.tolist()
            flist1 += my_loc_ffp
            #siglist1 += my_loc_ses
            self.locpathinfo[s][k].add(elist1, wlist1, flist1, siglist1)              


            nids, energies, weights = self.contour.distribute_nodes(4)
            elist2 = energies.tolist()
            wlist2 = weights.tolist()
            flist2 = []
            siglist2 = [[], []]
            for i in range(self.lead_num):
                flist2.append([[], []])
                for e in elist2:
                    flist2[i][0].append(fermidistribution(e - self.contour.leadfermi[i],
                                           kt) - fermidistribution(e -
                                          min_ef, kt)) 
            
                    flist2[i][1].append(fermidistribution(e - max_ef,
                                           kt) - fermidistribution(e -
                                            self.contour.leadfermi[i], kt))  
                    siglist2[i].append(self.selfenergies[i](e))
            self.nepathinfo[s][k].add(elist2, wlist2, flist2, siglist2)        
        self.contour.release()
       
    def calgfunc(self, zp, calcutype, flag='old'):
        #self.log('-> calgfunc()')
        #self.log('   zp: {0}'.format(zp))
        #self.log('   calcutype: {0}'.format(calcutype))
        #self.log('   flag: {0}'.format(flag))
        #calcutype = 
        #  - 'eqInt':  gfunc[Mx*Mx,nE] (default)
        #  - 'neInt':  gfunc[Mx*Mx,nE]
        #  - 'resInt': gfunc[Mx,Mx] = gr * fint
        #              fint = -2i*pi*kt
      
        contour = self.contour
        sgftol = 1e-10
        stepintcnt = 50
        nbmol = self.nbmol_inner
                
        if type(zp) == list:
            pass
        elif type(zp) == np.ndarray:
            pass
        else:
            zp = [zp]
        nume = len(zp)
        
        if calcutype == 'resInt':
            gfunc = np.zeros([nbmol, nbmol], complex)
        elif calcutype == 'neInt':
            gfunc = np.zeros([nume, 2, nbmol, nbmol], complex)            
        else:
            gfunc = np.zeros([nume, nbmol, nbmol], complex)
        for i in range(nume):
            sigma = []

            if self.cntint + 1 >= len(self.zint):
                self.zint += [0] * stepintcnt

            self.cntint += 1
            self.zint[self.cntint] = zp[i]

            for j in range(self.lead_num):
                tgt = self.selfenergies[j](zp[i])
                self.tgtint[j].append(tgt)
            
            for j in range(self.lead_num):
                ind = self.inner_lead_index[j]
                dim = len(ind)
                ind = np.resize(ind, [dim, dim])
                tgt = self.tgtint[j][self.cntint]
                sigma.append(tgt)
            #self.log('   zp[i]: {0}'.format(zp[i]))
            #self.log('   sigma: {0}'.format(sigma))
            gr = self.hsd.calculate_eq_green_function(zp[i], sigma, False)
            # --ne-Integral---
            kt = contour.kt
            ff = []
            if calcutype == 'neInt':
                ffocc = []
                ffvir = []
                for n in range(self.lead_num):
                    lead_ef = contour.leadfermi[n]
                    min_ef = contour.minfermi
                    max_ef = contour.maxfermi
                    self.fint[n][0].append(fermidistribution(zp[i] - lead_ef,
                                           kt) - fermidistribution(zp[i] -
                                          min_ef, kt))
                    self.fint[n][1].append(fermidistribution(zp[i] - max_ef,
                                           kt) - fermidistribution(zp[i] -
                                            lead_ef, kt))                    
                    ffocc.append(self.fint[n][0][self.cntint])
                    ffvir.append(self.fint[n][1][self.cntint])
                    
                gfunc[i, 0], gfunc[i, 1] = \
                                        self.hsd.calculate_ne_green_function(
                                                                 zp[i], sigma,
                                                          ffocc, ffvir, False)
                
            # --local-Integral--
            elif calcutype == 'locInt':
                # fmax-fmin
                max_ef = contour.maxfermi
                min_ef = contour.minfermi
                self.fint.append(fermidistribution(zp[i] - max_ef, kt) - 
                                 fermidistribution(zp[i] - min_ef, kt) )
                gfunc[i] = gr * self.fint[self.cntint]
 
            # --res-Integral --
            elif calcutype == 'resInt':
                self.fint.append(-2.j * np.pi * kt)
                gfunc += gr * self.fint[self.cntint]
            #--eq-Integral--
            else:
                if kt <= 0:
                    self.fint.append(1.0)
                else:
                    min_ef = contour.minfermi
                    self.fint.append(fermidistribution(zp[i] - min_ef, kt))
                gfunc[i] = gr * self.fint[self.cntint]    
        self.log('<- calgfunc()')
        if flag == 'old':
            return gfunc
        else:
            return gfunc, sigma

    def fock2den(self, s, k):
        self.log('fock2den()')
        self.hsd.s = s
        self.hsd.pk = k

        den = self.eq_fock2den(s, k)
        denocc, denvir = self.ne_fock2den(s, k)    
        den += denocc

        if not self.ground:
            denloc = self.eq_fock2den(s, k, el='loc')
            weight_mm = self.integral_diff_weight(denocc, denvir,
                                                                 'transiesta')
            diff = (denloc - (denocc + denvir)) * weight_mm
            den += diff
            percents = np.sum( diff * diff ) / np.sum( denocc * denocc )
            self.text('local percents %f' % np.real(percents))
        
        den = (den + den.T.conj()) / 2
        if self.wfs.dtype == float:
            den = np.real(den).copy()
        return den
    
    def ne_fock2den(self, s, k):
        self.log('ne_fock2den()')
        pathinfo = self.nepathinfo[s][k]
        nbmol = self.nbmol_inner
        denocc = np.zeros([nbmol, nbmol], complex)
        denvir = np.zeros([nbmol, nbmol], complex)
        ind = self.ne_par_energy_index[s][k]
        zp = pathinfo.energy

        self.timer.start('ne fock2den')
        for i in range(len(zp)):
            sigma = []
            for n in range(self.lead_num):
                sigma.append(pathinfo.sigma[n][i])
            ffocc = []
            ffvir = []
            for n in range(self.lead_num):
                ffocc.append(pathinfo.fermi_factor[n][0][i])
                ffvir.append(pathinfo.fermi_factor[n][1][i])
            glesser, ggreater = self.hsd.calculate_ne_green_function(zp[i],
                                                 sigma, ffocc, ffvir, False)
            weight = pathinfo.weight[i]            
            denocc += glesser * weight / np.pi / 2
            denvir += ggreater * weight / np.pi / 2
        self.energy_comm.sum(denocc)
        self.energy_comm.sum(denvir)
        self.timer.stop('ne fock2den')
        return denocc, denvir
    
    def eq_fock2den(self, s, k, el='eq'):
        self.log('eq_fock2den()')
        if el =='loc':
            pathinfo = self.locpathinfo[s][k]
        else:
            pathinfo = self.eqpathinfo[s][k]

        nbmol = self.nbmol_inner
        den = np.zeros([nbmol, nbmol], complex)
        zp = pathinfo.energy
        self.timer.start('eq fock2den')
        for i in range(len(pathinfo.energy)):
            sigma = []
            for n in range(self.lead_num):
                sigma.append(pathinfo.sigma[n][i])
            gr = self.hsd.calculate_eq_green_function(zp[i], sigma, False)
            fermifactor = pathinfo.fermi_factor[i]
            weight = pathinfo.weight[i]
            den += gr * fermifactor * weight
        self.energy_comm.sum(den)
        den = 1.j * (den - den.T.conj()) / np.pi / 2
        self.timer.stop('eq fock2den')
        return den
    
    def get_hamiltonian_matrix(self):
        self.log('get_hamiltonian_matrix()')
        self.update_density()
        if self.use_qzk_boundary:
            self.extended_calc.hamiltonian.update(self.extended_calc.density)
        else:
            self.update_hamiltonian()
        
        self.timer.start('record')        
        if self.foot_print:
            self.analysor.save_ele_step(self)
        self.timer.stop('record')
        self.record_time_cost = self.timer.timers['HamMM', 'record']
        
        self.timer.start('project hamiltonian')
        h_spkmm, s_pkmm = self.get_hs(self.extended_calc)

        if self.gate_mode == 'VM':
            ind = get_matrix_index(self.gate_basis_index)
            h_spkmm[:, :, ind.T, ind] += self.gate * s_pkmm[:, ind.T, ind]        

        self.timer.stop('project hamiltonian')                  
       
        for q in range(self.my_npk):
            if self.optimize:
                self.hsd.reset(0, q, s_pkmm[q], 'S')
            for s in range(self.my_nspins):
                self.hsd.reset(s, q, h_spkmm[s, q], 'H')
        #if self.ground:
        #    self.boundary_align_up()
        #    self.text('align shift-----' + str(self.align_shift))
  
    def get_forces(self, atoms):
        self.log('get_forces()')
        if self.non_sc:
            if not hasattr(self, 'contour'):
                self.contour = Contour(self.occupations.width * Hartree,
                            self.lead_fermi, self.bias, comm=self.wfs.gd.comm,
                             plot_eta=self.plot_eta, eta=self.eta,
                             plot_energy_range=self.plot_energy_range,
                             plot_energy_point_num=self.plot_energy_point_num)            
            if not hasattr(self, 'analysor'):
                self.analysor = Transport_Analysor(self, True)            
            if self.F_av is None:
                self.equivalent_atoms = self.atoms.copy()
                kwargs = self.gpw_kwargs.copy()
                kwargs['poissonsolver'] = PoissonSolver(nn=2)
                kpts = kwargs['kpts']
                kpts = kpts[:2] + (1,)
                kwargs['kpts'] = kpts
                if self.spinpol:
                    kwargs['mixer'] = MixerDif(self.density.mixer.beta, 5, weight=100.0)
                else:
                    kwargs['mixer'] = Mixer(self.density.mixer.beta, 5, weight=100.0)
                if 'txt' in kwargs and kwargs['txt'] != '-':
                    kwargs['txt'] = 'guess_' + kwargs['txt']            
                self.equivalent_atoms.set_calculator(gpaw.GPAW(**kwargs))
                calc = self.equivalent_atoms.calc
                calc.initialize(self.equivalent_atoms)
                calc.set_positions(self.equivalent_atoms)
                self.F_av = calc.get_forces(self.equivalent_atoms)                

            elif (atoms.positions != self.atoms.positions).any():
                self.atoms.set_positions(atoms.positions)
                self.equivalent_atoms.set_positions(atoms.positions)
                calc = self.equivalent_atoms.calc
                calc.density.reset()
                calc.set_positions(atoms)
                self.F_av = calc.get_forces(atoms)
            else:
                calc = self.equivalent_atoms.calc                
            self.extended_calc.hamiltonian = calc.hamiltonian
            self.analysor.save_bias_step(self)    
            self.analysor.save_ion_step(self)                
            return self.F_av
        
        else:            
            if (atoms.positions != self.atoms.positions).any():
                self.scf.converged = False
            if  hasattr(self.scf, 'converged') and self.scf.converged:
                pass
            else:
                self.negf_prepare(atoms)
                if np.sum(np.abs(self.bias)) < 1e-3:
                    self.ground = True
                self.get_selfconsistent_hamiltonian()
                self.analysor.save_ion_step(self)
                self.text('--------------ionic_step---' +
                          str(self.analysor.n_ion_step) + '---------------')
                self.F_av = None

            if not self.optimize:
                self.optimize = True

            if not self.use_qzk_boundary:    
                f = self.calculate_force()
                return f * Hartree / Bohr
            else:
                f = self.extended_calc.get_forces()[:len(self.atoms)]
                return f  

    def calculate_force(self):
        """Return the atomic forces."""
        self.log('calculate_force()')
        if self.F_av is not None:
            return self.F_av[:len(self.atoms)]
        natoms = len(self.extended_calc.wfs.setups)
        self.F_av = np.zeros((natoms, 3))

        hamiltonian = self.extended_calc.hamiltonian
        vt_sG = hamiltonian.vt_sG
        if len(vt_sG) == 2:
            vt_G = 0.5 * (vt_sG[0] + vt_sG[1])
        else:
            vt_G = vt_sG[0]

        wfs = self.extended_calc.wfs
        # Force from projector functions (and basis set):
        self.extended_calc.wfs.calculate_forces(hamiltonian, self.F_av)

        nn = self.surround.nn * 2
        vHt_g = self.surround.uncapsule(self, nn, hamiltonian.vHt_g,
                                                    self.finegd1, self.finegd)
        vt_G0 = self.surround.uncapsule(self, nn / 2, vt_G, self.gd1, self.gd)  
        if wfs.band_comm.rank == 0 and wfs.kd.comm.rank == 0:
            # Force from compensation charges:
            dF_aLv = self.density.ghat.dict(derivative=True)

            self.density.ghat.derivative(vHt_g, dF_aLv)
            for a, dF_Lv in dF_aLv.items():
                self.F_av[a] += np.dot(self.density.Q_aL[a], dF_Lv)
    
            # Force from smooth core charge:
            dF_av = self.density.nct.dict(derivative=True)
            self.density.nct.derivative(vt_G0, dF_av)
            for a, dF_v in dF_av.items():
                self.F_av[a] += dF_v[0]

            # Force from zero potential:
            dF_av = self.hamiltonian.vbar.dict(derivative=True)
            self.hamiltonian.vbar.derivative(self.density.nt_g, dF_av)
            for a, dF_v in dF_av.items():
                self.F_av[a] += dF_v[0]
    
            wfs.gd.comm.sum(self.F_av, 0)

        wfs.world.broadcast(self.F_av, 0)
        # Add non-local contributions:
        hamiltonian.xc.add_forces(self.F_av)
    
        if wfs.kd.symmetry:
            self.F_av = wfs.kd.symmetry.symmetrize_forces(self.F_av)

        self.forces.F_av = self.F_av[:len(self.atoms)]
        self.print_forces()
        return self.F_av[:len(self.atoms)]

    def calculate_to_bias(self, v_limit, num_v, gate=0, num_g=3, start=0):
        self.log('calculate_to_bias()')
        bias = np.linspace(0, v_limit, num_v)
        self.negf_prepare()
        # THa: initialize bias correct
        self.bias = [bias[0]/2., -bias[0]/2.]
        if abs(gate) > 0.001:
            gate = np.linspace(0, gate, num_g)
            for i in range(start, num_g):
                self.gate = gate[i]
                self.get_selfconsistent_hamiltonian()
            start = 0
        self.n_bias_step = start    
        for i in range(start, num_v):
            v = bias[i]
            self.bias = [v/2., -v /2.]
            self.get_selfconsistent_hamiltonian()        
    
    def calculate_to_gate(self, v_limit, num_v):
        self.log('calculate_to_gate()')
        gate = np.linspace(0, v_limit, num_v)
        self.negf_prepare() 
        for i in range(num_v):
            self.gate = gate[i]
            self.get_selfconsistent_hamiltonian()
            self.ground = True
        
    def get_potential_energy(self, atoms=None, force_consistent=False):
        self.log('get_potential_energy()')
        if self.non_sc:
            return self.equivalent_atoms.get_potential_energy()
        else:
            if hasattr(self.scf, 'converged') and self.scf.converged:
                pass
            else:
                self.negf_prepare()
                self.get_selfconsistent_hamiltonian()
            if force_consistent:
                # Free energy:
                return Hartree * self.hamiltonian.Etot
            else:
                # Energy extrapolated to zero Kelvin:
                return Hartree * (self.hamiltonian.Etot + 0.5 * self.hamiltonian.S)
      
    def induce_density_perturbation(self):
        self.log('induce_density_perturbation()')
        wfs = self.extended_calc.wfs
        basis_functions = wfs.basis_functions        
        density = self.density
        f_sM = np.empty((self.nspins, basis_functions.Mmax))
        D_asp = {}
        f_asi = {}
        for a in basis_functions.atom_indices:
            if a in self.perturbation_atoms:
                c = self.perturbation_charge / len(wfs.setups)  # distribute on all atoms
                f_si = wfs.setups[a].calculate_initial_occupation_numbers(
                        self.perturbation_magmom[a], True, charge=c, nspins=self.nspins)
                if a in basis_functions.my_atom_indices:
                    D_asp[a] = wfs.setups[a].initialize_density_matrix(f_si)
                f_asi[a] = f_si
            else:
                c = self.perturbation_charge / len(wfs.setups)
                f_si = wfs.setups[a].calculate_initial_occupation_numbers(
                        self.perturbation_magmom[a], True, charge=c, nspins=self.nspins)
                if a in basis_functions.my_atom_indices:
                    D_asp[a] = np.zeros_like(self.extended_D_asp[a])
                f_asi[a] = np.zeros_like(f_si)

        nt_sG = self.gd1.zeros(self.nspins)
        basis_functions.add_to_density(nt_sG, f_asi)
        nn = self.surround.nn
        self.perturbation_nt_sG = self.surround.uncapsule(self, nn, nt_sG, self.gd1,
                                                    self.gd)
        all_D_asp = collect_atomic_matrices(D_asp, wfs.setups, self.nspins,
                                            self.gd.comm, wfs.rank_a)
        D_asp = all_D_asp[:len(self.atoms)]
        self.perturbation_D_asp = D_asp
        #comp_charge = density.calculate_multipole_moments()
        #density.mix(comp_charge)
  
    def update_density(self):
        self.log('update_density()')
        self.timer.start('dmm recover')
        #self.fill_guess_with_leads()
        for kpt in self.extended_calc.wfs.kpt_u:
            if self.my_nspins == 2:
                kpt.rho_MM = self.hsd.D[kpt.s][kpt.q].recover(True)
            else:
                kpt.rho_MM = self.hsd.D[0][kpt.q].recover(True)
        self.timer.stop('dmm recover')        
        
        if not self.use_qzk_boundary:        
            density = self.density
        else:
            density = self.extended_calc.density
        self.timer.start('construct density')

        density.charge_eps = 1000
            
        nt_sG = self.gd1.zeros(self.nspins)
        self.extended_calc.wfs.calculate_density_contribution(nt_sG)
        
        if not self.use_qzk_boundary:
            nn = self.surround.nn
            density.nt_sG = self.surround.uncapsule(self, nn, nt_sG, self.gd1,
                                                    self.gd)
        else:
            density.nt_sG = nt_sG
            
        density.nt_sG += density.nct_G


        self.timer.stop('construct density')
        self.timer.start('atomic density')
        
        if not self.use_qzk_boundary:
            D_asp = self.extended_D_asp
        else:
            D_asp = self.extended_calc.density.D_asp
        self.extended_calc.wfs.calculate_atomic_density_matrices(D_asp)
        #all_D_asp = collect_D_asp2(D_asp, self.extended_calc.wfs.setups, self.nspins,
        #                    self.gd.comm, self.extended_calc.wfs.rank_a)
        
        if not self.use_qzk_boundary:    
            wfs = self.extended_calc.wfs
            all_D_asp = collect_atomic_matrices(D_asp, wfs.setups, self.nspins,
                                            self.gd.comm, wfs.rank_a)
            D_asp = all_D_asp[:len(self.atoms)]

            #distribute_D_asp(D_asp, density)
            distribute_atomic_matrices(D_asp, density.D_asp, density.setups)

        self.timer.stop('atomic density')
        comp_charge = density.calculate_multipole_moments()
        if self.neutral:
            self.normalize(comp_charge)
        density.mix(comp_charge)

        if self.step == 0 and self.analysor.n_bias_step in self.perturbation_steps:
            density.nt_sG += self.perturbation_nt_sG
            for a in range(len(density.setups)):
                if density.D_asp.get(a) is not None:
                    density.D_asp[a] += self.perturbation_D_asp[a]


    def normalize(self, comp_charge):
        self.log('normalize()')
        if not self.use_qzk_boundary:
            density = self.density
        else:
            density = self.extended_calc.density
        """Normalize pseudo density."""
        pseudo_charge = density.gd.integrate(density.nt_sG).sum()
        if pseudo_charge != 0:
            x = -(density.charge + comp_charge) / pseudo_charge
            density.nt_sG *= x - (x - 1) * self.alpha
            self.text('density scaling', x)        
           
    def update_hamiltonian(self):
        self.log('update_hamiltonian()')
        # only used by fixed bc
        
        ham = self.extended_calc.hamiltonian
        density = self.density
        self.timer.start('Hamiltonian')
        if ham.vt_sg is None:
            ham.vt_sg = ham.finegd.empty(ham.nspins)
            ham.vHt_g = ham.finegd.zeros()
            ham.vt_sG = ham.gd.empty(ham.nspins)
            #self.inner_poisson.initialize()
 
        nn = self.surround.nn * 2
        nt_sg = self.surround.capsule(self, nn, density.nt_sg, self.surround.nt_sg,
                                                    self.finegd1, self.finegd)
  
        Ebar = ham.finegd.integrate(ham.vbar_g, np.sum(nt_sg, axis=0),
                                     global_integral=False) 
        vt_g = ham.vt_sg[0]
        vt_g[:] = ham.vbar_g
        Eext = 0.0

        if ham.vext is not None:
            vt_g += ham.vext.get_potential(ham.finegd)
            Eext = np.vdot(vt_g, np.sum(nt_sg, axis=0)) * ham.finegd.dv - Ebar

        if ham.nspins == 2:
            ham.vt_sg[1] = vt_g
       
        Exc = ham.xc.calculate(ham.finegd, nt_sg, ham.vt_sg)
        Exc /= ham.gd.comm.size

        self.timer.start('Poisson')

        if self.hamiltonian.vHt_g is None:
            self.hamiltonian.vHt_g = self.finegd.zeros()

        actual_charge = self.finegd.integrate(density.rhot_g)
        self.text('actual_charge' + str(actual_charge))

        ham.npoisson = self.inner_poisson.solve(self.hamiltonian.vHt_g,
                                                  density.rhot_g,
                                                  charge = self.total_charge)
                                                  #charge=-density.charge)
        if self.fixed and self.gate_mode == 'VG':
            if self.gate_fun is None:
                self.hamiltonian.vHt_g += self.gate
            else:
                from scipy import interpolate                
                gate_vg = self.finegd.zeros(global_array=True)
                nz = self.gate_fun.shape[0]
                nzz = gate_vg.shape[2]
                xxx = np.linspace(0, 1, nzz)
                xx = np.linspace(0, 1, nz)
                f = interpolate.interp1d(xx, self.gate_fun)
                yyy = f(xxx)
                for i in range(nzz):
                    gate_vg[:, :, i] = yyy[i]
                local_gate_vg = self.finegd.zeros()
                self.finegd.distribute(gate_vg, local_gate_vg)
                self.hamiltonian.vHt_g += self.gate * local_gate_vg / Hartree
               
        self.surround.combine_vHt_g(self, self.hamiltonian.vHt_g)
        self.text('poisson iterations :' + str(ham.npoisson))
        self.timer.stop('Poisson')
        Epot = 0.5 * self.hamiltonian.finegd.integrate(self.hamiltonian.vHt_g,
                                                       density.rhot_g,
                                                        global_integral=False)
        Ekin = 0.0
        for vt_g, vt_G in zip(ham.vt_sg, ham.vt_sG):
            vt_g += ham.vHt_g
            ham.restrict(vt_g, vt_G)
        self.surround.refresh_vt_sG(self)
        
        nn = self.surround.nn
        vt_sG = self.surround.uncapsule(self, nn, ham.vt_sG, self.gd1, self.gd)
        for vt_G, nt_G in zip(vt_sG, density.nt_sG):    
            Ekin -= self.gd.integrate(vt_G, nt_G - density.nct_G,
                                                       global_integral=False)            
        self.timer.start('atomic hamiltonian')
        
        ham = self.hamiltonian
        W_aL = {}
        for a in density.D_asp:
            W_aL[a] = np.empty((ham.setups[a].lmax + 1)**2)
        density.ghat.integrate(ham.vHt_g, W_aL)
        ham.dH_asp = {}
        for a, D_sp in density.D_asp.items():
            W_L = W_aL[a]
            setup = ham.setups[a]

            D_p = D_sp.sum(0)
            dH_p = (setup.K_p + setup.M_p +
                    setup.MB_p + 2.0 * np.dot(setup.M_pp, D_p) +
                    np.dot(setup.Delta_pL, W_L))
            Ekin += np.dot(setup.K_p, D_p) + setup.Kc
            Ebar += setup.MB + np.dot(setup.MB_p, D_p)
            Epot += setup.M + np.dot(D_p, (setup.M_p + np.dot(setup.M_pp, D_p)))

            ham.dH_asp[a] = dH_sp = np.zeros_like(D_sp)
            Exc += ham.xc.calculate_paw_correction(setup, D_sp, dH_sp, a=a)

            if setup.HubU is not None:
                nspins = len(D_sp)
                
                l_j = setup.l_j
                l   = setup.Hubl
                nl  = np.where(np.equal(l_j,l))[0]
                nn  = (2*np.array(l_j)+1)[0:nl[0]].sum()
                
                for D_p, H_p in zip(D_sp, ham.dH_asp[a]):
                    [N_mm,V] = ham.aoom(unpack2(D_p),a,l)
                    N_mm = N_mm / 2 * nspins
                     
                    Eorb = setup.HubU / 2. * (N_mm - np.dot(N_mm,N_mm)).trace()
                    Vorb = setup.HubU * (0.5 * np.eye(2*l+1) - N_mm)
                    Exc += Eorb
                    if nspins == 1:
                        # add contribution of other spin manyfold
                        Exc += Eorb
                    
                    if len(nl)==2:
                        mm  = (2*np.array(l_j)+1)[0:nl[1]].sum()
                        
                        V[nn:nn+2*l+1,nn:nn+2*l+1] *= Vorb
                        V[mm:mm+2*l+1,nn:nn+2*l+1] *= Vorb
                        V[nn:nn+2*l+1,mm:mm+2*l+1] *= Vorb
                        V[mm:mm+2*l+1,mm:mm+2*l+1] *= Vorb
                    else:
                        V[nn:nn+2*l+1,nn:nn+2*l+1] *= Vorb
                    
                    Htemp = unpack(H_p)
                    Htemp += V
                    H_p[:] = pack2(Htemp)

            dH_sp += dH_p
            Ekin -= (D_sp * dH_sp).sum()

        self.timer.stop('atomic hamiltonian')

        ham.Enlxc = 0.0#xcfunc.get_non_local_energy()
        ham.Enlkin = ham.xc.get_kinetic_energy_correction()
        if ham.Enlxc != 0 or ham.Enlkin != 0:
            print('Where should we do comm.sum() ?')
        
        comm = ham.gd.comm
        ham.Ekin0 = comm.sum(Ekin)
        ham.Epot = comm.sum(Epot)
        ham.Ebar = comm.sum(Ebar)
        ham.Eext = comm.sum(Eext)
        ham.Exc = comm.sum(Exc)
        
        ham.Exc += ham.Enlxc
        ham.Ekin0 += ham.Enlkin
        
        #dH_asp = collect_D_asp3(ham, self.density.rank_a)
        dH_asp = collect_atomic_matrices(ham.dH_asp, ham.setups, ham.nspins,
                                         ham.gd.comm, self.density.rank_a)
        self.log('  collect_atomic_matrices()')
        self.surround.combine_dH_asp(self, dH_asp)
        if not self.analysis_mode:
            ham.get_energy(self.occupations)
        self.timer.stop('Hamiltonian')
        self.log('finish update_hamiltonian()')
        

    
    def print_boundary_charge(self):
        self.log('print_boundary_charge()')
        boundary_charge = []
        print_info = ''
        if self.hsd.S[0].extended:
            n = -2
        else:
            n = -1
        for i in range(self.lead_num):
            nb = self.nblead[i]
            qr_mm = np.zeros([nb, nb])
            for s in range(self.my_nspins):
                for pk in range(self.my_npk):
                    D = self.hsd.D[s][pk]
                    S = self.hsd.S[pk]
                    qr_mm+= np.real(dot(D.diag_h[i][n].recover(),
                                                     S.diag_h[i][n].recover()))
                    qr_mm += np.real(dot(D.dwnc_h[i][n], S.upc_h[i][n]))
                    if S.extended:
                        qr_mm += np.real(dot(D.dwnc_h[i][n + 1], S.upc_h[i][n + 1]))
                    else:
                        qr_mm += np.real(dot(D.dwnc_h[i][n], S.upc_h[i][n]))
            self.wfs.kd.comm.sum(qr_mm)
            boundary_charge.append(np.real(np.trace(qr_mm)))
            if i != 0:
                print_info += '******'
            print_info += str(boundary_charge[i])
        self.text(print_info)
    
    def set_buffer(self):
        self.log('set_buffer()')
        self.nbmol_inner = self.nbmol 
        if self.use_lead:
            self.nbmol_inner -= np.sum(self.buffer)
        ind = np.arange(self.nbmol)
        buffer_ind = []
        lead_ind = []

        for i in range(self.lead_num):
            buffer_ind += list(self.buffer_index[i])
            lead_ind += list(self.lead_index[i])

        ind = np.delete(ind, buffer_ind)
        self.inner_mol_index = ind
        #self.gate_mol_index = np.delete(ind, lead_ind)
        
        for i in range(self.lead_num):
             self.inner_lead_index[i] = np.searchsorted(ind,
                                                           self.lead_index[i])

    def integral_diff_weight(self, denocc, denvir, method='transiesta'):
        self.log('integral_diff_weight()')
        if method=='transiesta':
            eta = 1e-16
            weight = denocc * denocc.conj() / (denocc * denocc.conj() +
                                               denvir * denvir.conj() + eta)
        return weight

    def fill_lead_with_scat(self):
        self.log('fill_lead_with_scat()')
        assert self.hsd.extended
        n = -1
        m = -2
        for  i in range(self.lead_num):
            for s in range(self.my_nspins):
                for pk in range(self.my_npk):
                    self.lead_hsd[i].reset(s, pk,
                        self.hsd.H[s][pk].diag_h[i][m].recover(), 'H')
                    #if i == 0:
                    #    self.lead_couple_hsd[i].reset(s, pk,
                    #        self.hsd.H[s][pk].upc_h[i][n], 'H')
                    #elif i == 1:
                    #    self.lead_couple_hsd[i].reset(s, pk,
                    #        self.hsd.H[s][pk].dwnc_h[i][n], 'H')                        
                    #else:
                    #    raise NotImplementError()
       
    def estimate_transport_matrix_memory(self):
        self.log('estimate_transport_matrix_memory()')
        sum = 0
        ns = self.wfs.nspins
        if self.use_lead:
            nk = len(self.ibzk_qc_lead[0])
            nb = max(self.nblead)
            npk = len(self.wfs.kd.ibzk_qc)
            unit_real = np.array(1,float).itemsize
            unit_complex = np.array(1, complex).itemsize
            
            gamma = len(self.wfs.kd.bzk_kc) == 1 and not self.wfs.kd.bzk_kc[0].any()          
            if gamma:
                unit = unit_real
            else:
                unit = unit_complex
            sum += self.lead_num * (2 * ns + 1)* npk * nb**2 * unit
            
            if self.LR_leads:
                sum += ( 2 * ns + 1) * npk * nb ** 2 * unit
            sum += ns * npk * nb**2 * unit
            #print 'lead matrix memery  MB',  sum *1e-6
           
            ntgt = 200 // self.wfs.gd.comm.size
            tmp = self.lead_num * ns * npk * ntgt * nb**2 * unit_complex
            sum += tmp
            #print 'selfenergy memery  MB',  tmp *1e-6

        if gamma:
            unit = unit_real
        else:
            unit = unit_complex
        nk = len(self.wfs.kd.ibzk_qc)
        nb = self.wfs.setups.nao
        sum += (2*ns + 1) * nk * nb**2 * unit
        return tmp, (sum - tmp)
           
    def reset_lead_hs(self, s, k):
        self.log('reset_lead_hs()')
        if self.use_lead:    
            sg = self.selfenergies
            for i in range(self.lead_num):
                sg[i].s = s
                sg[i].pk = k

    def initialize_green_function(self):
        self.log('initialize_green_function()')
        self.selfenergies = []
        if self.use_lead:
            directions = ['left', 'right'] + ['left'] * 1000
            for i in range(self.lead_num):
                self.selfenergies.append(LeadSelfEnergy(self.lead_hsd[i],
                                            self.lead_couple_hsd[i],
                                           self.se_data_path, directions[i]))
                if not self.use_qzk_boundary:
                    self.selfenergies[i].set_bias(self.bias[i])
 
    def calculate_iv(self, v_limit=3, num_v=16, start=0):
        self.log('calculate_iv()')
        if self.non_sc:
            self.negf_prepare()
            self.non_sc_analysis()
        else:
            self.calculate_to_bias(v_limit, num_v, start=start)

    def calculate_transmission(self):
        self.log('calculate_transmission()')
        self.negf_prepare()
        self.non_sc_analysis()          
        
    def recover_kpts(self, calc):
        self.log('recover_kpts()')
        wfs = calc.wfs
        wfs.eigensolver.iterate(calc.hamiltonian, wfs)
        calc.occupations.calculate(wfs)

    def estimate_memory(self, mem):
        """Estimate memory use of this object."""
        self.log('estimate_memory()')
        mem_init = maxrss() # XXX initial overhead includes part of Hamiltonian
        mem.subnode('Initial overhead', mem_init)
        for name, obj in [('Density', self.density),
                          ('Hamiltonian', self.hamiltonian),
                          ('Wavefunctions', self.wfs)]:
            obj.estimate_memory(mem.subnode(name))
        #for i in range(self.lead_num):
        #    atoms = self.get_lead_atoms(i)
        #    atoms.calc.estimate_memory(mem.subnode('Leads' + str(i), 0))
        se_mem, mat_mem = self.estimate_transport_matrix_memory()
        mem.subnode('Matrix', mat_mem)
        mem.subnode('Selfenergy', se_mem)

    def get_extended_atoms(self):
        self.log('get_extended_atoms()')
        # for LR leads only
        if self.extended_atoms is None or (self.extended_atoms is not None
                                           and self.optimize):
            atoms = self.atoms.copy()
            cell = np.diag(atoms.cell)
            ex_cell = cell.copy()
            di = 2
            for i in range(self.lead_num):
                if self.leads is None:
                    atoms_l = self.atoms[self.pl_atoms[i]].copy()
                else:
                    atoms_l = self.leads[i].copy()
                    if i == 0:
                        j = 0
                    else:
                        j = 1
                    atoms_l.positions += self.atoms[self.pl_atoms[i]].positions[j]- \
                                              self.leads[i].positions[j]                        
                    
                cell_l = self.pl_cells[i]
                assert self.gd.orthogonal
                ex_cell[di] += self.gd.h_cv[2, 2] * Bohr * self.bnc[i]
                for atom in atoms_l:
                    if i == 0:
                        atom.position[di] -= cell_l[di]
                    else:
                        atom.position[di] += cell_l[di]
                atoms += atoms_l
            atoms.set_cell(ex_cell)
            atoms.set_pbc(self.atoms._pbc)
            atoms.positions[:, 2] += self.gd.h_cv[2, 2] * Bohr * self.bnc[0]
            self.extended_atoms = atoms
       
        if not self.optimize:
            p = self.gpw_kwargs.copy()
            p['h'] = None
            N_c = self.gd.N_c.copy()
            for i in range(self.lead_num):
                N_c[2] += self.bnc[i]
            p['gpts'] = N_c
            if 'mixer' in p:
                if hasattr(self.density.mixer, 'mixers'):
                    p['mixer'] = Mixer(self.density.mixer.beta, 5, weight=100.0)
                else:
                    p['mixer'] = MixerDif(self.density.mixer.beta, 5, weight=100.0)
            p['poissonsolver'] = PoissonSolver(nn=2)        
            if type(p['basis']) is dict and len(p['basis']) == len(self.atoms):
                basis = {}
                for i, btype in enumerate(range(len(self.atoms))):
                    basis[i] = p['basis'][i]
                for i, btype in enumerate(self.pl_atoms[0]):
                    basis[i + len(self.atoms)] = p['basis'][i]
                for i, btype in enumerate(self.pl_atoms[1]):
                    basis[i + len(self.atoms) + len(self.pl_atoms[0])] = p['basis'][i]
                p['basis'] = basis
            self.extended_atoms.set_calculator(Lead_Calc(**p))

    def get_linear_hartree_potential(self):
        self.log('get_linear_hartree_potential()')
        global_linear_vHt = self.finegd.zeros(global_array=True)
        dim = self.finegd.N_c[2]
        vt = np.linspace(self.bias[0]/Hartree, self.bias[1]/Hartree, dim)
        for i in range(dim):
            global_linear_vHt[:, :, i] = vt[i]
        self.linear_vHt_g = self.finegd.zeros()
        self.finegd.distribute(global_linear_vHt, self.linear_vHt_g)

    def get_inner_setups(self):
        self.log('get_inner_setups()')
        spos_ac0 = self.atoms.get_scaled_positions() % 1.0
        self.wfs.set_positions(spos_ac0)
        if self.hubbard_parameters is not None:
            element, U_ev, scale, store = self.hubbard_parameters
            for i, atom in enumerate(self.atoms):
                if atom.symbol == element:
                    self.hamiltonian.setups[i].set_hubbard_u(U_ev/Hartree,2,scale,store)
        self.inner_setups = self.wfs.setups
        self.inner_atom_indices = self.wfs.basis_functions.atom_indices
        self.inner_my_atom_indices = self.wfs.basis_functions.my_atom_indices
        self.inner_rank_a = self.wfs.rank_a
        
    def set_extended_positions(self):
        self.log('set_extended_positions()')
        spos_ac0 = self.atoms.get_scaled_positions() % 1.0
        spos_ac = self.extended_atoms.get_scaled_positions() % 1.0
        self.wfs.set_positions(spos_ac0)
        old_extended_rank_a = self.extended_calc.wfs.rank_a
        self.extended_calc.wfs.set_positions(spos_ac)

        self.density.set_positions(spos_ac0, self.wfs.rank_a)
        self.hamiltonian.set_positions(spos_ac0, self.wfs.rank_a)
        self.extended_calc.hamiltonian.set_positions(spos_ac,
                                                self.extended_calc.wfs.rank_a)

        if self.extended_D_asp is not None:
            requests = []
            D_asp = {}
            for a in self.extended_calc.wfs.basis_functions.my_atom_indices:
                if a in self.extended_D_asp:
                    D_asp[a] = self.extended_D_asp.pop(a)
                else:
                    # Get matrix from old domain:
                    ni = self.extended_calc.wfs.setups[a].ni
                    D_sp = np.empty((self.nspins, ni * (ni + 1) // 2))
                    D_asp[a] = D_sp
                    requests.append(self.gd.comm.receive(D_sp, old_extended_rank_a[a],
                                                         tag=a, block=False))
                
            for a, D_sp in self.extended_D_asp.items():
                # Send matrix to new domain:
                requests.append(self.gd.comm.send(D_sp, self.extended_calc.wfs.rank_a[a],
                                                  tag=a, block=False))
            for request in requests:
                self.gd.comm.wait(request)
            self.extended_D_asp = D_asp

        density = self.density
        wfs = self.extended_calc.wfs
        gd = self.gd
        gd1 = self.extended_calc.gd

        magmom_a = self.extended_atoms.get_initial_magnetic_moments()
        if density.nt_sG is None:
            if wfs.kpt_u[0].f_n is None or wfs.kpt_u[0].C_nM is None:
                f_sM = np.empty((self.nspins, wfs.basis_functions.Mmax))
                self.extended_D_asp = {}
                density.D_asp = {}
                f_asi = {}
                c = density.charge / len(density.setups)  
                for a in wfs.basis_functions.atom_indices:
                    f_si = wfs.setups[a].calculate_initial_occupation_numbers(
                        magmom_a[a], density.hund, charge=c,
                        nspins=self.nspins)
                    
                    if a in wfs.basis_functions.my_atom_indices:
                        self.extended_D_asp[a] = (
                            wfs.setups[a].initialize_density_matrix(f_si))
                    f_asi[a] = f_si

                for a in self.wfs.basis_functions.atom_indices:
                    setup = self.wfs.setups[a]
                    f_si = setup.calculate_initial_occupation_numbers(
                        density.magmom_av[a, 2], density.hund, charge=c,
                        nspins=self.nspins)
                    if a in self.wfs.basis_functions.my_atom_indices:
                        density.D_asp[a] = setup.initialize_density_matrix(
                            f_si)                    

                all_D_asp = []
                for a, setup in enumerate(wfs.setups):
                    D_sp = self.extended_D_asp.get(a)
                    if D_sp is None:
                        ni = setup.ni
                        D_sp = np.empty((density.nspins, ni * (ni + 1) // 2))
                    if density.gd.comm.size > 1:
                        density.gd.comm.broadcast(D_sp,
                                     self.extended_calc.hamiltonian.rank_a[a])
                    all_D_asp.append(D_sp)      
                
                D_asp = all_D_asp[:len(self.atoms)]
                #distribute_D_asp(D_asp, density)
                distribute_atomic_matrices(D_asp, density.D_asp,
                                           density.setups)

                nt_sG = gd1.zeros(self.nspins)
                wfs.basis_functions.add_to_density(nt_sG, f_asi)
                nn = self.surround.nn
                density.nt_sG = self.surround.uncapsule(self, nn, nt_sG, gd1, gd)
                density.nt_sG += density.nct_G
                density.normalize()

            else:
                density.nt_sG = self.gd.empty(self.nspins)
                density.calculate_pseudo_density(wfs)
                density.nt_sG += density.nct_G
                density.normalize() 
  
        comp_charge = density.calculate_multipole_moments()
        density.interpolate_pseudo_density(comp_charge)
        density.calculate_pseudo_charge()
        
        self.update_hamiltonian()
        self.log('  update_hamiltonian() called')
        self.scf.reset()
        self.log('  scf.reset() called')
        self.forces.reset()
        self.log('  forces.reset() called')
        self.print_positions()
        self.log('  print_positions() called')
        

    def analysis(self, n, n1=0, gate=False, gate_uplimit=None):
        self.log('analysis()')
        self.guess_steps = 1
        self.negf_prepare()
        flag = True
        self.contour = Contour(self.occupations.width * Hartree,
                            self.lead_fermi, self.bias, comm=self.wfs.gd.comm,
                             plot_eta=self.plot_eta, eta=self.eta,
                             plot_energy_range=self.plot_energy_range,
                             plot_energy_point_num=self.plot_energy_point_num)
        if not hasattr(self, 'analysor'):
            self.analysor = Transport_Analysor(self, True)
            
        for i in range(n1, n):
            if gate :
                self.gate = np.linspace(0, gate_uplimit, n)[i]
            else:
                self.gate = 0
            if i > n1:
                flag = False
            self.analysor.n_bias_step = i
            fd = file('bias_data' + str(i + 1), 'r')
            self.bias, vt_sG, dH_asp = cPickle.load(fd)
            fd.close()
          
            for j in range(self.lead_num):
                self.analysor.selfenergies[j].set_bias(self.bias[j])
            self.surround.combine_dH_asp(self, dH_asp)
            self.log('vt_sG: {0}'.format(vt_sG.shape))
            self.log('extcalc.ham.vt_sG: {0}'.format(self.extended_calc.hamiltonian.vt_sG))
            self.gd1.distribute(vt_sG, self.extended_calc.hamiltonian.vt_sG) 
            h_spkmm, s_pkmm = self.get_hs(self.extended_calc)
            if self.gate_mode == 'VM':
                ind = get_matrix_index(self.gate_basis_index)
                h_spkmm[:, :, ind.T, ind] += self.gate * s_pkmm[:, ind.T, ind]   
            
            nb = s_pkmm.shape[-1]
            dtype = s_pkmm.dtype
            for q in range(self.my_npk):
                self.hsd.reset(0, q, s_pkmm[q], 'S', flag)                
                for s in range(self.my_nspins):
                    self.hsd.reset(s, q, h_spkmm[s, q], 'H', flag)
                    self.hsd.reset(s, q, np.zeros([nb, nb], dtype), 'D', flag)
            if flag:
                self.log('  analysis() call append_buffer_hsd()')
                self.append_buffer_hsd()                    
                self.log('  analysis() append_buffer_hsd() finished')
            
            #self.analysor.save_ele_step(self)            
            self.log('  analysis() call save_bias_step()')
            self.analysor.save_bias_step(self)
            self.log('  analysis() save_bias_step() finished')

    def save_lead_hamiltonian_matrix(self):
        self.log('save_lead_hamiltonian_matrix()')
        print('assert self.nspins == 1')
        self.guess_steps = 1
        self.negf_prepare()
        kpt_comm = self.wfs.kd.comm
        for i in range(self.lead_num):
            nb = self.nblead[i]
            if kpt_comm.rank == 0:
                lead_overlap_matrix = np.zeros([self.npk, nb, nb], self.wfs.dtype)
                lead_coupling_overlap_matrix = np.zeros([self.npk, nb, nb], self.wfs.dtype)
                lead_hamiltonian_matrix = np.zeros([self.nspins, self.npk, nb, 
                                                              nb], self.wfs.dtype)
                lead_coupling_hamiltonian_matrix = np.zeros([self.nspins, self.npk, 
                                                          nb, nb], self.wfs.dtype)
            else:
                lead_overlap_matrix = None
                lead_coupling_overlap_matrix = None
                lead_hamiltonian_matrix = None
                lead_coupling_hamiltonian_matrix = None

            local_lead_overlap_matrix = np.zeros([self.my_npk, nb, nb], self.wfs.dtype)
            local_lead_coupling_overlap_matrix = np.zeros([self.my_npk, nb, nb], self.wfs.dtype)
            local_lead_hamiltonian_matrix = np.zeros([self.my_nspins, self.my_npk, nb, 
                                                              nb], self.wfs.dtype)
            local_lead_coupling_hamiltonian_matrix = np.zeros([self.my_nspins, self.my_npk, 
                                                          nb, nb], self.wfs.dtype)
            for pk in range(self.my_npk):
                local_lead_overlap_matrix[pk] = self.lead_hsd[i].S[pk].recover()
                local_lead_coupling_overlap_matrix[pk] = self.lead_couple_hsd[i].S[pk].recover()
                for s in range(self.my_nspins):
                    local_lead_hamiltonian_matrix[s, pk] = self.lead_hsd[i].H[s][pk].recover()
                    local_lead_coupling_hamiltonian_matrix[s, pk] = self.lead_couple_hsd[i].H[s][pk].recover()
            kpt_comm.gather(local_lead_overlap_matrix, 0, lead_overlap_matrix)
            kpt_comm.gather(local_lead_coupling_overlap_matrix, 0, lead_coupling_overlap_matrix)
            kpt_comm.gather(local_lead_hamiltonian_matrix, 0, lead_hamiltonian_matrix)
            kpt_comm.gather(local_lead_coupling_hamiltonian_matrix, 0, lead_coupling_hamiltonian_matrix)
            if world.rank == 0:
                fd = file('lead_hs_matrix' + str(i), 'wb')
                cPickle.dump((lead_overlap_matrix,
                             lead_coupling_overlap_matrix,
                             lead_hamiltonian_matrix,
                             lead_coupling_hamiltonian_matrix),
                             fd, 2)
                fd.close()

    def save_scat_hamiltonian_matrix(self,n, n1=0):
        self.log('save_scat_hamiltonian_matrix()')
        print('run save_lead_hamiltonian_matrix before')
        kpt_comm = self.wfs.kd.comm
        flag = True
        for i in range(n1, n):
           if i > n1:
               flag = False
           fd = file('bias_data' + str(i + 1), 'r')
           self.bias, vt_sG, dH_asp = cPickle.load(fd)
           fd.close()
           self.surround.combine_dH_asp(dH_asp)
           self.gd1.distribute(vt_sG, self.extended_calc.hamiltonian.vt_sG) 
           h_spkmm, s_pkmm = self.get_hs(self.extended_calc)
           if kpt_comm.rank == 0:
               scat_overlap_matrix = np.zeros([self.npk, self.nbmol, self.nbmol], self.wfs.dtype)
               scat_hamiltonian_matrix = np.zeros([self.nspins, self.npk, self.nbmol, self.nbmol], self.wfs.dtype)
           else:
               scat_overlap_matrix = None
               scat_hamiltonian_matrix = None
           local_scat_overlap_matrix = np.zeros([self.my_npk, self.nbmol, self.nbmol], self.wfs.dtype)
           local_scat_hamiltonian_matrix = np.zeros([self.my_nspins, self.my_npk, self.nbmol, self.nbmol], self.wfs.dtype)
           for q in range(self.my_npk):
               self.hsd.reset(0, q, s_pkmm[q], 'S', flag)                
               for s in range(self.my_nspins):
                   self.hsd.reset(s, q, h_spkmm[s, q], 'H', flag)
                   self.hsd.reset(s, q, np.zeros([self.nbmol, self.nbmol], self.wfs.dtype), 'D', flag)
           if flag:
                self.append_buffer_hsd()                    

           for q in range(self.my_npk):
               local_scat_overlap_matrix[q] = self.hsd.S[q].recover()
               for s in range(self.my_nspins):
                   local_scat_hamiltonian_matrix[s, q] = self.hsd.H[s][q].recover()
           kpt_comm.gather(local_scat_overlap_matrix, 0, scat_overlap_matrix)
           kpt_comm.gather(local_scat_hamiltonian_matrix, 0, scat_hamiltonian_matrix)
           if world.rank ==0:
               fd = file('scat_hs_matrix' + str(i), 'wb')
               cPickle.dump((scat_overlap_matrix,
                             scat_hamiltonian_matrix),
                             fd, 2)
               fd.close()
            
    def analysis_prepare(self, bias_step):
        self.log('analysis_prepare()')
        fd = file('lead_hs', 'r')
        lead_s00, lead_s01, lead_h00, lead_h01 = cPickle.load(fd)
        fd.close()
        plotter = Transport_Plotter('bias', 'bias_plot_data')
        s00 = plotter.bias_steps[bias_step].s00
        h00 = plotter.bias_steps[bias_step].h00
        lead_fermi = plotter.bias_steps[bias_step].lead_fermi    
        
        self.initialize_transport()
        flag = True
        if not hasattr(self, 'analysor'):
            self.analysor = Transport_Analysor(self, True)
        self.contour = Contour(self.occupations.width * Hartree,
                               self.lead_fermi, self.bias,
                               comm=self.wfs.gd.comm,
                               plot_eta=self.plot_eta, eta=self.eta,
                               plot_energy_range=self.plot_energy_range,
                             plot_energy_point_num=self.plot_energy_point_num)
        
        for j in range(self.lead_num):
                self.analysor.selfenergies[j].set_bias(self.bias[j])
        dtype = s00.dtype
        for q in range(self.npk):
            self.hsd.reset(0, q, s00[q], 'S', True)                
            for s in range(self.nspins):
                self.hsd.reset(s, q, h00[s, q], 'H', True)
                self.hsd.reset(s, q, np.zeros_like(h00[s, q]), 'D', True)            
        self.recover_lead_info(lead_s00, lead_s01, lead_h00, lead_h01, lead_fermi)      
        self.append_buffer_hsd()             

    def analysis_project_prepare(self):
        self.log('analysis_project_prepare()')
        spos_ac0 = self.atoms.get_scaled_positions() % 1.0
        spos_ac = self.extended_atoms.get_scaled_positions() % 1.0
        self.extended_calc.wfs.set_positions(spos_ac)        
        self.wfs.set_positions(spos_ac0)
        
             
