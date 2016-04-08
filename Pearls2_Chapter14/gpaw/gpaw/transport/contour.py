import numpy as np
from gpaw.mpi import world
from gpaw.transport.tools import fermidistribution, gather_ndarray_dict
from ase.parallel import parprint, paropen

'''The Contour class is used to get the integral contour parallelly and control
   the necessary parameters.
 
   To get a Keldysh Green Function, one need to to the integral like this:
 
                      /                          
                      |                         ---- 
      D=    1.j/(2*pi)| (Gr(Z) - Ga(Z)f(Z)dZ -  \   (Gr(E) - Ga(E)) * kt 
                      /                         /___
                   1, 2, 3                        4
                             /
                             |
                   + 1/(2*pi)| G<(E)dE
                             /
                             5
   The Path indexed by (1,2,3,4,5,6) looks like below:   
 
                              ^   
            |\                |
            | \               |
            |  \ 2            |
            |   \             |
          1 |    \            |
            |     \    3      |        6
            |      \---------------------
            |                 |4                 5
            |    --------------------------------------
           --------------------------------------------->
 
 
    Hierachy:
       NID(node ID) is the most important information of the energy node.
       With it you can get the corresponding energy and weight.
 
    To evaluate a function(named f) integral on a region[c-h, c+h],
    first, we use the Gauss-Lobatto formula to get a evaluation Q1
      / c+h
      |
      |  f(x) dx  = 1/6*f(c-h) + 5/6*f(c-h/sqrt(5)) + 5/6*f(c+h/sqrt(5)) + 1/6*f(c+h) = Q1
      |
      / c-h
 
    then, we use the cooresponding Kronrod formula to get another evalution Q2
      / c+h
      |
      |  f(x) dx  = 77/1470*f(c-h) + 432/1470*f(c-sqrt(2/3)*h) + 625/1470*f(c-h/sqrt(5)) +
      |
      / c-h
           672/1470*f(c) + 625/1470*f(c+h/sqrt(5)) + 432/1470*f(c+sqrt(2/3)*h) + 77/1470*f(c+h) = Q2

    we compare the two results, if the difference is less than tolerance, then divide
    the region into six parts, and for each of them, do the two evaluations until the difference
    less than the tolerance.''' 

class EnergyNode:
    def __init__(self, type, nlead):
        self.type = type
        self.num = 0
        self.lead_num = nlead
        self.energy = []
        self.weight = []
        self.nres = 0
        self.sigma = []
        for i in range(nlead):
            self.sigma.append([])
        if type == 'eq':
            self.fermi_factor = []
        elif type == 'ne':
            self.fermi_factor = []
            for i in range(nlead):
                self.fermi_factor.append([[], []])
        else:
            raise TypeError('unkown EnergyNode type')

    def add(self, elist, wlist, flist, siglist):
        self.num += len(elist)
        self.energy += elist
        self.weight += wlist
        if self.type == 'eq':
            self.fermi_factor += flist
        elif self.type == 'ne':
            for i in range(self.lead_num):
                for j in [0, 1]:
                    self.fermi_factor[i][j] += flist[i][j]
        else:
            raise TypeError('unkown EnergyNode type')
        for i in range(self.lead_num):
            self.sigma[i] += siglist[i]

    def set_nres(self, nres):
        self.nres = nres
        
class Path:
    poles_num = 4
    #int_step = 0.02
    bias_step = 0.1
    zone_sample = np.array([0, 0.55278640450004, 1.44721359549996, 2.0,
                            0.0]) / 2.
    zone_length = np.array([2., 0.55278640450004, 0.89442719099992007,
                                                     0.55278640450004]) / 2.0
    sample = np.array([0, 0.18350341907227,   0.55278640450004,   1.0,
         1.44721359549996,   1.81649658092773, 2.0]) / 2.
    zone_weight  = np.array([2., 0.55278640450004, 0.89442719099992007,
                                                     0.55278640450004]) / 2.0
    weights0 = np.array([6.0, 1.0, 0.0, 5.0, 0.0, 5.0, 0.0, 1.0]) / 6.0
    weights1 = np.array([1470., 77.0, 432.0, 625.0, 672.0, 625.0,
                          432.0, 77.0]) / 1470.
    weights2 = [3.0 / 8, 7.0 / 6, 23.0 / 24] 
    weights3 = [23.0 / 24, 7.0 / 6, 3.0 / 8]
    #these two parameters used to build selfenergy database are obselete
    bias_window_begin = -3
    bias_window_end = 3
    #------------
    def __init__(self, begin, end, index, maxdepth=7, type='Gaussian',
                 kt=0.1, int_step=0.02):
        self.begin = begin
        self.end = end
        self.type = type
        self.index = index
        self.maxdepth = maxdepth
        self.nids = []
        self.energies = []
        self.functions = []
        self.ses = []
        self.num = 10 ** self.maxdepth
        self.full_nodes = False
        self.kt = kt
        self.int_step = int_step
        self.initialize()

    def initialize(self):
        if self.type == 'linear':
            self.ne = int(max(np.ceil(np.abs((self.end - self.begin)
                                                  / self.int_step)) + 1, 6))
            self.int_step = np.abs(self.end - self.begin) / (self.ne - 1)
        
    def get_flags(self, num, path_flag=False):
        flags = []
        if self.type == 'Gaussian':
            if path_flag:
                assert num < self.num * 10
            else:
                assert num < self.num                
            digits = self.maxdepth
            
        elif self.type == 'linear':
            digits = int(np.ceil(np.log10(self.ne)))            
            
        elif self.type == 'poles':
            if self.full_nodes:
                nids_num = int((self.bias_window_end -
                                self.bias_window_begin) //
                                  self.bias_step + 1)                
                digits = int(np.ceil(np.log10(nids_num))) 
            elif path_flag:
                digits = 2
            else:
                digits = 3

        if path_flag:
            digits += 1

        for i in range(digits):
            unit = 10 ** (digits - i)
            digit = (num - (num // unit) * unit) // (unit // 10)
            flags.append(digit)            
        return np.array(flags)
        
    def get_full_nids(self):
        self.full_nodes = True
        self.nids = []
        self.energies = []
        depth = self.maxdepth    
        if self.type == 'Gaussian':
            base = self.num * self.index
            assert depth < 7
            dimension = [3] * (depth - 1) + [6]
            for i in range(self.num):
                flags = self.get_flags(i)
                this_node = True
                for j in range(depth):
                    if flags[j] >= dimension[j]:
                        this_node = False
                    else:
                        flags[j] += 1
                if this_node:
                    new_i = np.sum(flags * 10 ** np.arange(len(flags) - 1,
                                                           -1, -1))
                    self.nids.append(new_i + base)
                    self.energies.append(self.get_energy(flags))
        
        elif self.type == 'linear':
            num = self.ne
            vector = (self.end - self.begin) / (num - 1)
            digits = int(np.ceil(np.log10(num)))
            base = self.index * 10 ** digits
            for i in range(num):
                self.nids.append(i + base)
                self.energies.append(self.begin + i * self.int_step * vector)
        
        elif self.type == 'poles':
            num = int((self.bias_window_end - self.bias_window_begin) //
                                                      self.bias_step + 1)
            real_energies = np.linspace(self.bias_window_begin,
                                        self.bias_window_end, num)
            digits = int(np.ceil(np.log10(num * self.poles_num)))
            base = self.index * 10 ** digits
            for i in range(num):
                for k in range(self.poles_num):
                    self.nids.append(k + i * 10 + base)
                    self.energies.append(real_energies[i] +
                                        (2 * k + 1) * np.pi * self.kt * 1.j)
        else:
            raise RuntimeWarning('Wrong path type %s' % self.type)
            
    def get_new_nids(self, zone, depth=0):
        assert depth < 7
        nids = []
        if self.type == 'Gaussian':
            base = zone * 10 ** (self.maxdepth - depth)
            for i in range(2, 7):
                nids.append(base + i)
        elif self.type == 'linear':
            assert depth == 0
            nids_num = self.ne
            digits = int(np.ceil(np.log10(nids_num)))
            base = self.index * 10 ** digits
            for i in range(1, nids_num + 1):
                nids.append(i + base)
        elif self.type == 'poles':
            assert depth == 0
            base = zone * 100
            nids0 = base + 10 + np.arange(self.poles_num) + 1
            nids1 = base + 20 + np.arange(self.poles_num) + 1
            nids = np.append(nids0, nids1).tolist()
        else:
            raise RuntimeError('Wrong Path Type % s' % self.type)
        return nids
    
    def add_node(self, nid, energy, function, se):
        self.nids.append(nid)
        self.functions.append(function)
        self.energies.append(energy)
        self.ses.append(se)
  
    def get_ends_nids(self):
        nids = []
        num = self.num
        if self.index in [1, 2, 3, 6]:
            nids.append(self.index * num + 1)
        if self.index in [3, 6]:
            nids.append(self.index * num + 7)
        return nids
       
    def get_energy(self, flags):
        if self.type == 'Gaussian':
            pps = np.append([1], self.zone_length[np.array(flags[:-1])])
            lls = []
            for i in range(self.maxdepth):
                lls.append(np.product(pps[:i+1]))
            ss = np.append(self.zone_sample[np.array(flags[:-1])- 1],
                           self.sample[flags[-1] - 1])
            energy = np.sum(ss * lls) * (self.end - self.begin) + self.begin
            
        elif self.type == 'linear':
            ind = self.tonid(flags[1:])
            energy = np.linspace(self.begin, self.end, self.ne)[ind - 1]
       
        elif self.type == 'poles':
            if self.full_nodes:
                num = int((self.bias_window_end - self.bias_window_begin) //
                                                      self.bias_step + 1)
                real_energies = np.linspace(self.bias_window_begin,
                                        self.bias_window_end, num)
                tens = np.arange(len(flags[:-1]) - 1, -1) ** 10
                line_index = np.sum(flags[:-1] * tens)
                energy = real_energies[line_index] + (2 * flags[-1] -
                                                 1) * np.pi * self.kt * 1.j
            else:
                lines = [self.begin, self.end]
                energy = lines[flags[0] - 1] + (2 * flags[-1] -
                                                 1) * np.pi * self.kt * 1.j
        return energy
    
    def get_weight(self, flags):
        if self.type == 'Gaussian':
            wei0 = (self.end - self.begin) / 2.
            for i in range(self.maxdepth - 1):
                wei0 *= self.zone_weight[flags[i + 1]]
            wei1 = wei0
            wei0 *= self.weights0[flags[-1]]
            wei1 *= self.weights1[flags[-1]]
            return wei0, wei1
        elif self.type == 'linear':
            wei = self.weights2 + (self.ne - 6) * 1 + self.weights3
            ind = self.tonid(flags[1:]) - 1
            return wei[ind]
        elif self.type == 'poles':
            return 1

    def tonid(self, flags):
        dim = len(flags)
        nid = 0
        for i in range(dim):
            nid += flags[i] * 10 ** (dim - i - 1)
        return nid

class Contour:
    # see the file description of contour
    ne_err = 1e-4
    calcutype = ['eqInt', 'eqInt', 'eqInt', 'resInt', 'neInt', 'locInt']
    def __init__(self, kt, fermi, bias, maxdepth=7, comm=None, neint='linear',
                  plot_eta=1e-4, neintstep=0.02, eqinttol=1e-4, eta=0.01,
                  min_energy=-700, plot_energy_range=[-5.,5.],
                  plot_energy_point_num=201):
        self.kt = kt
        self.nkt = 8 * self.kt 
        self.dkt = 8 * np.pi * self.kt       
        self.fermi = fermi[0]
        self.bias = bias
        self.neint = neint
        self.leadfermi = []
        for i in range(len(bias)):
            self.leadfermi.append(fermi[i] + bias[i])
        self.minfermi = min(self.leadfermi)
        self.maxfermi = max(self.leadfermi)
        self.plot_eta = plot_eta
        self.eta = eta
        self.dtype = complex
        self.comm = comm
        if self.comm is None:
            self.comm = world
        self.maxdepth = maxdepth
        self.num = 10 ** (self.maxdepth - 1)
        self.converged_zones = []
        self.total_sum = 0
        self.plot_path = None
        self.neintstep = neintstep
        self.eqinttol = eqinttol
        self.min_energy = min_energy
        self.eq_err = eqinttol
        self.plot_energy_range = plot_energy_range
        self.plot_energy_point_num = plot_energy_point_num
        
    def get_dense_contour(self):
        self.paths = []
        depth = self.maxdepth
        self.paths.append(Path(self.min_energy, self.min_energy + 20. * 1.j,
                               1, 1, type='Gaussian', kt=self.kt))
        self.paths.append(Path(self.min_energy + 20. * 1.j, -4. + np.pi * 1.j,
                               2, 3, type='Gaussian', kt=self.kt))
        self.paths.append(Path(-4. + np.pi * 1.j, 4. + np.pi * 1.j, 3, depth,
                              type='linear', kt=self.kt))
        self.paths.append(Path(-3., 3., 4, depth, type='poles', kt=self.kt))
        self.paths.append(Path(-5. + self.eta * 1.j, 5. + self.eta * 1.j, 5,
                              depth, type='linear', kt=self.kt)) 
        for i in range(5):
            self.paths[i].get_full_nids()

    def get_optimized_contour(self, tp):
        tp.log('Contour.get_optimized_contour()')
        self.paths = []
        depth = self.maxdepth
        self.paths.append(Path(self.min_energy + self.minfermi,
                               self.min_energy + self.minfermi +
                               (10. + self.dkt) * 1.j,
                                1, depth,
                               type='Gaussian', kt=self.kt))
        
        self.paths.append(Path(self.min_energy + self.minfermi +
                                (10 + self.dkt) * 1.j,
                                self.minfermi - self.nkt + self.dkt * 1.j,
                                2, depth,
                                type='Gaussian', kt=self.kt))
        
        self.paths.append(Path(self.minfermi - self.nkt + self.dkt * 1.j,
                              self.minfermi + self.nkt + self.dkt * 1.j,
                              3, depth,
                              type='Gaussian', kt=self.kt))
        
        self.paths.append(Path(self.minfermi,
                               self.maxfermi,
                               4, depth,
                               type='poles', kt=self.kt))
        
        self.paths.append(Path(self.minfermi - self.nkt + self.eta * 1.j,
                              self.maxfermi + self.nkt + self.eta * 1.j,
                               5, depth,
                               type=self.neint, kt=self.kt,
                               int_step=self.neintstep))
        
        self.paths.append(Path(self.minfermi - self.nkt + self.dkt * 1.j,
                              self.maxfermi + self.nkt + self.dkt * 1.j,
                               6, depth,
                               type='Gaussian', kt=self.kt))
        
        zones = np.arange(1, 7)
        depth = 0
        converge = False
        while not converge and depth < self.maxdepth:
            nids, path_indices = self.collect(zones, depth)
            loc_nids, loc_path_indices = self.distribute(nids, path_indices)
            self.calculate(tp, loc_nids, loc_path_indices)
            if depth == 0:
                self.joint(zones)
            converge, zones = self.check_convergence(zones, depth, tp.nbmol)
            depth += 1
            self.transfer(zones, depth)

    def get_plot_path(self, ex=False):
        if ex:
            limit = 6.5
        else:    
            limits = self.plot_energy_range
        
        if self.plot_path is None:
            self.plot_path = Path(limits[0] + self.fermi + self.plot_eta * 1.j,
                              limits[1] + self.fermi + self.plot_eta * 1.j,
                               7, 1,
                               type='linear', kt=self.kt)
            path = self.plot_path  
            if ex:
                path.ne = 261
            else:    
                path.ne = self.plot_energy_point_num
            path.int_step = (limits[1] - limits[0]) / (
                                               self.plot_energy_point_num - 1)
           
            digits = int(np.ceil(np.log10(path.ne)))
            base = path.index * 10 ** digits
            energies = np.linspace(path.begin, path.end, path.ne)
            weights = path.weights2 + [1] * (path.ne - 6) + path.weights3
            weights = np.array(weights) * path.int_step
            nids = np.arange(path.ne) + base + 1

            loc_nids = np.array_split(nids, self.comm.size)[self.comm.rank]
            loc_energies = np.array_split(energies,
                                          self.comm.size)[self.comm.rank]
            loc_weights = np.array_split(weights,
                                         self.comm.size)[self.comm.rank]
            path.my_nids = loc_nids
            path.my_energies = loc_energies
            path.my_weights = loc_weights
            path.weights = weights
            path.energies = energies
            path.nids = nids
        return self.plot_path
           
    def collect(self, zones, depth):
        nids = []
        path_indices = []
        for zone in zones:
            path_index = zone // (10 ** depth) - 1
            new_nids = self.paths[path_index].get_new_nids(zone, depth)
            nids += new_nids
            path_indices += [path_index] * len(new_nids)
            
            if depth == 0:
                new_nids = self.paths[path_index].get_ends_nids()
                nids += new_nids
                path_indices += [path_index] * len(new_nids)
        return nids, path_indices
    
    def distribute(self, nids, path_indices):
        loc_nids = np.array_split(nids, self.comm.size)
        loc_path_indices = np.array_split(path_indices, self.comm.size)
        return loc_nids[self.comm.rank], loc_path_indices[self.comm.rank]
    
    def calculate(self, tp, loc_nids, path_indices):
        tp.log('Contour.calculate()')
        for nid, path_index in zip(loc_nids, path_indices):
            exp10 = int(np.floor(np.log10(nid)))
            flags = self.paths[path_index].get_flags(nid, True)
            energy = self.paths[path_index].get_energy(flags[1:])
            if  path_index in [0, 1, 2, 5]:
                calcutype = self.calcutype[path_index]
                green_function, se = tp.calgfunc(energy, calcutype,
                                                      'new')
                self.paths[path_index].functions.append(np.diag(
                                                          green_function[0]))
                self.paths[path_index].energies.append(energy)
                self.paths[path_index].ses.append(se)
                self.paths[path_index].nids.append(nid)
                
    def joint(self, zones):
        my_zones = np.array_split(zones, self.comm.size)[self.comm.rank]
        my_info_dict = {}
        num = 0
        for zone in my_zones:
            if zone in [3, 4, 5, 6]:
                pass
            else:
                path_index = zone - 1
                order = 10 ** self.maxdepth
                base = zone * order
                nid = base + 7
                link_nid = nid + order - 6
                link_path_index = link_nid // (10 ** self.maxdepth) - 1
                flag = str(self.comm.rank) + '_' + str(num)
                my_info_dict[flag] = np.array([nid, link_nid,
                                           path_index, link_path_index], int)
                num += 1
                  
        info_dict = gather_ndarray_dict(my_info_dict, self.comm,
                                                             broadcast=True)
        for name in info_dict:
            nid, link_nid, path_index, link_path_index = info_dict[name]
            rank = self.get_rank(link_path_index, link_nid)
            
            if self.comm.rank == rank:
                link_path = self.paths[link_path_index]                
                index = link_path.nids.index(link_nid)
                function = link_path.functions[index]
                se = link_path.ses[index]
                energy = link_path.energies[index]
                self.paths[path_index].add_node(nid, energy, function, se)

    def transfer(self, zones, depth):
        my_zones = np.array_split(zones, self.comm.size)[self.comm.rank]
        my_info_dict = {}
        num = 0
        maximum = 100000
        name_flags = []
        for zone in my_zones:
            path_index = zone // (10 ** depth) - 1
            order = 10 ** (self.maxdepth - depth)
            base = zone * order
            node_index = zone % 10
            nid = zone * order + 1

            link_nid = (zone - node_index) * order + 2 * node_index - 1
            flag = str((self.comm.rank  + 1) * maximum + num)                
            my_info_dict[flag] = np.array([nid, link_nid, path_index], int)
            num += 1
                
            nid = zone * order + 7
            link_nid = (zone - node_index) * order + 2 * node_index + 1
            flag = str((self.comm.rank + 1) * maximum + num)
            my_info_dict[flag] = np.array([nid, link_nid, path_index], int)
            num += 1                
            
        info_dict = gather_ndarray_dict(my_info_dict, self.comm,
                                        broadcast=True)
        for name in info_dict:
            name_flags.append(eval(name))
        name_flags = np.sort(name_flags)
        #sort is necessrary because sometime the dict sequency is different for
        #different processors!!

        for name_flag in name_flags:
            nid, link_nid, path_index = info_dict[str(name_flag)]
            rank = self.get_rank(path_index, link_nid)
            if self.comm.rank == rank:
                path = self.paths[path_index]                
                index = path.nids.index(link_nid)
                function = path.functions[index]
                energy = path.energies[index]
                se = path.ses[index]
                self.paths[path_index].add_node(nid, energy, function, se)
            
    def get_rank(self, path_index, nid):
        info_array = np.zeros([self.comm.size], int)
        if nid in self.paths[path_index].nids:
            info_array[self.comm.rank] = 1
        self.comm.sum(info_array)
        assert np.sum(info_array) == 1
        return np.argmax(info_array)
       
    def check_convergence(self, zones, depth, nbmol):
        new_zones = []
        converged = True
        errs = [0]
        err = 0
        for zone in zones:
            indices = []
            nids = []
            if zone in [4, 5]:
                pass
            else:
                order = 10 ** (self.maxdepth - depth)
                base = zone * order
                original_nids = np.arange(1, 8) + base
                path_index = zone // 10 ** depth - 1
                gr_sum0 = np.zeros([1, nbmol, nbmol], dtype=self.dtype)
                gr_sum1 = np.zeros([1, nbmol, nbmol], dtype=self.dtype)
                path = self.paths[path_index]            
                for nid in original_nids:
                    if nid in path.nids:
                        flags = path.get_flags(nid, True)
                        weight0, weight1 = path.get_weight(flags)
                        index = path.nids.index(nid)
                        indices.append(index)
                        nids.append(nid)
                        gr_sum0 += path.functions[index] * weight0
                        gr_sum1 += path.functions[index] * weight1
                self.comm.sum(gr_sum0)
                self.comm.sum(gr_sum1)
                err = np.max(abs(gr_sum0 - gr_sum1))
                if err > self.eq_err:
                    converged = False
                    new_zones += range(zone * 10 + 1, zone * 10 + 4)
                    errs.append(err)
            if err < self.eq_err:
                self.converged_zones.append(zone)
                for ind, nid in zip(indices, nids):
                    digit = nid % 10
                    if digit > 1 and digit < 7:            
                        path.functions[ind] = None
        return converged, new_zones        

    def sort_contour(self, nbmol, cal_den=False):
        weights = []
        energies = []
        nids = []
        ses = []
        for zone in self.converged_zones:
            if zone in [4, 5]:
                if cal_den:
                    den_ne = np.zeros([nbmol, nbmol], dtype=self.dtype)
                else:    
                    pass
            else:
                depth = int(np.floor(np.log10(zone)))
                order = 10 ** (self.maxdepth - depth)
                base = zone * order
                original_nids = np.arange(1, 8, 2) + base
                path_index = zone // 10 ** depth - 1
                if cal_den:
                    den_eq = np.zeros([nbmol, nbmol], dtype=self.dtype)
                path = self.paths[path_index]            
                for nid in original_nids:
                    if nid in path.nids:
                        flags = path.get_flags(nid, True)
                        weight0, weight1 = path.get_weight(flags)
                        energy = path.get_energy(flags[1:])   
                        index = path.nids.index(nid)
                        se = path.ses[index]
                        weights.append(weight0)
                        energies.append(energy)
                        nids.append(nid)
                        ses.append(se)
                        if cal_den:
                            den_eq += path.functions[index] * weight0
                if cal_den:
                    self.comm.sum(den_eq)
                    self.comm.sum(den_ne)
        
        seq = np.argsort(nids)
        del_seq = []
        tol = 1e-9
        for i, j in zip(seq[:-1], seq[1:]):
            if np.abs(energies[i] - energies[j]) < tol:
                weights[j] += weights[i]     
                del_seq.append(i)
    
        del_seq.sort()
        del_seq.reverse()
        for i in del_seq:
            del nids[i]
            del energies[i]
            del weights[i]
            del ses[i]
        del self.converged_zones[:]
        
        if cal_den:
            return nids, energies, weights, ses, (den_eq + den_ne)
        else:
            return nids, energies, weights, ses

    def release(self):
        for path in self.paths:
            for function in path.functions:
                function = None
            for se in path.ses:
                se = None
       
    def distribute_nodes(self, path_index):
        path = self.paths[path_index]
        if path.type == 'linear':
            digits = int(np.ceil(np.log10(path.ne)))
            base = path.index * 10 ** digits
            energies = np.linspace(path.begin, path.end, path.ne)
            weights = path.weights2 + [1] * (path.ne - 6) + path.weights3
            weights = np.array(weights) * path.int_step
            nids = np.arange(path.ne) + base + 1
        
        elif path.type == 'poles':
            base = path.index * 100
            nids0 = base + 10 + np.arange(path.poles_num) + 1
            nids1 = base + 20 + np.arange(path.poles_num) + 1
            nids = np.append(nids0, nids1)
            energies0 = path.begin + (np.arange(path.poles_num) * 2
                                                        - 1) * np.pi * 1.j
            energies1 = path.end + (np.arange(path.poles_num) * 2
                                                        - 1) * np.pi * 1.j
            weights0 = [-1] * path.poles_num
            weights1 = [1] * path.poles_num
            weights = np.append(weights0, weights1)
        
        loc_nids = np.array_split(nids, self.comm.size)[self.comm.rank]
        loc_energies = np.array_split(energies,
                                         self.comm.size)[self.comm.rank]
        loc_weights = np.array_split(weights, self.comm.size)[self.comm.rank]
        return loc_nids, loc_energies, loc_weights


