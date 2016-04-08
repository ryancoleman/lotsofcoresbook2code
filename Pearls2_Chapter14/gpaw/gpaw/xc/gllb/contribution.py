class Contribution:
    def __init__(self, nlfunc, weight):
        self.weight = weight
        self.nlfunc = nlfunc
        nlfunc.add_contribution(self)

    def initialize(self):
        raise NotImplementedError

    def initialize_1d(self):
        raise NotImplementedError

    def get_name(self):
        raise NotImplementedError

    def set_positions(self, spos_ac):
        pass

    def get_desc(self):
        raise NotImplementedError

    def calculate_spinpaired(self, e_g, n_g, v_g):
        raise NotImplementedError

    def calculate_spinpolarized(self, e_g, na_g, va_g, nb_g, vb_g):
        raise NotImplementedError

    def calculate_energy_and_derivatives(self, setup, D_sp, H_sp):
        raise NotImplementedError
    
    def add_smooth_xc_potential_and_energy_1d(self, vt_g):
        raise NotImplementedError

    def add_extra_setup_date(self, dict):
        raise NotImplementedError
        
