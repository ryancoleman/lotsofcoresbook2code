class Density:
    def __init__(self):
        self.gd = None
        self.mixer = None
        self.nt_sG = None
        self.nct = None

    def read(self, reader):
        self.gd = reader.get_grid_descriptor()
        self.mixer = reader.read('mixer')
        self.nt_sG = reader.get_array('...')

    def update(self, atoms, **kwargs):
        # Initialize stuff:
        if self.gd is None:
            self.gd = ...
        if self.mixer is None:
            self.mixer = Mixer(self.gd)
        if self.nct is None:
            self.nct = LFC(atoms)

        # Change stuff:
        if 'mixer' in kwargs:
            self.mixer = kwargs['mixer']

        # Update stuff:
        self.nct.set_positions(atoms)

    def allocate(self, lfc=True):
        if lfc:
            self.nct.allocate()

    def memory_estimate(self):
        ...

    def write(self, ...):
        ...
