class PAW(ASECalculator):
    def __init__(self, restart=None, **kwargs):
        ASECalculator.__init__(self)
        self.density = Density()
        self.hamiltonian = Hamiltonian()
        self.wfs = WaveFunctions()

        if restart:
            self.read(Reader(restart))
            self.update(self.atoms, **kwargs)

    def read(self, reader):
        self.atoms = reader.read_atoms()
        self.density.read(reader)
        self.hamiltonian.read(reader)
        self.wfs = self.wfs.read(reader)
        
    def update(self, atoms, **kwargs):
        """Lazy update."""
        self.density.update(self.atoms, kwargs)
        self.hamiltonian.update(self.atoms, kwargs)

        # If we change mode, we could get a completely new type of
        # wave function object:
        self.wfs = self.wfs.update(self.atoms, kwargs)
    
    def set(self, **kwargs):
        self.update(self,atoms, **kwargs)

    def allocate(self, wfs=False, lfc=True):
        self.density.allocate(lfc)
        self.hamiltonian.allocate(lfc)
        self.wfs.allocate(wfs, lfc)

    def calculate(self, atoms):
        self.update(atoms)
        self.allocate(wfs=True)
        ...

    def get_potential_energy(self, atoms):
        ASECalculator.get_potential_energy(self, atoms)
        return self.hamiltonian.energy

