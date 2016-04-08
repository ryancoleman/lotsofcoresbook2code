class ASECalculator:
    def __init__(self):
        self.atoms = None
        
    def get_potential_energy(self, atoms):
        if self.calculation_required(atoms, 'energy'):
            self.calculate(atoms)
            self.atoms = atoms.copy()  # store copy of last configuration
        return None

    ...
