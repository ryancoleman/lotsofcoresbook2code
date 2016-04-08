import numpy as np

from ase.calculators.calculator import Calculator, all_properties


class SinglePointCalculator(Calculator):
    """Special calculator for a single configuration.

    Used to remember the energy, force and stress for a given
    configuration.  If the positions, atomic numbers, unit cell, or
    boundary conditions are changed, then asking for
    energy/forces/stress will raise an exception."""
    
    name = 'unknown'
    
    def __init__(self, *args, **results):
        """Save energy, forces, stress, ... for the current configuration."""
        if args and isinstance(args[0], float):
            # Old interface:
            assert not results
            for key, value in zip(['energy', 'forces', 'stress', 'magmoms'],
                                  args):
                if value is not None:
                    results[key] = value
            atoms = args[-1]
        else:
            if args:
                atoms = args[0]
            else:
                atoms = results.pop('atoms')
            
        Calculator.__init__(self)
        self.results = {}
        for property, value in results.items():
            assert property in all_properties
            if value is None:
                continue
            if property in ['energy', 'magmom']:
                self.results[property] = value
            else:
                self.results[property] = np.array(value, float)
        self.atoms = atoms.copy()

    def get_property(self, name, atoms):
        if name not in self.results or self.check_state(atoms):
            raise NotImplementedError
        return self.results[name]

    
class SinglePointKPoint:
    def __init__(self, weight, s, k, eps_n=[], f_n=[]):
        self.weight = weight
        self.s = s  # spin index
        self.k = k  # k-point index
        self.eps_n = eps_n
        self.f_n = f_n


class SinglePointDFTCalculator(SinglePointCalculator):
    def __init__(self, *args, **results):
        if args and isinstance(args[0], float):
            # Old interface:
            assert not results
            for key, value in zip(['energy', 'forces', 'stress', 'magmoms'],
                                  args):
                if value is not None:
                    results[key] = value
            atoms = args[4]
            if len(args) > 5:
                eFermi = args[5]
                if len(args) > 6:
                    energies = args[6]
        else:
            if args:
                atoms = args[0]
            else:
                atoms = results.pop('atoms')
            eFermi = results.pop('eFermi', None)
            energies = results.pop('energies', None)
            eref = results.pop('Eref', None)

        SinglePointCalculator.__init__(self, atoms, **results)

        if eFermi is not None:
            self.eFermi = eFermi
        if energies is not None:
            self.energies = energies
        if eref is not None:
            self.eref = eref
        self.kpts = None

    def get_fermi_level(self):
        """Return the Fermi-level(s)."""
        return self.eFermi

    def get_bz_k_points(self):
        """Return the k-points."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return np.zeros((1, 3))
        return None

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return len(self.kpts)
        return None

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        nos = self.get_number_of_spins()
        if nos is not None:
            return nos == 2
        return None
    
    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone."""
        return self.get_bz_k_points()

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.f_n
        return None

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.eps_n
        return None
