import numpy as np

from gpaw.xc.hybrid import HybridXCBase


class ForceCalculator:
    def __init__(self, timer):
        self.timer = timer
        self.reset()
        
    def reset(self):
        self.F_av = None

    def calculate(self, wfs, dens, ham):
        """Return the atomic forces."""

        assert not isinstance(ham.xc, HybridXCBase)

        if self.F_av is not None:
            return self.F_av

        self.timer.start('Force calculation')

        natoms = len(wfs.setups)
        self.F_av = np.zeros((natoms, 3))

        # Force from projector functions (and basis set):
        wfs.calculate_forces(ham, self.F_av)
        
        try:
            # ODD functionals need force corrections for each spin
            correction = ham.xc.setup_force_corrections
        except AttributeError:
            pass
        else:
            correction(self.F_av)
        
        if wfs.bd.comm.rank == 0 and wfs.kd.comm.rank == 0:
            ham.calculate_forces(dens, self.F_av)

        wfs.world.broadcast(self.F_av, 0)
        
        self.F_av = wfs.kd.symmetry.symmetrize_forces(self.F_av)

        self.timer.stop('Force calculation')

        return self.F_av
