import numpy as np

from ase.parallel import parprint

from gpaw.lrtddft import LrTDDFT
from gpaw.utilities.folder import Folder

from gpaw.inducedfield.inducedfield_base import BaseInducedField


class LrTDDFTInducedField(BaseInducedField):
    """Induced field class for linear response TDDFT.
    
    Attributes (see also ``BaseInducedField``):
    -------------------------------------------
    lr: LrTDDFT object
        LrTDDFT object
    """
    
    def __init__(self, filename=None, paw=None, lr=None, ws='all',
                  frequencies=None, folding='Gauss', width=0.08,
                  kickdir=0
                  ):
        """
        Parameters (see also ``BaseInducedField``):
        -------------------------------------------
        paw: GPAW object
            GPAW object for ground state
        lr: LrTDDFT object
            LrTDDFT object
        kickdir: int
            Kick direction: 0, 1, 2 for x, y, z
        """
        # TODO: change kickdir to general kick = [1e-3, 1-e3, 0] etc.

        if type(lr) is not LrTDDFT:
            raise RuntimeError('Provide LrTDDFT object.')
        
        # Check that lr is diagonalized
        if len(lr) == 0:
            raise RuntimeError('Diagonalize LrTDDFT first.')
        
        self.lr = lr
        
        # "Kick" direction
        self.kickdir = kickdir
        
        self.readwritemode_str_to_list = \
        {'': ['Frho', 'atoms'],
         'all': ['Frho', 'Fphi', 'Fef', 'Ffe', 'atoms'],
         'field': ['Frho', 'Fphi', 'Fef', 'Ffe', 'atoms']}
        
        BaseInducedField.__init__(self, filename, paw, ws,
                                  frequencies, folding, width)
    
    def initialize(self, paw, allocate=True):
        BaseInducedField.initialize(self, paw, allocate)

        # Artificial background electric field (linear response)
        # TODO: change kickdir to general kick = [1e-3, 1-e3, 0] etc.
        Fbgef_v = np.zeros((self.nv,), dtype=float)
        Fbgef_v[self.kickdir] = 1.0
        self.Fbgef_v = Fbgef_v

        # Initialize PAW object
        self.density.ghat.set_positions(self.atoms.get_scaled_positions())
        paw.converge_wave_functions()
        
    def set_folding(self, folding, width):
        BaseInducedField.set_folding(self, folding, width)
        
        if self.folding is None:
            self.folder = Folder(None, None)
        else:
            if self.folding == 'Gauss':
                self.folder = Folder(self.width, 'ComplexGauss')
            elif self.folding == 'Lorentz':
                self.folder = Folder(self.width, 'ComplexLorentz')
            else:
                raise RuntimeError('unknown folding "' + self.folding + '"')

    def get_induced_density(self, from_density='comp', gridrefinement=2):
        
        paw = self.paw
        lr = self.lr
        omega_w = self.omega_w

        if gridrefinement == 1:
            gd = self.gd
        elif gridrefinement == 2:
            gd = self.density.finegd
        else:
            raise NotImplementedError
        
        nkss = len(lr.kss)
        nexcs = len(lr)
        omega_I = np.empty((nexcs), dtype=float)
        
        # C_Iij: weights for each excitation I and KS-pair ij
        C_Iij = np.zeros((nexcs, nkss), dtype=float)
        
        # Calculate C_Iij
        FIsqfe_ij = np.zeros((nkss), dtype=float)
        for I, exc in enumerate(lr):
            omega_I[I] = exc.energy
        
            B = 0.0
            for ij, ks in enumerate(lr.kss):
                FIsqfe_ij[ij] = exc.f[ij] * np.sqrt(ks.fij * ks.energy)
                B += ks.mur[self.kickdir] * FIsqfe_ij[ij]
            
            C_Iij[I] = 2 * B * FIsqfe_ij

        # Fold C_Iij to C_wij
        # C_wij: weights for each requested frequency w and KS-pair ij
        self.omega_w, C_wij = self.folder.fold_values(omega_I, C_Iij, omega_w)
        
        assert (self.omega_w == omega_w).all()  # TODO: remove
        
        Fn_wg = gd.zeros((self.nw), dtype=complex)
        
        kpt_u = paw.wfs.kpt_u
        
        for ij, ks_ in enumerate(lr.kss):
            # TODO: better output of progress
            parprint('%d/%d' % (ij, nkss))
            
            # TODO: speedup:
            # for each w, check magnitude of C_wij and
            # take C_wij for those ij that contribute most
            
            # Take copy so that memory for pair density is released
            # after each loop iteration
            ks = ks_.copy()
            
            # Initialize pair density
            kpt = kpt_u[ks.spin]
            ks.set_paw(paw)
            ks.initialize(kpt, ks.i, ks.j)
            
            # Get pair density
            # TODO: gridrefinements
            yes_finegrid = gridrefinement == 2
            if from_density == 'pseudo':
                nij_g = ks.get(finegrid=yes_finegrid)
            elif from_density == 'comp':
                nij_g = ks.with_compensation_charges(finegrid=yes_finegrid)
            elif from_density == 'ae':
                nij_g = ks.with_ae_corrections(finegrid=yes_finegrid)
            # TODO: better to split pair density in pseudo density part
            # and FD_asp etc. coefficients (like in time-propagation).
            # Then do gridrefinement and add AE/comp corrections afterwards
            # -> speedup (no need to do summations on fine grid)
            
            # TODO: speedup: check magnitude of C_wij (see above)
            for w in range(self.nw):
                Fn_wg[w] += nij_g * C_wij[w, ij]
        
        # Return charge density (electrons = negative charge)
        return - Fn_wg, gd
