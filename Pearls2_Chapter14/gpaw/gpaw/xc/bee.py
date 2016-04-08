import numpy as np
from types import FloatType
from ase.parallel import rank
from ase.units import Hartree

import _gpaw
from gpaw.xc import XC
from gpaw.xc.kernel import XCKernel
from gpaw.xc.libxc import LibXC
from gpaw.xc.vdw import VDWFunctional
from gpaw import debug


class BEE2(XCKernel):
    """GGA exchange expanded in Legendre polynomials."""
    def __init__(self, parameters=None):
        """BEE2.

        parameters: array
            [transformation,0.0,[orders],[coefs]].

        """

        if parameters is None:
            parameters = (
                [4.0, 0.0] + range(30) +
                [1.516501714304992365356, 0.441353209874497942611,
                 -0.091821352411060291887, -0.023527543314744041314,
                 0.034188284548603550816, 0.002411870075717384172,
                 -0.014163813515916020766, 0.000697589558149178113,
                 0.009859205136982565273, -0.006737855050935187551,
                 -0.001573330824338589097, 0.005036146253345903309,
                 -0.002569472452841069059, -0.000987495397608761146,
                 0.002033722894696920677, -0.000801871884834044583,
                 -0.000668807872347525591, 0.001030936331268264214,
                 -0.000367383865990214423, -0.000421363539352619543,
                 0.000576160799160517858, -0.000083465037349510408,
                 -0.000445844758523195788, 0.000460129009232047457,
                 -0.000005231775398304339, -0.000423957047149510404,
                 0.000375019067938866537, 0.000021149381251344578,
                 -0.000190491156503997170, 0.000073843624209823442])
        else:
            assert len(parameters) > 2
            assert np.mod(len(parameters), 2) == 0
            assert parameters[1] == 0.0

        parameters = np.array(parameters, dtype=float).ravel()
        self.xc = _gpaw.XCFunctional(17, parameters)
        self.type = 'GGA'
        self.name = 'BEE2'


class BEEVDWKernel(XCKernel):
    """Kernel for BEEVDW functionals."""
    def __init__(self, bee, xcoefs, ldac, ggac):
        """BEEVDW kernel.

        parameters:

        bee : str
            choose BEE1 or BEE2 exchange basis expansion.
        xcoefs : array
            coefficients for exchange.
        ldac : float
            coefficient for LDA correlation.
        pbec : float
            coefficient for PBE correlation.

        """

        if bee == 'BEE2':
            self.BEE = BEE2(xcoefs)
            self.GGAc = LibXC('GGA_C_PBE')
            self.xtype = 'GGA'
            self.type = 'GGA'
        elif bee == 'BEE3':
            self.BEE = LibXC('MGGA_X_MBEEFVDW')
            self.GGAc = LibXC('GGA_C_PBE_SOL')
            self.xtype = 'MGGA'
            self.type = 'MGGA'
        else:
            raise ValueError('Unknown BEE exchange: %s', bee)

        self.LDAc = LibXC('LDA_C_PW')
        self.ldac = ldac
        self.ggac = ggac
        if bee in ['BEE1', 'BEE2']:
            self.ggac -= 1.0
        self.name = 'BEEVDW'

    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg=None, dedsigma_xg=None,
                  tau_sg=None, dedtau_sg=None):
        if debug:
            self.check_arguments(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg,
                                 tau_sg, dedtau_sg)

        if self.xtype == 'GGA':
            self.BEE.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg)
        elif self.xtype == 'MGGA':
            self.BEE.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg,
                               tau_sg, dedtau_sg)
        else:
            raise ValueError('Unexpected value of xtype:', self.xtype)

        e0_g = np.empty_like(e_g)
        dedn0_sg = np.empty_like(dedn_sg)
        dedsigma0_xg = np.empty_like(dedsigma_xg)
        for coef, kernel in [(self.ldac, self.LDAc),
                             (self.ggac, self.GGAc)]:
            dedn0_sg[:] = 0.0
            kernel.calculate(e0_g, n_sg, dedn0_sg, sigma_xg, dedsigma0_xg)
            e_g += coef * e0_g
            dedn_sg += coef * dedn0_sg
            if kernel.type == 'GGA':
                dedsigma_xg += coef * dedsigma0_xg


class BEEFEnsemble:
    """BEEF ensemble error estimation."""
    def __init__(self, calc):
        """BEEF ensemble

        parameters:

        calc: object
            Calculator holding a selfconsistent BEEF type electron density.
            May be BEEF-vdW or mBEEF.
        """

        self.calc = calc

        self.e_dft = None
        self.e0 = None

        # determine functional and read parameters
        self.xc = self.calc.get_xc_functional()
        if self.xc == 'BEEF-vdW':
            self.bee_type = 1
        elif self.xc == 'mBEEF':
            self.bee_type = 2
            self.max_order = 8
            self.trans = [6.5124, -1.0]
            self.calc.converge_wave_functions()
            if rank == 0:
                print('wave functions converged')
        elif self.xc == 'mBEEF-vdW':
            self.bee_type = 3
            self.max_order = 5
            self.trans = [6.5124, -1.0]
            self.calc.converge_wave_functions()
            if rank == 0:
                print('wave functions converged')
        else:
            raise NotImplementedError('xc = %s not implemented' % self.xc)

    def create_xc_contributions(self, type):
        """General function for creating exchange or correlation energies"""
        assert type in ['exch', 'corr']
        err = 'bee_type %i not implemented' % self.bee_type

        if type == 'exch':
            if self.bee_type == 1:
                out = self.beefvdw_energy_contribs_x()
            elif self.bee_type in [2, 3]:
                out = self.mbeef_exchange_energy_contribs()
            else:
                raise NotImplementedError(err)
        else:
            if self.bee_type == 1:
                out = self.beefvdw_energy_contribs_c()
            elif self.bee_type == 2:
                out = np.array([])
            elif self.bee_type == 3:
                out = self.mbeefvdw_energy_contribs_c()
            else:
                raise NotImplementedError(err)
        return out

    def get_non_xc_total_energies(self):
        """Compile non-XC total energy contributions"""
        if self.e_dft is None:
            self.e_dft = self.calc.get_potential_energy()
        if self.e0 is None:
            from gpaw.xc.kernel import XCNull
            xc_null = XC(XCNull())
            self.e0 = self.e_dft + self.calc.get_xc_difference(xc_null)
        assert isinstance(self.e_dft, FloatType)
        assert isinstance(self.e0, FloatType)

    def mbeef_exchange_energy_contribs(self):
        """Legendre polynomial exchange contributions to mBEEF Etot"""
        self.get_non_xc_total_energies()
        e_x = np.zeros((self.max_order, self.max_order))
        for p1 in range(self.max_order):  # alpha
            for p2 in range(self.max_order):  # s2
                pars_i = np.array([1, self.trans[0], p2, 1.0])
                pars_j = np.array([1, self.trans[1], p1, 1.0])
                pars = np.hstack((pars_i, pars_j))
                x = XC('2D-MGGA', pars)
                e_x[p1, p2] = (self.e_dft +
                               self.calc.get_xc_difference(x) - self.e0)
                del x
        return e_x

    def beefvdw_energy_contribs_x(self):
        """Legendre polynomial exchange contributions to BEEF-vdW Etot"""
        self.get_non_xc_total_energies()
        e_pbe = (self.e_dft + self.calc.get_xc_difference('GGA_C_PBE') -
                 self.e0)

        exch = np.zeros(30)
        for p in range(30):
            pars = [4, 0, p, 1.0]
            bee = XC('BEE2', pars)
            exch[p] = (self.e_dft + self.calc.get_xc_difference(bee) -
                       self.e0 - e_pbe)
            del bee
        return exch

    def beefvdw_energy_contribs_c(self):
        """LDA and PBE correlation contributions to BEEF-vdW Etot"""
        self.get_non_xc_total_energies()
        e_lda = self.e_dft + self.calc.get_xc_difference('LDA_C_PW') - self.e0
        e_pbe = self.e_dft + self.calc.get_xc_difference('GGA_C_PBE') - self.e0
        corr = np.array([e_lda, e_pbe])
        return corr

    def mbeefvdw_energy_contribs_c(self):
        """LDA, PBEsol, and nl2 correlation contributions to mBEEF-vdW Etot"""
        self.get_non_xc_total_energies()
        e_lda = self.e_dft + self.calc.get_xc_difference('LDA_C_PW') - self.e0
        e_sol = (self.e_dft + self.calc.get_xc_difference('GGA_C_PBE_SOL') -
                 self.e0)
        vdwdf2 = VDWFunctional('vdW-DF2')
        self.calc.get_xc_difference(vdwdf2)
        e_nl2 = vdwdf2.get_Ecnl() * Hartree
        corr = np.array([e_lda, e_sol, e_nl2])
        return corr
