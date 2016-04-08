"""Python wrapper for GPAW's native XC functionals."""

import _gpaw
from gpaw import debug

codes = {
    'LDA': -1,
    'PBE': 0,
    'revPBE': 1,
    'RPBE': 2,
    'PW91': 14,
    'TPSS': 20,
    'M06L': 21,
    'revTPSS': 22
    }
# NOTE: when adding MGGA functionals to the above
# list, self.type must be set to MGGA in XCKernel:__init__
        
class XCNull:
    type = 'LDA'
    name = 'null'
    def calculate(self, e_g, n_sg, dedn_sg):
        e_g[:] = 0.0


class XCKernel:
    def __init__(self, name):
        self.name = name
        if name == 'LDA':
            self.type = 'LDA'
        elif name == 'TPSS' or name == 'M06L' or name == 'revTPSS':
            self.type = 'MGGA'
        else:
            self.type = 'GGA'
        self.xc = _gpaw.XCFunctional(codes[name])
        
    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg=None, dedsigma_xg=None,
                  tau_sg=None, dedtau_sg=None):
        if debug:
            self.check_arguments(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg,
                                 tau_sg, dedtau_sg)
        self.xc.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg, tau_sg, dedtau_sg)

    def check_arguments(self, e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg,
                        tau_sg, dedtau_sg):
        S = n_sg.shape[0]
        G = n_sg.shape[1:]
        assert e_g.shape == G
        assert e_g.flags.contiguous and e_g.dtype == float
        assert dedn_sg.shape == (S,) + G
        assert dedn_sg.flags.contiguous
        assert dedn_sg.dtype == float
        if self.type != 'LDA':
            assert sigma_xg.shape == (2 * S - 1,) + G
            assert dedsigma_xg.shape == (2 * S - 1,) + G
            assert sigma_xg.flags.contiguous and sigma_xg.dtype == float
            assert (dedsigma_xg.flags.contiguous and
                    dedsigma_xg.dtype == float)
            if self.type == 'MGGA':
                assert tau_sg.shape == (S,) + G
                assert dedtau_sg.shape == (S,) + G
                assert tau_sg.flags.contiguous and tau_sg.dtype == float
                assert (dedtau_sg.flags.contiguous and
                        dedtau_sg.dtype == float)
