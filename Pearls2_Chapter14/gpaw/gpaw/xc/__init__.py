from __future__ import print_function
from gpaw.xc.libxc import LibXC
from gpaw.xc.lda import LDA
from gpaw.xc.gga import GGA
from gpaw.xc.mgga import MGGA


def XC(kernel, parameters=None):
    """Create XCFunctional object.

    kernel: XCKernel object or str
        Kernel object or name of functional.
    parameters: ndarray
        Parameters for BEE functional.

    Recognized names are: LDA, PW91, PBE, revPBE, RPBE, BLYP, HCTH407,
    TPSS, M06L, revTPSS, vdW-DF, vdW-DF2, EXX, PBE0, B3LYP, BEE,
    GLLBSC.  One can also use equivalent libxc names, for example
    GGA_X_PBE+GGA_C_PBE is equivalent to PBE, and LDA_X to the LDA exchange.
    In this way one has access to all the functionals defined in libxc.
    See xc_funcs.h for the complete list.  """
    
    if isinstance(kernel, str):
        name = kernel
        if name in ['vdW-DF', 'vdW-DF2', 'optPBE-vdW', 'optB88-vdW',
                    'C09-vdW', 'mBEEF-vdW', 'BEEF-vdW']:
            from gpaw.xc.vdw import VDWFunctional
            return VDWFunctional(name)
        elif name in ['EXX', 'PBE0', 'B3LYP']:
            from gpaw.xc.hybrid import HybridXC
            return HybridXC(name)
        elif name in ['HSE03', 'HSE06']:
            from gpaw.xc.exx import EXX
            return EXX(name)
        elif name == 'BEE1':
            from gpaw.xc.bee import BEE1
            kernel = BEE1(parameters)
        elif name == 'BEE2':
            from gpaw.xc.bee import BEE2
            kernel = BEE2(parameters)
        elif name.startswith('GLLB'):
            from gpaw.xc.gllb.nonlocalfunctionalfactory import \
                NonLocalFunctionalFactory
            xc = NonLocalFunctionalFactory().get_functional_by_name(name)
            xc.print_functional()
            return xc
        elif name == 'LB94':
            from gpaw.xc.lb94 import LB94
            kernel = LB94()
        elif name.startswith('ODD_'):
            from ODD import ODDFunctional
            return ODDFunctional(name[4:])
        elif name.endswith('PZ-SIC'):
            try:
                from ODD import PerdewZungerSIC as SIC
                return SIC(xc=name[:-7])
            except:
                from gpaw.xc.sic import SIC
                return SIC(xc=name[:-7])
        elif name == 'TPSS' or name == 'M06L' or name == 'revTPSS':
            from gpaw.xc.kernel import XCKernel
            kernel = XCKernel(name)
        elif name.startswith('old'):
            from gpaw.xc.kernel import XCKernel
            kernel = XCKernel(name[3:])
        elif name == 'PPLDA':
            from gpaw.xc.lda import PurePythonLDAKernel
            kernel = PurePythonLDAKernel()
        elif name in ['pyPBE', 'pyPBEsol', 'pyRPBE', 'pyzvPBEsol']:
            from gpaw.xc.gga import PurePythonGGAKernel
            kernel = PurePythonGGAKernel(name)
        elif name == '2D-MGGA':
            from gpaw.xc.mgga import PurePython2DMGGAKernel
            kernel = PurePython2DMGGAKernel(name, parameters)
        elif name[0].isdigit():
            from gpaw.xc.parametrizedxc import ParametrizedKernel
            kernel = ParametrizedKernel(name)
        else:
            kernel = LibXC(kernel)
    if kernel.type == 'LDA':
        return LDA(kernel)
    elif kernel.type == 'GGA':
        return GGA(kernel)
    else:
        return MGGA(kernel)

        
def xc(filename, xc, ecut=None):
    """Calculate non self-consitent energy.
    
    filename: str
        Name of restart-file.
    xc: str
        Functional
    ecut: float
        Plane-wave cutoff for exact exchange.
    """
    name, ext = filename.rsplit('.', 1)
    assert ext == 'gpw'
    if xc in ['EXX', 'PBE0', 'B3LYP']:
        from gpaw.xc.exx import EXX
        exx = EXX(filename, xc, ecut=ecut, txt=name + '-exx.txt')
        exx.calculate()
        e = exx.get_total_energy()
    else:
        from gpaw import GPAW
        calc = GPAW(filename, txt=None)
        e = calc.get_potential_energy() + calc.get_xc_difference(xc)
    print(e, 'eV')
