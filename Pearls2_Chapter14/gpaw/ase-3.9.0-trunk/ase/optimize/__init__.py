"""Structure optimization. """

from ase.optimize.optimize import NDPoly, polyfit
from ase.optimize.mdmin import MDMin
from ase.optimize.lbfgs import HessLBFGS, LineLBFGS
from ase.optimize.fire import FIRE
from ase.optimize.lbfgs import LBFGS, LBFGSLineSearch
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.optimize.bfgs import BFGS
from ase.optimize.oldqn import GoodOldQuasiNewton

try:
    import scipy.optimize
    from ase.optimize.sciopt import Converged
    from ase.optimize.sciopt import SciPyFminCG
    from ase.optimize.sciopt import SciPyFminBFGS
    from ase.optimize.sciopt import SciPyGradientlessOptimizer
    from ase.optimize.sciopt import SciPyFmin
    from ase.optimize.sciopt import SciPyFminPowell
except ImportError:
    pass

QuasiNewton = BFGSLineSearch
