from __future__ import print_function
from numpy import around
from gpaw import GPAW
from gpaw.wannier import LocFun

calc = GPAW('CO.gpw', txt=None)
locfun = LocFun()
locfun.localize(calc, M=8, ortho=True, verbose=True)
print(around(locfun.U_nn, 1))

# non ortho
#        O_s  O_py O_pz O_px C_s  C_py C_pz C_px
# s1  [[-0.9  0.   0.4  0.  -0.5  0.  -0.5  0. ]
# s2   [ 0.5  0.   0.8  0.  -0.6  0.  -0.2  0. ]
# pi   [ 0.  -0.4  0.   0.7  0.  -0.3  0.   0.5]
# pi   [ 0.  -0.7  0.  -0.4  0.  -0.5  0.  -0.3]
# s3   [ 0.   0.   0.4  0.   0.6  0.  -0.8  0. ]
# pi*  [ 0.  -0.2  0.  -0.5  0.   0.2  0.   0.8]
# pi*  [ 0.  -0.5  0.   0.2  0.   0.8  0.  -0.2]
# s4   [-0.2  0.  -0.1  0.  -0.2 -0.   0.2  0. ]]

# ortho
#        O_s  O_py O_pz O_px C_s  C_py C_pz C_px
# s1  [[-0.7  0.   0.1  0.  -0.5  0.  -0.4  0. ]
# s2   [ 0.5  0.   0.7  0.  -0.5  0.  -0.1  0. ]
# pi   [ 0.  -0.4  0.   0.7  0.  -0.3  0.   0.5]
# pi   [ 0.  -0.7  0.  -0.4  0.  -0.5  0.  -0.3]
# s3   [-0.1  0.   0.5  0.   0.7  0.  -0.5  0. ]
# pi*  [ 0.  -0.2  0.  -0.5  0.   0.2  0.   0.8]
# pi*  [ 0.  -0.5  0.   0.2  0.   0.8  0.  -0.2]
# s4   [-0.5  0.   0.5  0.   0.2  0.   0.7  0. ]]
