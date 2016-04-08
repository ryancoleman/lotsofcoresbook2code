from __future__ import print_function
from gpaw.lrtddft import d2Excdnsdnt, d2Excdn2
import numpy as np

dup=np.array([.1,.3])
ddn=np.array([.2,.005])

def equiv(xl,yl):
    for x,y in zip(xl,yl):
#        print x,y
        assert abs(x-y) < 1e-8

out=0

# unpolarised ......................................

expected=[ -0.99599492, -12.51604312]
res=d2Excdn2(ddn)
if out: print(res)
equiv(xl=res,yl=expected)

# polarised ........................................

isp=1
ksp=1
expected=[-1.6505808,  -0.93437538]
res=d2Excdnsdnt(dup,ddn)
if out: print(res)
equiv(res[isp][ksp],expected)

isp=1
ksp=0
expected=[-0.18324442, -0.11948321]
res=d2Excdnsdnt(dup,ddn)
if out: print(res)
equiv(res[isp][ksp],expected)

isp=0
ksp=0
expected=[ -1.14670206, -11.02164441]
res=d2Excdnsdnt(dup,ddn)
if out: print(res)
equiv(res[isp][ksp],expected)

