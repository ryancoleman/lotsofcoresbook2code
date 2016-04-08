from __future__ import print_function
from gpaw.utilities import erf, cerf

# Mathematica says:
#   z             erf(z)
# 1           0.842701
# 3           0.999978
# I           0. + 1.65043 I
# 1 + I       1.31615 + 0.190453 I
# 0.3 + 3*I   1467.69 - 166.561 I
# 3 - 0.3*I   1.00001 - 0.0000228553 I
values = [
    [ 1, 0.84270079295+0j],
    [ 3, 0.999977909503+0j],
    [ 1j, 1.6504257588j],
    [ 1+1j, 1.3161512816979477+0.19045346923783471j ],
    [ 0.3 + 3j, 1467.69028322-166.560924526j ],
    [ 3 + 0.3j, 0.99997602085736015+2.1863701577230078e-06j]]

maxerror = 1.e-10

for test in values:
    z, res = test
    try:
        r = z.real
    except:
        r = z
    error = abs(res / cerf(z) - 1.)
    if error < maxerror:
        print('z=', z, ' ok (error=', error, ')')
    else:
        print(z, res, cerf(z), erf(r), error)
        string = ('error for erf(' + str(z) +') = ' + str(error) + 
                  ' > ' + str(maxerror)) 
        assert(error < maxerror), string
