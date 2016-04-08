"""Small test for local parts using yukawa potential"""

from numpy import exp, sqrt, pi, absolute
from gpaw.utilities import fact
from gpaw.atom.radialgd import AERadialGridDescriptor

#Values from Rico, Lopez, Ramirez, Ema, Theor Chem Acc (2013) 132:1304
#   Table 1 and 2.
#
#   Normalized STO: N = n + L : lambda (gamma) = 0.34
#   n    n'   L   Zeta     Zeta'      I(n,n',L,Zeta,Zeta')
#   8    6    4   8.4      3.5        0.793087136533672
#   2    1    5   8.4      3.5        0.271337849250231
#   3    4    6   8.4      3.5        0.170346973782779
#   3    2    8   8.4      3.5        0.141826311951003
#   2    1   10   8.4      3.5        9.005003835363871(-2)

#
#   Zeta = lambda (gamma)
#   n = 3, L=5, n'=2, L'=5, zeta'=3.5
#   0.10                              9.503545229396430(-6)
#   0.25                              5.869262722600324(-4)
#   0.50                              1.194846721531269(-2)
#   0.75                              5.829196062489732(-2)
#   1.00                              0.153003818405446


basic_sto_list = [  # See table 1 of Rico, Lopez, Ramirez, Ema,
                    # Theor Chem Acc (2013) 132:1304
#        {
#            'n': [8, 6],            # n, n'
#            'Zeta': [8.4, 3.5],     # Z, Z'
#            'L': 4,                 # L
#            'I': 0.793087136533672  # [X_LM^n|X_LM^n']
#            },
#        {
#            'n': [2, 1],
#            'Zeta': [8.4, 3.5],
#            'L': 5,
#            'I': 0.271337849250231
#            },
#        {          # values are completely of
#            'n': [3, 4],
#            'Zeta': [8.4, 3.5],
#            'L': 6,
#            'I': 0.170346973782779
#            },
#        {
#            'n': [3, 2],            # n, n'
#            'Zeta': [8.4, 3.5],     # Z, Z'
#            'L': 8,                 # L
#            'I': 0.141826311951003  # [X_LM^n|X_LM^n']
#            },
        {
            'n': [2, 1],
            'Zeta': [8.4, 3.5],
            'L': 10,
            'I': 9.005003835363871e-2
            }
        ]
        
gamma_sto_list = [  # See table 2 of Rico, Lopez, Ramirez, Ema,
                    # Theor Chem Acc (2013) 132:1304
#        {
#            'zeta': 0.10,
#            'I': 9.503545229396430e-6
#        },
#        {
#            'zeta': 0.25,
#            'I': 5.869262722600324e-4
#        },
#        {
#            'zeta': 0.50,
#            'I': 1.194846721531269e-2
#        },
#        {
#            'zeta': 0.75,
#            'I': 5.829196062489732e-2
#        },
        {
            'zeta': 1.00,
            'I': 0.153003818405446
        },
        ]


def radial_sto(n, zeta, l, r):
    """Build radial part of slater type orbital"""
#       Stolen from gpaw all_electron.py (intialize_wave_functions)
#
#       STOs are defined as
#
#           N Y_l^m r^{n-1} exp(-zeta * r)
#
#           N = (2 zeta)^n sqrt(2 zeta / (2n)!)
    assert n > 0
    radial = r**(n + l + 0.5)
    radial *= exp(-zeta * r)
    N = (2 * zeta)**(n + l) * sqrt(2 * zeta / fact(2 * (n + l)))
    radial *= N  # perhaps also do Y_l^m normalization
    return radial


def test_different_sto(rgd):
    """Check integration of two different STO."""
    gamma = 0.34
    for test_sto in basic_sto_list:
        wave_functions = []
        l_sto = test_sto['L']
        for n, zeta in zip(test_sto['n'], test_sto['Zeta']):
            radial = radial_sto(n, zeta, l_sto, rgd.r_g)
            norm = sqrt((2 * l_sto + 1) / (4 * pi))  # m = 0
            radial *= norm
            wave_functions.append([radial])
        I = rgd.integrate_yukawa(wave_functions[0], wave_functions[1], 
                l_sto, gamma)
        scale = 16 * pi**2 / (2 * l_sto + 1)
        I *= scale
#        print(u"{:7.5e}||{:7.5e}||{:7.5e}".format(test_sto['I'], I,
#                absolute(I - test_sto['I']) * 100.0 / test_sto['I']))
        assert (absolute(I - test_sto['I'])/test_sto['I'] < 0.001)


def test_same_sto(rgd):
    """Check integration of nearly identical STO."""
    n = 3
    l_sto = 5
    np = 2
    zetap = 3.5
    for test_sto in gamma_sto_list:
        zeta = test_sto['zeta']
        radial = radial_sto(n, zeta, l_sto, rgd.r_g)
        norm = sqrt((2 * l_sto + 1) / (4 * pi))  # m = 0
        radial *= norm
        radial2 = radial_sto(np, zetap, l_sto, rgd.r_g)
        radial2 *= norm
        I = rgd.integrate_yukawa(radial, radial2, l_sto, zeta)
        scale = 16 * pi**2 / (2 * l_sto + 1)
        I *= scale
#        print(u"{:7.5e}||{:7.5e}||{:7.5e}".format(test_sto['I'], I,
#                absolute(I - test_sto['I']) * 100.0 / test_sto['I']))
        assert (absolute(I - test_sto['I'])/test_sto['I'] < 0.001)


def main():
    """Work on it."""
    N = 1200
    beta = 0.4 * 600 / N
    rgd = AERadialGridDescriptor(beta / N, 1.0 / N, N)
    test_same_sto(rgd)
    test_different_sto(rgd)

if __name__ == '__main__':
    main()
