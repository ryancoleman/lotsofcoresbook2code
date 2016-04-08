from gpaw.test import equal
from gpaw.gauss import test_derivatives

for i in range(3):
    i1, i2 = test_derivatives((1.0, -3.4, 1.2),
                             (0, 0, 0), (1, 0, 0), 1.4, 3.0, i)
    equal(i1, i2, 4e-10)
    i1, i2 = test_derivatives((1.0, -3.4, 1.2),
                             (0, 1, 0), (0, 0, 1), 1.4, 3.0, i)
    equal(i1, i2, 2e-10)
    i1, i2 = test_derivatives((1.0, -3.4, 1.2),
                             (0, 1, 0), (1, 0, 1), 1.4, 3.0, i)
    equal(i1, i2, 4e-11)
    i1, i2 = test_derivatives((1.0, -3.4, 1.2),
                             (0, 2, 0), (1, 0, 1), 1.4, 3.0, i)
    equal(i1, i2, 6e-10)


