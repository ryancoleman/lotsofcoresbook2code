import sys
import numpy as np
import matplotlib.pyplot as plt
for filename in sys.argv[1:]:
    a = np.loadtxt(filename, delimiter=',').T
    x = a[0]
    for y in a[1:]:
        plt.plot(x, y, '-')
plt.show()
