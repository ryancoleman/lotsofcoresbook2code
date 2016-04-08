import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('eels.csv', delimiter=',')
plt.plot(data[:,0], data[:,2])
plt.savefig('aluminum_EELS.png', bbox_inches='tight')
plt.show()
