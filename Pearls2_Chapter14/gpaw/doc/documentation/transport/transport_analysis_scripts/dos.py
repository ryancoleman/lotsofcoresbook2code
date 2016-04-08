from gpaw.transport.analysor import Transport_Plotter
import numpy as np
import sys
from pylab import *

if '*' in sys.argv[1]:
    fd=0
    bias_step = int(sys.argv[1].split('*')[0])
else:
    fd=1
    bias_step = int(sys.argv[1])

plotter=Transport_Plotter(fd)
plotter.plot_setup()
dos = plotter.dos(bias_step)
ee=np.linspace(-5,5,201)
plot(ee, dos, 'b-o')
dense_level=1
if dense_level>1:
    from scipy import interpolate
    tck = interpolate.splrep(ee, dos, s=0)
    numb = len(ee)
    newee = np.linspace(ee[0], ee[-1], numb * (dense_level))
    newtc = interpolate.splev(newee, tck, der=0)
    ee = newee
    dos = newdos
    plot(ee, dos, 'r-o')

eye = np.zeros([10, 1]) + 1
bias = plotter.get_info('bias', bias_step)
f1 = bias[0] * eye
f2 = bias[1] * eye        
a1 = np.max(dos)
l1 = np.linspace(0, a1, 10)

plot(f1, l1, 'r--')
plot(f2, l1, 'r--')
show()


