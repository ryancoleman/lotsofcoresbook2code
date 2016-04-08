from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys
if '*' in sys.argv[1]:
    fd=0
    nbias = int(sys.argv[1].split('*')[0])
else:
    fd=1
    nbias = int(sys.argv[1])

plotter = Transport_Plotter(fd)
dense_level=1
plotter.plot_setup()
if len(sys.argv) > 1:
    bias, current = plotter.iv(nbias)
else:
    bias, current = plotter.iv()

bias=np.abs(bias)
plot(bias, current, 'r-o')

if dense_level>1:
    from scipy import interpolate
    tck = interpolate.splrep(bias, current, s=0)
    numb = len(bias)
    newbias = np.linspace(bias[0], bias[-1], numb * (dense_level))
    newcurrent = interpolate.splev(newbias, tck, der=0)
    bias=newbias
    current = newcurrent
plot(np.abs(bias), current, 'b-o')
xlabel('Bias(V)')
ylabel('Current($\mu$A)')
show()
