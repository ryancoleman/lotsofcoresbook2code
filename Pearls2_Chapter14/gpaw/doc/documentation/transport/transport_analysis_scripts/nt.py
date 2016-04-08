from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys

plotter = Transport_Plotter()
plotter.plot_setup()
if len(sys.argv) <= 2:
    if len(sys.argv[1]) <= 2:
        nt = plotter.get_info('nt', int(sys.argv[1]), 0)
    else:
        tmp = sys.argv[1].split('-')
        sam = int(tmp[0])
        ref = int(tmp[1])
        nt = plotter.get_info('nt', sam, 0) - plotter.get_info('nt', ref, 0)
else:
    nt = plotter.get_info('nt', int(sys.argv[1]), int(sys.argv[2]))
   

plot(nt[0], 'b-o')
xlabel('Transport Direction')
ylabel('Pseudo Density')
show()
