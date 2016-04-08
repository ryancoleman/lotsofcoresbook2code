from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys

plotter = Transport_Plotter()
plotter.plot_setup()

if len(sys.argv)<=2:
    if len(sys.argv[1]) <= 2:
        nt = plotter.get_info('ntx', int(sys.argv[1]), 0)
    else:
        tmp = sys.argv[1].split('-')
        sam = int(tmp[0])
        ref = int(tmp[1])
        nt = plotter.get_info('ntx', sam, 0) - plotter.get_info('ntx', ref, 0)
else:
    nt = plotter.get_info('ntx', int(sys.argv[1]), int(sys.argv[2]))

matshow(nt[0])
colorbar()
xlabel('Transport Direction')
ylabel('Pseudo Density')
show()
