from gpaw.transport.analysor import Transport_Plotter
import numpy as np
from pylab import *
import sys

plotter = Transport_Plotter()
plotter.plot_setup()
if len(sys.argv)<= 2:
    
    if len(sys.argv[1]) <= 2:
        vt = plotter.get_info('vtx', int(sys.argv[1]), 0)
    else:
        tmp = sys.argv[1].split('-')
        sam = int(tmp[0])
        ref = int(tmp[1])
        vt = plotter.get_info('vtx', sam, 0) - plotter.get_info('vtx', ref, 0)
else:
    vt = plotter.get_info('vtx', int(sys.argv[1]), int(sys.argv[2]))


matshow(vt[0])
colorbar()
xlabel('Transport Direction')
ylabel('Effective Potential')
show()
