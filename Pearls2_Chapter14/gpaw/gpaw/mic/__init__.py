import pymic
from gpaw.mpi import rank
from gpaw import use_mic
import os

# Distribute multiple devices to ranks. Assumes that ranks are assigned to
# nodes sequencely i.e. with 4 ranks and 2 nodes ranks 0,1 reside in node 0
# and ranks 2,3 in node 1
if use_mic:
    ppn=int(os.getenv("GPAW_PPN", 1))
    rpd=ppn / pymic.number_of_devices()
    if rpd == 0:
        rpd = 1
    dev_id = (rank / rpd) % 2
    stream = pymic.devices[dev_id].get_default_stream()
    print "GPAW: rank {0:02d} will be using offload device {1}".format(rank, dev_id)
else:
    dev_id = -1
    stream = None
