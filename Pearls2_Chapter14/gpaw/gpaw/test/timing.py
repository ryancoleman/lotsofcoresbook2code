# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.
from __future__ import print_function
from gpaw.grid_descriptor import GridDescriptor
from gpaw.transformers import Transformer
import time


n = 6
gda = GridDescriptor((n,n,n))
gdb = gda.refine()
gdc = gdb.refine()
a = gda.zeros()
b = gdb.zeros()
c = gdc.zeros()

inter = Transformer(gdb, gdc, 2).apply
restr = Transformer(gdb, gda, 2).apply

t = time.clock()
for i in range(8*300):
    inter(b, c)
print(time.clock() - t)

t = time.clock()
for i in range(8*3000):
    restr(b, a)
print(time.clock() - t)
