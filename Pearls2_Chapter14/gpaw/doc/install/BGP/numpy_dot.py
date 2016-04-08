from __future__ import print_function
num_string = "numpy"
#num_string = "Numeric"
if num_string == "numpy":
    import numpy as num
elif num_string == "Numeric":
    import Numeric as num
print(num.__file__)

from time import time

import random

N = 1700

A = num.array(num.ones((N,N)))
Al = A.tolist()
for item in Al:
    for n,value in enumerate(item):
        if (n % 2) == 0:
            item[n] = random.random()
Anew = num.array(Al)

t = time()
num.dot(Anew, Anew)
print(num_string, time()-t)
