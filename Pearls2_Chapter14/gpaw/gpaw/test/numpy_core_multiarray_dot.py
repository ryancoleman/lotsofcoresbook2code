# Segmentation fault with acml's _dotblas.so
import numpy as np
from numpy.core.multiarray import dot
b = np.ones(13, np.complex); dot(b, b)
