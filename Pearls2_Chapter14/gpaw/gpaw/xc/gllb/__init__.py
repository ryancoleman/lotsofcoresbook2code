import numpy as np

def safe_sqr(u_j):
    return np.where(abs(u_j) < 1e-160, 0, u_j)**2
