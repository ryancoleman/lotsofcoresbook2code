import numpy as np
from gpaw.response.bse import BSE

w_grid = np.linspace(0, 15, 1001)

# It stores the four-points kernel used for building the two-particles
# Hamiltonian in LiF_W_qGG.
bse = BSE('LiF_fulldiag.gpw',
          w=w_grid,
          q=np.array([0.0001, 0., 0.]),
          optical_limit=True,
          ecut=30,
          nbands=60,
          eta=0.1,
          kernel_file='LiF_W_qGG',
          txt='LiF_BSE_out.txt')

# Calculate the dielectric function calculated at the BSE level:
df_BSE = bse.get_dielectric_function()
