from gpaw.response.df import DielectricFunction as DF

df = DF('LiF_fulldiag.gpw',
        domega0=0.01,  # grid-spacing at 0 eV
        omega2=10.0,   # frequency where grid-spacing has doubled to 0.02 eV
        eta=0.1,       # broadening parameter
        nbands=60,     # number of bands to consider for building chi
        ecut=30,       # energy cutoff for planewaves
        txt='LiF_RPA_out2.txt')

# Calculate the dielectric function without and with local field effects:
df.get_dielectric_function()
