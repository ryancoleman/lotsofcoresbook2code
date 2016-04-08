from gpaw.xc.rpa import RPACorrelation

rpa = RPACorrelation('H.ralda.pbe_wfcs.gpw',
                     txt='H.ralda_05_rpa_at_pbe.output.txt')
rpa.calculate(ecut=300)
