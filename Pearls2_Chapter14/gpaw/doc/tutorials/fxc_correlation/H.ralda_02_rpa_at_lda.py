from gpaw.xc.rpa import RPACorrelation

rpa = RPACorrelation('H.ralda.lda_wfcs.gpw',
                     txt='H.ralda_02_rpa_at_lda.output.txt')
rpa.calculate(ecut=300)
