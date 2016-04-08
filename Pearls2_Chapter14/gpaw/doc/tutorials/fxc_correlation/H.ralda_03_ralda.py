from gpaw.xc.fxc import FXCCorrelation

fxc = FXCCorrelation('H.ralda.lda_wfcs.gpw',
                     xc='rALDA', txt='H.ralda_03_ralda.output.txt')
fxc.calculate(ecut=300)
