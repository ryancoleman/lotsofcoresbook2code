from gpaw.xc.rpa import RPACorrelation

#This calculation is too heavy to run as an exercise!!

rpa1 = RPACorrelation('si.rpa.isolated.gpw', txt='si.atom.rpa_output.txt')

E1_i = rpa1.calculate(ecut=400.0)
