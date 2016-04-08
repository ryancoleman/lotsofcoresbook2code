code = open('eels.py').read().replace('Ag_GLLBSC.gpw',
                                      '../band_structure/Ag_GLLBSC.gpw')
with open('eels2.py', 'w') as fd:
    fd.write(code)

execfile('eels2.py')
