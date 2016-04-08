from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw.atom.basis import BasisMaker

symbol = 'Au'
args = parameters[symbol] # Dictionary of default setup parameters
args['rcut'] = 2.6 # Set cutoff of augmentation sphere

generator = Generator(symbol, 'RPBE')
generator.N *= 2 # Increase grid resolution
generator.run(write_xml=False, **args)

bm = BasisMaker(generator, name='special', run=False)

# Create double-zeta RPBE basis where p orbital is considered a valence state
# (ordinary dzp basis would use a smaller p-type Gaussian for polarization)
# The list jvalues indicates which states should be included, and the
# ordering corresponds to the valence states in the setup.
basis = bm.generate(zetacount=2, polarizationcount=0,
                    energysplit=0.1, jvalues=[0, 1, 2],
                    rcutmax=12.0)

basis.write_xml() # Dump to file 'Au.special.basis'
