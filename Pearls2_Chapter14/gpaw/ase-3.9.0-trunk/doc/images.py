import os.path
from urllib import urlretrieve

url = 'http://wiki.fysik.dtu.dk/ase-files/'


def setup(app):
    pass
    
for file in ['ase/ag.png',
             'ase/dft/water_divide_surf.png',
             'ase/ase-talk.pdf']:
    if os.path.isfile(file):
        continue
    urlretrieve(url + os.path.basename(file), file)
