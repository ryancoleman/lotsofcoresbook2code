setenv APPS "/share/home/01067/tg803307/CAMd"
setenv MODULEFILES "${APPS}/modulefiles"

# download packages
set nose_version=0.11.3
wget http://python-nose.googlecode.com/files/nose-${nose_version}.tar.gz
set numpy_version=1.5.0
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
set ase_version=3.4.1.1765
wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase_version}.tar.gz
set gpaw_version=0.7.2.6974
wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw_version}.tar.gz
set gpaw_setups_version=0.6.6300
wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${gpaw_setups_version}.tar.gz
