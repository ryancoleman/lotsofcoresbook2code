export APPS="~/CAMd"
export MODULEFILES="${APPS}/modulefiles"

cd ${APPS}
# download packages
nose_version=1.1.2
http://pypi.python.org/packages/source/n/nose/nose-${nose_version}.tar.gz
numpy_version=1.5.1
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
ase_version=3.6.0.2515
wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase_version}.tar.gz
gpaw_version=0.9.0.8965
wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw_version}.tar.gz
gpaw_setups_version=0.8.7929
wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${gpaw_setups_version}.tar.gz
