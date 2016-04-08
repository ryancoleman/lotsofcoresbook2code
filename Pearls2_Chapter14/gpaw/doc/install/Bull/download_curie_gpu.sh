export APPS="readlink -f ~/CAMd"
export MODULEFILES="${APPS}/modulefiles"

# warning - firewall blocks that, so download on other machine and scp!

cd ${APPS}
# download packages
nose_version=1.1.2
wget http://pypi.python.org/packages/source/n/nose/nose-${nose_version}.tar.gz
numpy_version=1.5.1
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
scipy_version=0.9.0
wget https://downloads.sourceforge.net/project/scipy/scipy/${scipy_version}/scipy-${scipy_version}.tar.gz
ase_version=3.7.0.3168
wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase_version}.tar.gz
gpaw_version=0.9.0.8965
wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw_version}.tar.gz
gpaw_setups_version=0.9.9672
wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${gpaw_setups_version}.tar.gz
# OK, curie does not allow svn!
svn co https://svn.fysik.dtu.dk/projects/gpaw/branches/rpa-gpu-expt
