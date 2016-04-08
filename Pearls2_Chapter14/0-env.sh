oldpwd=`pwd`

echo Loading pyMIC environment
cd pyMIC
source ./env.sh
cd $oldpwd

echo Loading GPAW environment
cd gpaw
source ./0-env.sh
cd $oldpwd

echo The environment is all set now
