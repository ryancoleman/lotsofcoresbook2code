# Location of cross compiled libxml2
#
# Set this so we can add it to the front of the path and to the 
# --with-libxml2
#
LIBXML=/home/ed/bj/Devel/SciDAC/install/libxml2-2.6.6

# Location of the QOS library 
# This is used to find the setup scripts and to set the 
# QMP include and link paths
QOS=/home/ed/bj/v2.5.3/aix5.2

# Location of the crosscompile tools and extra packages
CROSSCOMPILE_TOOLS=/qcdoc/sfw/packages/v1.1/local/bin

#### Hacks to make things work:
#
# Create a directory, create a symlink to powerpc-gnu-elf-ar
# to work around autotools stupidity of not letting me set it.
# I have to add it to the front of my path so that it gets picked up
# first

export HEREDIR=`pwd`
if [ ! -d ./build ]; 
then 
   mkdir build;
fi
if [ ! -e ./build/ar ];
then 
   echo Creating Symlink to ${CROSSCOMPILE_TOOLS} powerpc-gnu-elf-ar
   ln -s ${CROSSCOMPILE_TOOLS}/powerpc-gnu-elf-ar ./build/ar
fi

echo sourcing setup scripts
source ${QOS}/scripts/setup.sh

echo prepending ${LIBXML}/bin to path so native one doesnt get picked up
export PATH=${LIBXML}/bin:$PATH

echo prepending full path of ./build to path for the ar override
export PATH=${HEREDIR}/build:${PATH}
