##### 
# SET UP ENVIRONMENT

source /opt/intel/composerxe/bin/compilervars.sh intel64
TOPDIR=$PWD/..

### OpenMP
OMPFLAGS="-openmp -D_REENTRANT"
OMPENABLE="--enable-openmp"

# Enable this for MIC
ARCHFLAGS="-mmic -qopt-report-phase=vec -qopt-report=1 -restrict -mGLOB_default_function_attrs=\"use_gather_scatter_hint=off\""

# If you want to compile with VTune ITT API use the line below
#  ITTFLAGS="-I/opt/intel/vtune_amplifier_xe/include"

IFFTLAGS=""
PK_CXXFLAGS=${OMPFLAGS}" -g -O2 -finline-functions -fno-alias -std=c++0x "${ARCHFLAGS}" "${ITTFLAGS}

PK_CFLAGS=${OMPFLAGS}" -g  -O2 -fno-alias -std=c99 "${ARCHFLAGS}

### Make
MAKE="make -j 8"
PK_CC=icc
PK_CXX=icpc
