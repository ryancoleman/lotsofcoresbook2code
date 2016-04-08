# libxc
export LIBXCDIR=$PWD/libxc-2.2.0/install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBXCDIR/lib

# ASE
export PYTHONPATH=$PWD/ase-3.9.0-trunk:$PYTHONPATH
export PATH=$PWD/ase-3.9.0-trunk/tools:$PATH

# GPAW setups
export GPAW_SETUP_PATH=$PWD/gpaw-setups-0.9.9672

# GPAW
export PYTHONPATH=$PWD:$PYTHONPATH
export PATH=$PWD/build/bin.linux-x86_64-intel64-2.7:$PWD/tools/:$PATH

# pyMIC
export PYMIC_LD_LIBRARY_PATH=$PYMIC_LD_LIBRARY_PATH
export PYMIC_LIBRARY_PATH=$PYMIC_LIBRARY_PATH:$PWD/gpaw/mic

# modify the prompt to show nwchem environment
if [ "x$NWCHEM_OLDPS1" == "x" ]; then
        export NWCHEM_OLDPS1="$PS1"
else
        export PS1="$NWCHEM_OLDPS1"
fi
_red="$(path tput setaf 1 2> /dev/null)"
_sgr0="$(path tput sgr0 2> /dev/null)"
export PS1="\[$_bred\]gpaw\[$_sgr0\] $PS1"
