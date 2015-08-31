echo
echo "-----------------------------------------------------------------"
echo " Sourcing Parallel Studio XE 2015 (Cluster Edition)"
echo "-----------------------------------------------------------------"
echo

#______________________________________________________________________
#NOTE: Update this with the appropriate script for your versions of the
#      Intel compilers and Intel MPI. Multi-threaded versions of Intel
#      MPI are REQUIRED.

source /opt/intel/parallel_studio_xe_2015/bin/psxevars.sh intel64

echo
echo "-----------------------------------------------------------------"
echo " Downloading and Extracting Uintah"
echo "-----------------------------------------------------------------"
echo

CURR="$(pwd)"

cd $HOME

if [ -d "uintah" ]; then
    rm -r uintah
    echo "Removed directory: \"\$HOME/uintah\""
    echo
fi

cd $CURR

tar xvf pearls-example.tar.gz -C $HOME

echo
echo "-----------------------------------------------------------------"
echo " Creating File Structure for 3rd Party Libraries"
echo "-----------------------------------------------------------------"
echo

cd $HOME

if [ ! -d "installs" ]; then
    mkdir installs
    echo "Created directory: \"\$HOME/installs\""
    echo
fi

cd $HOME/installs

if [ ! -d "tarballs" ]; then
    mkdir tarballs
    echo "Created directory: \"\$HOME/installs/tarballs\""
    echo
fi

cd $HOME/installs/tarballs

echo
echo "-----------------------------------------------------------------"
echo " Downloading and Extracting libxml2-2.7.8"
echo "-----------------------------------------------------------------"
echo

if [ -f "libxml2-2.7.8.tar.gz" ]; then
    rm libxml2-2.7.8.tar.gz
    echo "Removed file: \"\$HOME/installs/tarballs/libxml2-2.7.8.tar.gz\""
    echo
fi

if [ -d "libxml2-2.7.8" ]; then
    rm -r libxml2-2.7.8
    echo "Removed directory: \"\$HOME/installs/tarballs/libxml2-2.7.8\""
    echo
fi

wget http://xmlsoft.org/sources/libxml2-2.7.8.tar.gz
tar xvf libxml2-2.7.8.tar.gz

echo
echo "-----------------------------------------------------------------"
echo " Downloading and Extracting zlib-1.2.8"
echo "-----------------------------------------------------------------"
echo

if [ -f "zlib-1.2.8.tar.gz" ]; then
    rm zlib-1.2.8.tar.gz
    echo "Removed file: \"\$HOME/installs/tarballs/zlib-1.2.8.tar.gz\""
    echo
fi

if [ -d "zlib-1.2.8" ]; then
    rm -r zlib-1.2.8
    echo "Removed directory: \"\$HOME/installs/tarballs/zlib-1.2.8\""
    echo
fi

wget http://zlib.net/zlib-1.2.8.tar.gz
tar xvf zlib-1.2.8.tar.gz

echo
echo "-----------------------------------------------------------------"
echo " Installing zlib-1.2.8 (coprocessor-build)"
echo "-----------------------------------------------------------------"
echo

if [ -d "$HOME/installs/mic/zlib-1.2.8" ]; then
    rm -r $HOME/installs/mic/zlib-1.2.8
    echo "Removed directory: \"\$HOME/installs/mic/zlib-1.2.8\""
    echo
fi

cd $HOME/installs/tarballs/zlib-1.2.8
make clean
make cleanreally
make reallyclean
CC=icc CXX=icpc CFLAGS="-mmic" CXXFLAGS="-mmic" LDFLAGS="-mmic" \
    ./configure \
    --prefix=$HOME/installs/mic/zlib-1.2.8 \
    --static
make -j
make install

export LD_LIBRARY_PATH=$HOME/installs/mic/zlib-1.2.8/lib:$LD_LIBRARY_PATH

echo 
echo "-----------------------------------------------------------------"
echo " Installing libxml2-2.7.8 (coprocessor-build)"
echo "-----------------------------------------------------------------"
echo

if [ -d "$HOME/installs/mic/libxml2-2.7.8" ]; then
    rm -r $HOME/installs/mic/libxml2-2.7.8
    echo "Removed directory: \"\$HOME/installs/mic/libxml2-2.7.8\""
    echo
fi

cd $HOME/installs/tarballs/libxml2-2.7.8
./configure \
    --prefix=$HOME/installs/mic/libxml2-2.7.8 \
    --host=x86_64-k1om-linux \
    --enable-static \
    --without-python \
    CC=icc \
    CXX=icpc \
    CFLAGS="-mmic" \
    CXXFLAGS="-mmic" \
    LDFLAGS="-mmic"
make -j
make install

export LD_LIBRARY_PATH=$HOME/installs/mic/libxml2-2.7.8/lib:$LD_LIBRARY_PATH
export PATH=$HOME/installs/mic/libxml2-2.7.8/bin:$PATH

echo
echo "-----------------------------------------------------------------"
echo " Installing Uintah"
echo "-----------------------------------------------------------------"
echo

cd $HOME/uintah/branches/pearls

if [ ! -d "example" ]; then
    mkdir example
    echo
    echo "Created directory: \"example\""
fi

cd example

#______________________________________________________________________
#NOTE: Update this with the appropriate paths for your versions of the
#      Intel compilers and Intel MPI. Multi-threaded versions of Intel
#      MPI are REQUIRED.

../src/configure \
    --host=x86_64-k1om-linux \
    --enable-64bit \
    --enable-optimize="-O2 -mmic -mt_mpi" \
    --enable-assertion-level=0 \
    --enable-static \
    --with-libxml2=$HOME/installs/mic/libxml2-2.7.8 \
    --with-mpi=/opt/intel/impi_5.0.1/mic \
    --with-zlib=$HOME/installs/mic/zlib-1.2.8 \
    CC=mpiicc \
    CXX=mpiicpc \
    F77=mpiifort

make -j sus

cd $HOME/uintah/branches/pearls

if [ ! -d "tmp.example" ]; then
    mkdir tmp.example
    echo "Created directory: \"tmp.example\""
    echo
fi

cd example/StandAlone
cp sus $HOME/uintah/branches/pearls/tmp.example/sus.mic-example

echo
echo "-----------------------------------------------------------------"
echo " Transferring Intel Libraries to mic0"
echo "-----------------------------------------------------------------"
echo

ssh mic0 "mkdir -p /$HOME/intel/bin"
ssh mic0 "mkdir -p /$HOME/intel/lib64"

#______________________________________________________________________
#NOTE: Update this with the appropriate paths for your versions of the
#      Intel compilers and Intel MPI. Multi-threaded versions of Intel
#      MPI are REQUIRED.

scp /opt/intel/impi_5.0.1/mic/bin/* mic0:$HOME/intel/bin/
scp /opt/intel/impi_5.0.1/mic/lib/libmpicxx.* mic0:$HOME/intel/lib64/
scp /opt/intel/impi_5.0.1/mic/lib/libmpifort.* mic0:$HOME/intel/lib64/
scp /opt/intel/impi_5.0.1/mic/lib/libmpi_mt.* mic0:$HOME/intel/lib64/
scp /opt/intel/composer_xe_2015.0.090/compiler/lib/mic/* mic0:$HOME/intel/lib64/

echo
echo "-----------------------------------------------------------------"
echo " Transferring UCF-Related Files to mic0"
echo "-----------------------------------------------------------------"
echo

ssh mic0 "mkdir -p /$HOME/uintah/branches/pearls/src/StandAlone/inputs"

scp -r $HOME/uintah/branches/pearls/src/StandAlone/inputs/UPS_SPEC/ mic0:$HOME/uintah/branches/pearls/src/StandAlone/inputs/

scp -r $HOME/uintah/branches/pearls/tmp.example/ mic0:$HOME

echo
echo "-----------------------------------------------------------------"
echo " Going Down Successfully"
echo "-----------------------------------------------------------------"
echo
echo " To launch your first Uintah Computational Framework simulations:"
echo "    ssh mic0"
echo '    export LD_LIBRARY_PATH=$HOME/intel/lib64:$LD_LIBRARY_PATH'
echo '    export PATH=$HOME/intel/bin:$PATH'
echo "    cd tmp.example"
echo "    ./run-pearls-example.sh"
echo "-----------------------------------------------------------------"
echo
