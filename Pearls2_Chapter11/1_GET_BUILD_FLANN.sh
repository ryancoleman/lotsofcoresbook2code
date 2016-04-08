#set path to icc and icpc

#get flann and unzip it
pushd ..
wget http://www.cs.ubc.ca/research/flann/uploads/FLANN/flann-1.8.4-src.zip
unzip flann-1.8.4-src.zip
cp -r flann-1.8.4-src flann-1.8.4-src-mic

cd flann-1.8.4-src
mkdir build
cd build
cmake ../ -DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS=" -O3" -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS=" -O3 "
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
cd ../..

cd flann-1.8.4-src-mic
mkdir build
cd build
cmake ../ -DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS=" -O3 -mmic" -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS=" -O3 -mmic "
make
export MIC_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH:$PWD/lib
cd ../..

popd





