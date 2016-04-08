#!/bin/bash

basedir=$PWD

rm $basedir/build.log
cd $basedir/libxc-2.2.0
if [ ! -e libxc.build.done ]; then
	./configure CC=icc CFLAGS="-fPIC -xHOST" --prefix=$PWD/install -enable-shared |& tee $basedir/build.log
	(make -j 10 && make install && touch libxc.build.done) |& tee -a $basedir/build.log
fi
cd $basedir
 
python setup.py build_ext |& tee -a $basedir/build.log
cd $basedir/gpaw/mic
make |& tee -a $basedir/build.log
cd $basedir 

