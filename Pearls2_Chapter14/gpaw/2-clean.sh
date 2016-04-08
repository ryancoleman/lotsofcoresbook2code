#!/bin/bash

python setup.py clean

basedir=$PWD

cd $basedir/libxc-2.2.0
if [ -e libxc.build.done ]; then
	make clean
	rm -f libxc.build.done
fi
cd $basedir

