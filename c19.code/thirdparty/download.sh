#!/bin/bash

# Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# get time utility from HAM
mkdir -p ham/include/ham/util/time
wget https://raw.githubusercontent.com/noma/ham/master/include/ham/util/time.hpp -P ham/include/ham/util/time/
wget https://raw.githubusercontent.com/noma/ham/master/LICENSE_1_0.txt -P ham/

# get CLU
git clone https://github.com/Computing-Language-Utility/CLU CLU.git

# get Vc
git clone git://code.compeng.uni-frankfurt.de/vc.git vc.git

# get VCL
wget http://www.agner.org/optimize/vectorclass.zip
unzip vectorclass.zip -d vcl
unzip vectorclass.zip -d vcl_mic # second copy for merging with VCLKNC
rm -f vectorclass.zip

# get VLCKNC
git clone https://bitbucket.org/veclibknc/vclknc.git vclknc.git
# merge with vcl 
cp -rf vclknc.git/* vcl_mic

