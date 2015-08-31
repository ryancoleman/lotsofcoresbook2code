#!/bin/bash

# Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

scp -r bin.mic mic0:
ssh mic0 env LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH:$LIBRARY_PATH KMP_PLACE_THREADS=61c,4t OMP_NUM_THREADS=244 KMP_AFFINITY=granularity=fine,balanced KMP_BLOCK_TIME=0 $1

