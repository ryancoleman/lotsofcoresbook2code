#!/bin/bash

export KMP_AFFINITY=compact
export OMP_NUM_THREADS=56
bin/MD -s 4
