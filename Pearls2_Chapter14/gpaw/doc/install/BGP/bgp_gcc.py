#!/usr/bin/env python
"""bgp_gcc.py is a wrapper for the BGP gcc compiler,
   converting/removing incompatible gcc args.   """
from __future__ import print_function
import sys
from subprocess import call
from glob import glob

args2change = {"-fno-strict-aliasing":"",
               "-fmessage-length=0":"",
               "-Wall":"",
               "-std=c99":"",
               "-fPIC":"",
               "-g":"",
               "-D_FORTIFY_SOURCE=2":"",
               "-DNDEBUG":"",
               "-UNDEBUG":"",
               "-pthread":"",
               "-shared":"",
               "-Xlinker":"",
               "-export-dynamic":"",
               "-Wstrict-prototypes":"",
               "-dynamic":"",
               "-O3":"",
               "-O2":"",
               "-O1":""}

fragile_files = ["test.c"]
non_c99files = glob('c/libxc/src/*.c')

cmd = ""
opt = 1

for arg in sys.argv[1:]:
    cmd += " "
    t = arg.strip()
    if t in fragile_files:
        opt = 2
    if t in non_c99files:
        opt = 3
    if t in args2change:
        cmd += args2change[t]
    else:
        cmd += arg

flags_list = {1: "-g -O3 -std=c99 -fPIC",
              2: "-g -O2 -std=c99 -fPIC",
              3: "-g -O3 -fPIC",
              }

flags = flags_list[opt]     
cmd = "mpicc %s %s"%(flags, cmd)

print("\nexecmd: %s\n"%cmd)
call(cmd, shell=True)
