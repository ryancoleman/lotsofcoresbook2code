#!/usr/bin/python
"""icc.py is a wrapper for the Intel compiler,
   converting/removing incompatible gcc args.   """

import sys
from subprocess import call

args2change = {"-fno-strict-aliasing":"",
               "-fmessage-length=0":"",
               "-Wall":"",
               "-std=c99":"-qlanglvl=extc99",
               "-fPIC":"",
               "-g":"",
               "-D_FORTIFY_SOURCE=2":"",
               "-DNDEBUG":"",
               "-UNDEBUG":"",
               "-pthread":"",
               "-shared":"-qmkshrobj",
               "-Xlinker":"",
               "-export-dynamic":"",
               "-Wstrict-prototypes":"",
               "-dynamic":"",
               "-O3":"",
               "-O3":"",
               "-O2":"",
               "-O1":""}

fragile_files = ["test.c"]

cmd = ""
fragile = False
for arg in sys.argv[1:]:
    cmd += " "
    t = arg.strip()
    if t in fragile_files:
        fragile = True
    if t in args2change:
        cmd += args2change[t]
    else:
        cmd += arg

flags = "-w -O3 -std=c99"
cmd = "mpicc %s %s"%(flags, cmd)

call(cmd, shell=True)

