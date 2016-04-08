#!/usr/bin/env python

import os
from optparse import OptionParser
from ase.io.trajectory import print_trajectory_info
from ase.io.bundletrajectory import print_bundletrajectory_info

description = 'Print summary of information from trajectory files.'

def main():
    p = OptionParser(usage='%prog file.traj [file2.traj ...]',
                     description=description)

    opts, args = p.parse_args()

    if len(args) == 0:
        p.error('Incorrect number of arguments')

    for f in args:
        if os.path.isfile(f):
            print_trajectory_info(f)
        elif os.path.isdir(f):
            print_bundletrajectory_info(f)
        else:
            p.error('%s is neither a file nor a directory!' % f)
