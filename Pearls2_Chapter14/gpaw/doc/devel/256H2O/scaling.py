#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Emacs: treat this as -*- python -*-
from __future__ import print_function
from optparse import OptionParser

parser = OptionParser(usage='%prog [options] output_prefix.\nExample of call:\n'+
                      'python %prog --dir=. --pattern="b256H2O_120_04x04m64.grid_*_" b256H2O\n',
                      version='%prog 0.1')
parser.add_option('--dir', dest="dir",
                  default='.',
                  help='results directory.')
parser.add_option('--iter', dest="iter",
                  default=5,
                  help='number of SCF steps.')
parser.add_option('--pattern', dest="pattern",
                  default='',
                  help='pattern for directories to search.')
parser.add_option('-v', '--verbose', action='store_true',
                  default=False,
                  help='verbose mode.')

opt, args = parser.parse_args()

import os

def T(s):
    """Return time in seconds from hh:mm:ss"""
    t = 0
    for x in s.split(':'):
        t = t * 60 + int(x)
    return t

def analyse_benchmark(dir, pattern, output_prefix, iter, verbose=False):

    from os.path import abspath, exists, join

    output = output_prefix+'.txt'

    root_abspath = abspath(dir)
    # length of directory name
    rootlen = len(root_abspath) + 1

    assert exists(root_abspath)

    from glob import glob

    f_list = glob(join(root_abspath, pattern, output))
    f_list.sort()
    assert not f_list == [], 'Error: list of output files is empty'

    import re

    time = {}
    processes = []

    for f in f_list:
        # extract the number of processes p
        d = f.split('/')[-2]
        p = d.split('_')[-2]
        p = int(p)
        processes.append(p)
        #
        time[p] = {
        'start': 0.0,
        'fixdensity_start_estimate': 0.0,
        'fixdensity_end': 0.0,
        'SCF_start': 0.0,
        'SCF_end': 0.0,
        'forces_start': 0.0,
        'forces_end': 0.0,
        'total': 0.0,
        }
        #
        lines = open(f).readlines()
        # extract gpaw version
        for n, l in enumerate(lines):
            if l.startswith(' |__ |  _|___|_____|'):
                gpaw_version = lines[n + 0].strip().split()[3].split('.')[-1]
                break
        if gpaw_version[-1] == 'M':
            gpaw_version = gpaw_version[:-1]
        if gpaw_version.rfind(':') != -1:
            gpaw_version = gpaw_version[:gpaw_version.rfind(':')]
        gpaw_version = int(gpaw_version)
        # assume old version (< 3172)
        start_iter = 0
        # old style
        #
        #            Atomic orbitals used for initialization: 1536
        #                     log10-error:    Total        Iterations:
        #           Time      WFS    Density  Energy       Fermi  Poisson
        #iter:   0  12:10:09                  -3568.07860  0      11
        #                     log10-error:    Total        Iterations:
        #           Time      WFS    Density  Energy       Fermi  Poisson
        #iter:   0  12:13:29  +0.2            -3663.00195  0      1
        #iter:   1  12:15:18  -0.8            -4417.57264  0      1
        #iter:   2  12:17:07  -1.3            -4469.68829  0      1
        #iter:   3  12:18:57  -0.9   -0.8     -4091.42827  0      7
        #iter:   4  12:20:48  -1.1   -1.0     -4055.26110  0      7
        #iter:   5  12:22:40  -1.4   -1.3     -4106.38102  0      7
        #
        # new style
        #                     log10-error:    Total        Iterations:
        #           Time      WFS    Density  Energy       Fermi  Poisson
        #iter:   1  16:16:59  +0.4            -3663.00195  0      1
        #iter:   2  16:19:11  -0.8            -4417.57264  0      1
        #iter:   3  16:21:21  -1.3            -4469.68829  0      1
        #iter:   4  16:23:34  -0.9   -0.8     -4091.42827  0      7
        #iter:   5  16:25:49  -1.1   -1.0     -4055.26110  1      7
        #
        if len(str(gpaw_version)) == 1:
            # stable release found
            # assume new version - this is wrong but nothing else can be done
            start_iter = 1
        if len(str(gpaw_version)) > 4:
            # more then 4 digits in svnversion found
            # assume new version (can't compare strings here)
            start_iter = 1
        else:
            if gpaw_version >= 3172:
                # new version
                start_iter = 1
        # extract start time
        for n, l in enumerate(lines):
            if l.startswith('Date: '):
                #print l, n, f
                t = T(lines[n + 0].split()[4])
                break
        time[p]['start'] = t
        # extract SCF begining time estimate and end time (constant potential steps (fixdensity))
        for n, l in enumerate(lines):
            if l.startswith('iter:   1'):
                #print l, n, f
                if start_iter == 0:
                    fixdensity_start = n-1
                else:
                    fixdensity_start = n+0
                t1 = T(lines[fixdensity_start + 0].split()[2])
                t2 = T(lines[fixdensity_start + 1].split()[2])
                t3 = T(lines[fixdensity_start + 2].split()[2])
                break
        # estimate the begining of fixdensity based on 3 fixdensity steps
        time[p]['fixdensity_start_estimate'] = t1-(t3-t1)/2.0
        time[p]['fixdensity_end'] = t3
        # extract SCF begining and end time
        time[p]['SCF_start'] = time[p]['fixdensity_end']
        for n, l in enumerate(lines):
            if l.startswith('iter: '+"%3s" % iter):
                #print l, n, f
                t = T(lines[n + 0].split()[2])
                break
        time[p]['SCF_end'] = t
        for n, l in enumerate(lines):
            if l.startswith('Total:') and (l.find('MB') == -1) and (l.find('GB') == -1):
                #print l, n, f
                t = l.split()[1]
                break
        time[p]['total'] = t
    #
    # results
    #
    speedup = {}
    efficiency = {}
    #
    if verbose:
        print("# p - processes, p0 - reference processes, t - time [sec], s - speedup, e - efficiency")
        print("# GPAW version "+str(gpaw_version)+": stages: 1 - initialization, 2 - fixdensity, 3 - SCF, 4 - forces, 5 - total")
        print("# p    "+" p/p0  "+" t1     "+" s1     "+" e1   "+" t2     "+" s2     "+" e2   "+" t3     "+" s3     "+" e3   "+" t4     "+" s4     "+" e4   "+" t5     "+" s5     "+" e5")
    for p in processes:
        time[p]['init'] = time[p]['fixdensity_start_estimate'] - time[p]['start']
        time[p]['fixdensity'] = time[p]['fixdensity_end'] - time[p]['fixdensity_start_estimate']
        time[p]['SCF'] = time[p]['SCF_end'] - time[p]['SCF_start']
        time[p]['forces'] = time[p]['forces_end'] - time[p]['forces_start']
        tot = max(time[p]['fixdensity_end'], time[p]['SCF_end'], time[p]['forces_end'])-time[p]['start']
        sum_of_entries = time[p]['init'] + time[p]['fixdensity']+ time[p]['SCF'] + time[p]['forces']
        #print time[p]['init'], time[p]['fixdensity'], time[p]['SCF'], time[p]['forces']
        if verbose:
            if abs(float(tot)-float(time[p]['total'])) > 5.0:
                print('Warning: Sum of time entries: '+str(tot)+' does not match total time in the output: '+str(time[p]['total']))
        time[p]['total'] = sum_of_entries
        # calculate
        speedup[p] = {}
        efficiency[p] = {}
        for stage in ['init', 'fixdensity', 'SCF', 'forces', 'total']:
            if time[p][stage] > 0.0:
                speedup[p][stage] = float(time[processes[0]][stage])/float(time[p][stage])*processes[0]
            else:
                speedup[p][stage] = 0.0
            efficiency[p][stage] = speedup[p][stage]/p
        # print results
        print('  %5d %6.2f %7.1f %7.1f %5.2f %7.1f %7.1f %5.2f %7.1f %7.1f %5.2f %7.1f %7.1f %5.2f %7.1f %7.1f %5.2f' %(
            p, float(p)/processes[0],
            float(time[p]['init']), speedup[p]['init'], efficiency[p]['init'],
            float(time[p]['fixdensity']), speedup[p]['fixdensity'], efficiency[p]['fixdensity'],
            float(time[p]['SCF']), speedup[p]['SCF'], efficiency[p]['SCF'],
            float(time[p]['forces']), speedup[p]['forces'], efficiency[p]['forces'],
            float(time[p]['total']), speedup[p]['total'], efficiency[p]['total'],
            ))

if __name__ == '__main__':
    from os import environ

    assert len(args) == 1, 'Error: Only one argument allowed: output prefix'
    output_prefix = args[0]

    analyse_benchmark(opt.dir, opt.pattern, output_prefix, opt.iter, opt.verbose)
