#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Emacs: treat this as -*- python -*-
from __future__ import print_function
from optparse import OptionParser

parser = OptionParser(usage='%prog [options]',
                      version='%prog 0.1')
parser.add_option('--dir', dest="dir",
                  default='.',
                  help='Results directory')
parser.add_option("--runs", dest="runs",
                  default=5,
                  help='use that many runs to calculate the average.')
parser.add_option("--startcores", dest="startcores",
                  default=1,
                  help='use at lease that many cores.')

opt, args = parser.parse_args()

import os

import datetime

from math import sqrt

import numpy as np

colors = [
    'black',
    'brown',
    'red',
    'orange',
    'yellow',
    'green',
    'blue',
    'violet',
    'gray',
    'gray']

from ase.data import atomic_numbers as numbers

def plot(xdata, ydata, std,
         title,
         xlabel, ylabel,
         label, color, num=1):
    import matplotlib
    #matplotlib.use('Agg')
    import pylab

    # all goes to figure num
    pylab.figure(num=num, figsize=(7, 5.5))
    pylab.gca().set_position([0.10, 0.20, 0.85, 0.60])
    # let the plot have fixed y-axis scale
    miny = min(ydata)
    maxy = max(ydata)
    ywindow = maxy - miny
    pylab.gca().set_ylim(miny-ywindow/4.0, maxy+ywindow/3.0)
    #pylab.plot(xdata, ydata, 'b.', label=label, color=color)
    #pylab.plot(xdata, ydata, 'b-', label='_nolegend_', color=color)
    pylab.bar(xdata, ydata, 0.3, yerr = std, label=label, color=color)
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    #pylab.legend(loc='upper right')
    #pylab.savefig(directory_name + os.path.sep + out_prefix +'.png')

def plot_save(directory_name, out_prefix):
    from os.path import exists
    assert exists(directory_name)
    import pylab

    pylab.savefig(directory_name + os.path.sep + out_prefix +'.png')

def analyse_benchmark(ncores=8, startcores=1, machine='TEST', runs=7):
    #system = ['carbon_py']
    #system = ['carbon']
    #system = ['niflheim_py']
    #system = ['niflheim']
    #system = ['TEST_py']
    system = machine+'_py'

    systems_string = {
        'carbon_py' : 'gpaw 1865 on carbon',
        'carbon' : 'mkl 10.0.2.018 dsyev on carbon',
        'niflheim_py' : 'gpaw 1865 on niflheim',
        'niflheim' : 'acml 4.0.1 dsyev on niflheim',
        #'TEST_py' : 'gpaw on TEST',
        }.get(system, False)

    processes = {
        'carbon_py' : [1, 2, 4, 6, 8],
        'carbon' : [1, 2, 4, 8],
        'niflheim_py' : [1, 2, 3, 4],
        'niflheim' : [1, 2, 4],
        #'TEST_py' : [1, 2, 4, 6, 8],
        }.get(system, False)


    if not systems_string:
        systems_string = 'gpaw on '+machine
    if not processes:
        processes = [startcores]
        for n in range(startcores+1, ncores+1):
            if n%2==0:
                processes.append(n)

    timer_entries_all = []
    if system.find('_py') == -1:
        for i in range(runs):
            timer_entries_all.append('run: '+str(i))
    else:
        for i in range(runs):
            timer_entries_all.append('Run:  '+str(i))

    import re

    # Select timer entries
    selected_entries = range(runs)

    height = {}

    gpaw_versions = []

    pre_results = {}
    results = {}

    timer_entries = []
    timer_entries_re = {}
    for entry in selected_entries:
        height[entry] = []
        timer_entries.append(timer_entries_all[entry])
        timer_entries_re[timer_entries_all[entry]] = re.compile(timer_entries_all[entry])

    # absolute path to directory
    root_abspath = os.path.abspath(opt.dir)
    # lenght of directory name
    rootlen = len(root_abspath) + 1

    ref_value_3300 = -44.85826
    ref_SCF_3300 = 19
    ref_value_3301 = -44.85709
    ref_SCF_3301 = 31
    ref_value_3721 = -44.85666
    ref_SCF_3721 = 35
    ref_value_5147 = -44.83504
    ref_SCF_5147 = 30

    ref_value_6383 = -44.84197
    ref_SCF_6383 = 28

    ref_value = ref_value_6383
    ref_SCF = ref_SCF_6383

    tolerance = 0.0001

    ref_failed = False
    h_failed = False
    for run in [str(p)+'_01' for p in processes]:
        # extract results
        rundir = os.path.join(root_abspath, system+run)
        file = os.path.join(rundir, 'out.txt')
        try:
            f = open(file, 'r')
            #
            print('Analysing '+file, end=' ')
            #
            lines = f.readlines()
        except: pass
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
        if len(str(gpaw_version)) > 1:
            if gpaw_version <= 6383:
                ref_value = ref_value_6383
                ref_SCF = ref_SCF_6383
            if gpaw_version <= 5147:
                ref_value = ref_value_5147
                ref_SCF = ref_SCF_5147
            if gpaw_version <= 3720:
                ref_value = ref_value_3301
                ref_SCF = ref_SCF_3301
            if gpaw_version <= 3300:
                ref_value = ref_value_3300
                ref_SCF = ref_SCF_3300
        elif len(str(gpaw_version)) == 1:
            if gpaw_version <= 4:
                ref_value = ref_value_3300
                ref_SCF = ref_SCF_3300
        gpaw_versions.append(gpaw_version)
        # search for timings
        print('gpaw version %d' % gpaw_version)
        for entry in selected_entries:
            h = []
            ref = []
            for line in lines:
                m = timer_entries_re[timer_entries_all[entry]].search(line)
                if m is not None:
                    h.append(float(line.split(':')[-1]))
                   #break # stop after the first match
            for h_entry in h:
                if float(h_entry) < 0.0:
                    h_failed = True
                    break
            height[entry].append(h)
            for line in lines:
                m = re.compile('Zero').search(line)
                if m is not None:
                    ref.append(float(line.split(':')[-1]))
                   #break # stop after the first match
            for ref_entry in ref:
                if abs(float(ref_entry)-ref_value) > tolerance:
                    ref_failed = True
                    break
    #
    assert len(processes) == len(gpaw_versions)
    for p in range(len(processes)):
        assert gpaw_versions[p] == max(gpaw_versions), 'incompatible gpaw versions across cores'
    #
    if h_failed:
        print('Panic: negative time in '+file)
        assert not h_failed
    if ref_failed:
        print('Panic: wrong Zero Kelvin: value in '+file+' - should be '+str(ref_value)+' +- '+str(tolerance))
        assert not ref_failed
    # arrange results
    for p in range(len(processes)):
        pre_results[processes[p]] = []
        for i in range(len(height)):
            pre_results[processes[p]].append(height[i][p])
    #
    # arrange results - calculate statistics
    for p in processes:
    #for p in range(len([1])):
        #print pre_results[p]
        results[p] = []
        temp = []
        for q in range(p):
            temp_q = []
            for i in range(len(pre_results[p])):
                #print pre_results[p][i][q]
                temp_q.append(pre_results[p][i][q])
                temp.append(pre_results[p][i][q])
            # averages for a given core q
            results[p].append((np.average(temp_q), np.std(temp_q)))
        # max, avrg, and std across all cores
        results[p].append((np.average(temp), np.std(temp), min(temp), max(temp)))
    #for p in processes:
    #    #N = len(pre_results[p])
    #    #avg = sum(pre_results[p])/N
    #    #q = sqrt(sum([(x-avg)**2/(N) for x in pre_results[p]]))
    #    avg.append(np.average(pre_results[p]))
    #    q.append(np.std(pre_results[p]))
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pylab, ticker
    from twiny import twiny
    # from http://matplotlib.sourceforge.net/examples/dashtick.py
    ROTATION=75
    DASHROTATION=115
    DASHBASE=5
    DASHLEN=25
    DASHSTAGGER=3
    FONTSIZE=10
    def dashlen(step):
        return DASHBASE+(DASHLEN*(step%DASHSTAGGER))
    # print scaling results
    parameters = processes
    zero = [0.0 for i in range(len(parameters))]
    pylab.plot(parameters, zero, 'k-', label='_nolegend_')
    ay1=pylab.gca()
    ay1.xaxis.set_ticks(parameters)
    ay1.xaxis.set_ticklabels([str(x) for x in parameters])
    for p in processes:
        parameters = []
        avg = []
        std = []
        for i in range(len(results[p])-1):
            parameters.append(p+0.3*i)
            # avg and std across processes
            avg.append(results[p][i][0])
            std.append(results[p][i][1])
        # height
        #print parameters, avg, std
        print('No. of processes '+str(int(parameters[0]))+': time [sec]: avg '+str(round(results[p][-1][0],1))+', stddev '+str(round(results[p][-1][1],1))+', min '+str(round(results[p][-1][2],1))+', max '+str(round(results[p][-1][3],1)))
        plot(
            parameters, avg, std,
            systems_string+' version '+str(gpaw_version),
            'processes per node',
            'time [s]',
            'gpaw',
            (colors[p%10]),
            num=1)
    # from two_scales.py
    plot_save(".", 'memory_bandwidth_'+system)
    pylab.close(1)
#

if __name__ == '__main__':
    from os import environ

    NCORES = int(environ.get('NCORES', 8))
    MACHINE = environ.get('MACHINE', 'TEST')
    assert NCORES >= 1, str(NCORES)+' must be >= 1'

    runs = int(opt.runs)
    assert runs >= 1, runs+' must be >= 1'
    startcores = int(opt.startcores)
    assert startcores >= 1, startcores+' must be >= 1'

    analyse_benchmark(NCORES, startcores, MACHINE, runs=runs)
