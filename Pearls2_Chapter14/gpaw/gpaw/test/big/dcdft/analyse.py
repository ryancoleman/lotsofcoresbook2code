import os

import sys

import urllib
import urllib2

import tarfile
import zipfile

import csv

import numpy as np

from ase.test import NotAvailable

from ase import units

from ase.test.tasks.dcdft import DeltaCodesDFTTask as Task

dir = 'Delta'

if len(sys.argv) == 1:
    tag = None
    reffile = os.path.join(dir, 'WIEN2k.txt')
else:
    if len(sys.argv) == 3:
        tag = sys.argv[1]
        reffile = sys.argv[2]
    else:
        tag = sys.argv[1]
        reffile = os.path.join(dir, 'WIEN2k.txt')

src = 'https://molmod.ugent.be/sites/default/files/Delta_v3-0.zip'
name = os.path.basename(src)
if not os.path.exists(dir): os.makedirs(dir)
os.chdir(dir)
if not os.path.exists('calcDelta.py'):
    try:
        resp = urllib2.urlopen(src)
        urllib.urlretrieve(src, filename=name)
        z = zipfile.ZipFile(name)
        try:  # new in 2.6
            z.extractall()
        except AttributeError:
            # http://stackoverflow.com/questions/7806563/how-to-unzip-a-zip-file-with-python-2-4
            for f in z.namelist():
                fd = open(f, "w")
                fd.write(z.read(f))
                fd.close()
        # AttributeError if unzip not found
    except (urllib2.HTTPError, AttributeError):
        raise NotAvailable('Retrieval of zip failed')
os.chdir('..')

task = Task(
    tag=tag,
    use_lock_files=True,
    )

# header
h = ['#element', 'V0', 'B0', 'B1']

if not os.path.exists('%s_raw.csv' % tag):
    # read calculated results from json file and write into csv
    task.read()
    task.analyse()

    f1 = open('%s_raw.csv' % tag, 'wb')
    csvwriter1 = csv.writer(f1)
    csvwriter1.writerow(h)

    for n in task.collection.names:
        row = [n]
        if n in task.data.keys():
            try:
                v = task.data[n]['dcdft volume']
                b0 = task.data[n]['dcdft B0'] / (units.kJ * 1e-24)
                b1 = task.data[n]['dcdft B1']
                row.extend([v, b0, b1])
            except KeyError: # completely failed to find eos minimum
                row.extend(['N/A', 'N/A', 'N/A'])
        else:
            # element not calculated
            row.extend(['N/A', 'N/A', 'N/A'])
        if 'N/A' not in row:
            csvwriter1.writerow(row)
    f1.close()

# read raw results
csvreader1 = csv.reader(open('%s_raw.csv' % tag, 'r'))
data = {}
for row in csvreader1:
    if '#' not in row[0]:
        data[row[0]] = {'dcdft volume': float(row[1]),
                        'dcdft B0': float(row[2]),
                        'dcdft B1': float(row[3])}

csvwriter2 = csv.writer(open('%s.csv' % tag, 'wb'))
h2 = h + ['%' + h[1], '%' + h[2], '%' + h[3]]
csvwriter2.writerow(h2)

refs = np.loadtxt(reffile,
                  dtype={'names': ('element', 'V0', 'B0', 'BP'),
                         'formats': ('S2', np.float, np.float, np.float)})
# convert into dict
refsd = {}
for e, v, b0, b1 in refs:
    refsd[e] = [v, b0, b1]

rows = []
rowserr = []
for n in task.collection.names:
    row = [n]
    if n in data.keys():
        if 0:
            ref = task.collection.ref[n]  # don't use collection data
        else:
            ref = refsd[n]
        try:
            v = round(data[n]['dcdft volume'], 3)
            b0 = round(data[n]['dcdft B0'], 3)
            b1 = round(data[n]['dcdft B1'], 3)
            row.extend([v, b0, b1])
        except KeyError: # completely failed to find eos minimum
                row.extend(['N/A', 'N/A', 'N/A'])
    else:
        # element not calculated
        row.extend(['N/A', 'N/A', 'N/A'])
    if 'N/A' not in row:
        v0, b00, b10 = ref
        ve = round((v - v0) / v0 * 100, 1)
        b0e = round((b0 - b00) / b00 * 100, 1)
        b1e = round((b1 - b10) / b10 * 100, 1)
        rows.append(row)
        #print row + ref + [ve, b0e, b1e]
        csvwriter2.writerow(row + [ve, b0e, b1e])

# calculate Delta
f = open('%s.txt' % tag, 'wb')
csvwriter3 = csv.writer(f, delimiter='\t')
for r in rows:
    csvwriter3.writerow(r)
f.close()
cmd = 'python ' + os.path.join(dir, 'calcDelta.py')
cmd += ' ' + '%s.txt ' % tag + reffile + ' --stdout'
cmd += ' > ' + '%s_Delta.txt' % tag
os.system(cmd)
