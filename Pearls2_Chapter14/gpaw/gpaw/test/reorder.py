from optparse import OptionParser
import numpy as np

usage = '%prog [OPTION] FILE'
description = 'Parse test results in FILE and reorder tests by duration.'

p = OptionParser(usage=usage,
                 description=description)
p.add_option('--pretty', action='store_true',
             help='print nicely with timings.  '
             'Otherwise print Python syntax')
opts, args = p.parse_args()
if len(args) != 1:
    p.error('One file please.')

fd = open(args[0])
for line in fd:
    if line.startswith('===='):
        break

tests = []
durations = []
results = []

for line in fd:
    if line.startswith('===='):
        break
    t1, t2, t3 = line.split()
    tests.append(t1)
    if t3 == 'OK':
        durations.append(float(t2))
    else:
        durations.append(np.inf)
    results.append(t3)

# Use stable sort for minimum shuffling
args = np.argsort(durations, kind='mergesort')
tests = np.array(tests)[args]
durations = np.array(durations)[args]
results = np.array(results)[args]

if opts.pretty:
    for test, duration, result in zip(tests, durations, results):
        print('%30s %10s %10s' % (test, duration, result))
else:
    maxnamelen = max([len(test) for test in tests])
    for test, duration in zip(tests, durations):
        comment = ''
        if duration > 1.0:
            if np.isfinite(duration):
                comment = '# ~%ds' % int(duration)
            else:
                comment = '# duration unknown'
        print(''.join(['    ',
                       ("'%s'," % test).ljust(maxnamelen + 5),
                       comment]).rstrip())
