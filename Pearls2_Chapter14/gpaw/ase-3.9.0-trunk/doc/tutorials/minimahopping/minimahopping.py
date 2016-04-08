# creates:  summary.png

execfile('Cu2_Pt110.py')
execfile('mhsummary.py')

# Clean up directory.
import os
os.remove('hop.log')
os.remove('minima.traj')
os.remove('qn00000.traj')
os.remove('qn00000.log')
for index in range(1, 10):
    os.remove('qn%05i.traj' % index)
    os.remove('qn%05i.log' % index)
    os.remove('md%05i.traj' % index)
    os.remove('md%05i.log' % index)
