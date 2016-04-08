# creates: h2.txt
import os
from sys import executable
assert os.system('%s h2.py' % executable) == 0
os.system('cp h2.txt ../_build')
