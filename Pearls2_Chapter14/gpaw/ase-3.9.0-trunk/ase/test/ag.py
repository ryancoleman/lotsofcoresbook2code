import os
import sys

from ase import Atoms
from ase.io import write


write('x.json', Atoms('X'))

# Make sure ase-gui can run in terminal mode without $DISPLAY and gtk:
sys.argv = ['ase-gui', '--terminal', 'x.json']
display = os.environ.pop('DISPLAY', None)
error = False
try:
    from ase.gui.ag import main
    main()
    assert 'gtk' not in sys.modules
except:
    error = True

if display is not None:
    os.environ['DISPLAY'] = display

if error:
    raise
