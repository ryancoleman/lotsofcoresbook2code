# creates:  nitrogen.txt
import os
import sys
from StringIO import StringIO

def output_to_string(pythonfile):
    """Returns the stdout of executing the code in pythonfile
    as a string."""
    buffer = StringIO()
    sys.stdout = buffer
    execfile(pythonfile)
    sys.stdout = sys.__stdout__
    return buffer.getvalue()

# Only save the parts relevant to thermochemistry
nitrogen = output_to_string('nitrogen.py')
nitrogen = nitrogen[nitrogen.find('Enthalpy'):]
with open('nitrogen.txt', 'w') as f:
    f.write(nitrogen)
gold = output_to_string('gold.py')
gold = gold[gold.find('Internal'):]
with open('gold.txt', 'w') as f:
    f.write(gold)

# Clean up.
vibfiles = [file for file in os.listdir(os.getcwd()) if
            file.startswith('vib.') or file.startswith('phonon.')]
for file in vibfiles:
    os.remove(file)


