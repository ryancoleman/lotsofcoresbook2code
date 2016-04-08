import os
import shutil

machine = os.environ.get('MACHINE', 'TEST')
shutil.rmtree(machine, ignore_errors=True)
os.system('sh prepare.sh')
