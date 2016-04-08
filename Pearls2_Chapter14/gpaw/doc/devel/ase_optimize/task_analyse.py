import os

import sys

import glob

from ase.test.tasks.analyse import AnalyseOptimizersTask as Task

def get_optimizer_shortname(optimizer):
    name =  {
        'FIREE19': 'FIRE\neigenstates 1.e-9',
        'GoodOldQuasiNewtonD16': 'GoodOldQuasiNewton\ndensity 1.e-6',
        'BFGSLineSearch': 'BFGSLineSearch',
        'LBFGSLineSearch': 'LBFGSLineSearch',
        'SciPyFminCG': 'SciPyFminCG',
        'SciPyFminBFGS': 'SciPyFminBFGS',
        'GoodOldQuasiNewton': 'GoodOldQuasiNewton'}.get(optimizer, None)
    if name is None:
        return optimizer
    else:
        return name

if __name__ == '__main__':

    assert len(sys.argv) > 1
    if len(sys.argv) == 2:
        taskname = sys.argv[1]
        tag = None
        runs = None
    if len(sys.argv) == 3:
        taskname = sys.argv[1]
        tag = sys.argv[2]
        runs = None
    if len(sys.argv) == 4:
        taskname = sys.argv[1]
        tag = sys.argv[2]
        runs = sys.argv[3]

    if runs is None:  # use all json files as runs
        runs = []
        for f in glob.glob(taskname + '-' + tag + '*.json'):
            runs.append(os.path.splitext(f)[0].split('_')[-1])
    else:
        runs = runs.split(',')
    labels = [get_optimizer_shortname(r) for r in runs]

    steps = 50
    t = Task(taskname, ','.join(runs), labels=labels, tag=tag, steps=steps,
             tunit='h')
    t.analyse()
