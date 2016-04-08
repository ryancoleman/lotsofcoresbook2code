def agts(queue):
    si = queue.add('submit.agts.py', ncpus=4, walltime=60)
    bs = queue.add('bs.py', deps=si, ncpus=4, walltime=60)
    queue.add('bs_plot.py', deps=bs, creates='bs-PBE.png')

    w = queue.add('wannier.py', deps=si, ncpus=1, walltime=60)
    bs = queue.add('bs_pbe0.py', deps=w, ncpus=1, walltime=10 * 60)
    queue.add('bs_plot2.py', deps=bs, creates='bs-PBE0.png')

    deps = []
    for k in range(3, 8):
        deps.append(queue.add('eos.py %d' % k, ncpus=4, walltime=10 * 60))
    for k in range(8, 13):
        deps.append(queue.add('eos.py %d' % k, ncpus=16, walltime=10 * 60))
    queue.add('results.py', deps=deps, creates=['a.png', 'B.png'])


if __name__ == '__main__':
    from si_pbe import groundstate
    si = groundstate(5.43, 8)
    si.calc.write('Si-PBE.gpw', mode='all')
