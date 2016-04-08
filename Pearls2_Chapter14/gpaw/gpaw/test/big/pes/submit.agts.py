
def agts(queue):
    calc_co = queue.add('PES_CO.py', ncpus=8, walltime=2 * 45)
    calc_h2o = queue.add('PES_H2O.py', ncpus=8, walltime=2 * 45)
    calc_nh3 = queue.add('PES_NH3.py', ncpus=8, walltime=55)
    queue.add('PES_plot.py', ncpus=1, walltime=5,
              deps=[calc_co, calc_h2o, calc_nh3], creates=['PES_fig.png'])

