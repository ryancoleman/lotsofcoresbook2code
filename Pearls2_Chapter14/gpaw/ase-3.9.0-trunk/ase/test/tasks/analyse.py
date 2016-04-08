import csv

import heapq

import numpy as np

from ase.tasks.io import read_json

ann_fontsize = 'small'
label_fontsize = 12

def most_common_value(values, precision=None):
    import heapq
    # find the most common value
    if precision is not None:
        vs = [round(v, precision) for v in values]
    else:
        vs = values[:]
    if len(vs) > 0:
        nlargest = heapq.nlargest(1, [vs.count(e) for e in vs])[0]
        index = [vs.count(e) for e in vs].index(nlargest)
        v = vs[index]
    else:
        v = None
    return v


def get_differences(reference, values, precision=None):
    # diffrences with respect to the reference value within precision
    values = values[:]
    for n, v in enumerate(values):
        if v is not None:
            de = v - reference
            if precision is not None:
                if abs(de) < 1. / precision:
                    # value within precision
                    values[n] = 0.0
                else:
                    values[n] = round(de, precision)
            else:
                values[n] = de
    return values


def get_key_summary_list(key, name, runs, data, precision=3, relative=False):
    # if relative, then first entry is a common value, other entries relative
    # if not relative, then first entry is average value, other entries verbatim
    nruns = len(runs)
    d = data.copy()
    l = [name]
    values = [d[r][key] for r in runs if r in d and key in d[r] and d[r][key] is not None]
    if relative:
        # find the most common value
        e = most_common_value(values, precision)
    else:
        e = np.mean(values)
    if e is None:
        l.extend(['None' for i in range(nruns + 1)])
    else:
        # print the found common value with the precision
        l.append(("%." + "%df" % (precision + 0)) % e)
        es = [d[r].get(key, None) for r in runs if r in d]
        if relative:
            des = get_differences(e, es)
            # and the corresponding differences
            for n, de in enumerate(des):
                if de is not None:
                    if abs(de) < 1. / precision:
                        # value within precision
                        des[n] = '-'
                    else:
                        des[n] = ("%." + "%df" % precision) % de
                else:
                    des[n] = 'None'
            l.extend(des)
        else:
            des = es[:]
            for n, de in enumerate(es):
                if de is not None:
                    des[n] = str(de)
                else:
                    des[n] = 'None'
            l.extend(des)
    return l


def get_key_stats(data, key, nnlargest):
    keys = []
    values = []
    stats = {'absaverage': None,
             'average': None,
             'std': None,
             'max': None,
             'N': None,
             'nlargest': []}
    for k, v in data.iteritems():
        if key in v:
            if v.get(key) is not None:
                keys.append(k)
                values.append(v[key])
    if len(keys) > 0:
        stats['absaverage'] = np.average(abs(np.array(values)))
        stats['average'] = np.average(np.array(values))
        stats['std'] = np.std(np.array(values))
        stats['max'] = np.max(np.array(values))
        stats['N'] = len(values)
        nlargest = heapq.nlargest(nnlargest, np.abs(np.array(values)))
        # already found elements must be removed in order to avoid repetition
        k = keys[:]
        v = values[:]
        nl = []
        for i in nlargest:
            ind = v.index(i)
            nl.append((k[ind], i))
            v.pop(ind)
            k.pop(ind)
        stats['nlargest'] = nl
    return stats

def plot_single(xdata, ydata, std,
         title,
         xlabel, ylabel,
         label, color, alpha,
         miny, maxy,
         num=1,
         ):
    import matplotlib
    matplotlib.rc('text', usetex=False)
    import pylab
    import matplotlib.font_manager

    # all goes to figure num
    pylab.figure(num=num, figsize=(9.5, 9))
    pylab.gca().set_position([0.10, 0.20, 0.85, 0.60])
    # let the plot have fixed y-axis scale
    ywindow = maxy - miny
    pylab.bar(xdata, ydata, 0.9, yerr=std,
              label=label, color=color, alpha=alpha)
    pylab.gca().set_ylim(miny, maxy)
    t = pylab.title(title)
    # http://old.nabble.com/More-space-between-title-and-secondary-x-axis-td31722298.html
    t.set_y(1.05)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    prop = matplotlib.font_manager.FontProperties(size=12)
    leg = pylab.legend(loc='upper right', fancybox=True, prop=prop)
    leg.get_frame().set_alpha(0.5)

def plot(runs, data,
         labels, # labels coresponding to runs
         key, # the main quantity to plot
         failkey, # the key that informs about failure/success in run
         nnlargest,
         plot_label,
         plot_xlabel,
         plot_ylabel,
         miny, maxy,
         tag=None,
         tunit='min'):
    # horror function
    d = data.copy()
    # collect stats for each run
    stats = {'average': [],
             'std': [],
             'sdom': [],
             'nlargest': [],
             'time': [],
             'converged': []}
    for n, r in enumerate(runs):
        s = get_key_stats(d[r], key, nnlargest)
        assert s['average'] is not None, r + ': no key ' + "'" + key + "'" + ' in results?'
        stats['average'].append(s['average'])
        stats['std'].append(s['std'])
        stats['sdom'].append(s['std']/np.sqrt(s['N']))
        stats['nlargest'].append(s['nlargest'])
        # total run time
        t = [d[r][k]['time'] for k in d[r].keys() if 'time' in d[r][k]]
        if tunit == 'sec':
            stats['time'].append(sum(t)) # sec
        elif tunit == 'min':
            stats['time'].append(sum(t) / 60) # min
        elif tunit == 'h':
            stats['time'].append(sum(t) / 3600) # hours
        # number of converged systems
        c = [d[r][k][failkey] for k in d[r].keys() if failkey in d[r][k]]
        stats['converged'].append(len([s for s in c if s is not None]))
    import matplotlib
    #matplotlib.use('Agg')
    matplotlib.rc('text', usetex=False)
    from matplotlib import pylab, ticker
    num=1
    scale = [i for i in range(len(runs))]
    # plot skeleton
    plot_single(
        scale, stats['average'], stats['sdom'],
        plot_label,
        plot_xlabel,
        plot_ylabel,
        'Average (avg)',
        'blue',
        0.5,
        miny, maxy,
        num=1,
        )
    zero = [0.0 for v in scale]
    pylab.plot(scale, zero, 'k-', label='_nolegend_')
    ay1=pylab.gca()
    ay1.xaxis.set_ticks([n + 0.5 for n in scale])
    ay1.xaxis.set_ticklabels(labels)
    ay1.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    for label in ay1.get_xticklabels() + ay1.get_yticklabels():
        label.set_fontsize(label_fontsize)
    # rotate labels http://old.nabble.com/Rotate-x-axes-%28xticks%29-text.-td3141258.html
    for label in ay1.get_xticklabels():
        label.set_rotation(75)
    # plot stats
    for n, r in enumerate(runs):
        label = ''
        for n1, v in enumerate(['average', 'std', 'time']):
            l = {
                'average': '\navg\n',
                'absaverage': '\nabs\n',
                'std': '\nstd\n',
                'time': '\nt [%s]\n' % tunit,
                }[v]
            value = {
                'average': stats['average'][n],
                #'absaverage': stats['absaverage'][n],
                'std': stats['std'][n],
                'time': stats['time'][n],
                }[v]
            label += l + ' ' + str(round(value, 1)) + '\n'
        label += '\nconv.\n' + ' ' + str(int(stats['converged'][n])) + '\n'
        pylab.annotate(label,
                       xy=(n + 0.0, 0.0),
                       xytext=(n + 0.2, miny / 2 - 0.5),
                       arrowprops=None,
                       horizontalalignment='left',
                       verticalalignment='center',
                       fontsize=ann_fontsize,
                       )
    # plot compounds with largest errors
    for n, l in enumerate(stats['nlargest']):
        l.reverse()
        for n1, (c, e) in enumerate(l):
            name = c + '\n'
            label = name + ' ' + str(int(e))
            xytext=(n + 0.3, maxy / 2 + maxy / 10 * n1)
            if e > maxy:
                # labels exceeding y-scale
                xy = xytext
            else:
                xy=(n + 0.05, e)
            pylab.annotate(label,
                           xy=xy,
                           xytext=xytext,
                           arrowprops=dict(width=0.05,
                                           headwidth=5.0,
                                           facecolor='black',
                                           shrink=1.00),
                           horizontalalignment='left',
                           verticalalignment='center',
                           fontsize=ann_fontsize,
                           )
    # axis limits after the last plot!
    pylab.ylim(miny, maxy)
    kname = key.replace(' ', '_')
    if tag is not None:
        plotname = '%s_%s.png' % (tag, kname)
    else:
        plotname = '%s.png' % kname
    pylab.savefig(plotname, bbox_inches='tight')

class AnalyseOptimizersTask:

    def __init__(self, taskname, runs,
                 labels=None, tag=None, steps=100, precision=3, tunit='min'):
        """Analyse optimizers runs.  """

        self.taskname = taskname
        self.runs = runs
        if labels is None: # run labels for plotting
            self.labels = self.runs.split(',')
        else:
            self.labels = labels
        self.tag = tag
        self.steps = steps # only for purpose of limit plot y-axis
        self.precision = precision
        assert tunit in ['sec', 'min', 'h']
        self.tunit = tunit

        self.key_summary = 'relaxed energy'
        self.key_plot = 'optimizer force calls'

        self.plot_label = 'Summary of optimizers'
        self.plot_xlabel = 'optimizer'
        self.plot_ylabel = self.key_plot

    def read_run(self, taskname, tag, run):
        if tag is not None:
            tag += '_%s' % run
        else:
            tag = run

        return read_json(taskname + '-' + tag + '.json')

    def get_data(self, taskname, tag, runs):
        d = {}
        for r in runs:
            d[r] = self.read_run(taskname, tag, r)
        return d

    def analyse(self):
        runs = self.runs.split(',')
        # dict of systems for each run
        datarun = self.get_data(self.taskname, self.tag, runs)
        # the number of systems (based on results from json)
        nsystems = max([len(datarun[n]) for n in datarun.keys()])
        # dict of runs for each system
        datasys = {}
        for k, v in datarun.iteritems():
            for k1, v1 in v.iteritems():
                if k1 not in datasys:
                    datasys[k1] = {k: v1}
                else:
                    datasys[k1].update({k: v1})
        # csv summary of self.key_plot
        key_name = self.key_plot.replace(' ', '_')
        row = ['formula', self.key_plot]
        row.extend([r for r in range(len(runs))])
        rows = [row]
        for name, data in datasys.items():
            if not data:
                continue
            row = get_key_summary_list(self.key_plot,
                                       name,
                                       runs,
                                       data,
                                       precision=self.precision,
                                       relative=False)
            row = [r.replace('None', 'N/A') for r in row]
            # only failed or non-common runs
            for k in row[2:]:
                if k == 'N/A' or k != '-':
                    rows.append(row)
                    break
        if len(rows) > 0: # always create csv file
            if self.tag is not None:
                csvwriter = csv.writer(
                    open('%s_%s.csv' % (self.tag, key_name), 'wb'))
            else:
                csvwriter = csv.writer(open('%s.csv' % key_name, 'wb'))
            for r in rows:
                csvwriter.writerow(r)
        # csv summary of self.key_summary
        summary_name = self.key_summary.replace(' ', '_')
        row = ['formula', self.key_summary]
        row.extend([r for r in range(len(runs))])
        rows = [row]
        for name, data in datasys.items():
            if not data:
                continue
            row = get_key_summary_list(self.key_summary,
                                       name,
                                       runs,
                                       data,
                                       precision=self.precision,
                                       relative=True)
            row = [r.replace('None', 'N/A') for r in row]
            # only failed or non-common runs
            for k in row[2:]:
                if k == 'N/A' or k != '-':
                    rows.append(row)
                    break
        if len(rows) > 0: # always create csv file
            if self.tag is not None:
                csvwriter = csv.writer(
                    open('%s_%s.csv' % (self.tag, summary_name), 'wb'))
            else:
                csvwriter = csv.writer(open('%s.csv' % summary_name, 'wb'))
            for r in rows:
                csvwriter.writerow(r)
        # plot
        maxy = self.steps - 5
        miny = - ((maxy + 5) / 2 + 5)
        # plot only runs with data
        plotruns = runs[:]
        labels = self.labels[:]
        for n, r in enumerate(runs):
            nd = [1 for k in datarun[r] if datarun[r][k]]
            if len(nd) == 0:
                print 'skipped plotting empty data', r
                # no data for this run
                i = plotruns.index(r)
                plotruns.pop(i)
                labels.pop(i)
            continue
        plot(plotruns, datarun,
             labels,
             self.key_plot, # the main quantity to plot
             self.key_summary,
             4,
             self.plot_label + ': ' + str(nsystems) + ' systems',
             self.plot_xlabel,
             self.plot_ylabel,
             miny, maxy,
             tag=self.tag,
             tunit=self.tunit)

class AnalyseSCFTask(AnalyseOptimizersTask):

    def __init__(self, taskname, runs,
                 labels=None, tag=None, steps=100, precision=3, tunit='min'):
        """Analyse SCF runs.  """

        AnalyseOptimizersTask.__init__(self, taskname, runs,
                                       labels, tag, steps, precision, tunit)

        self.key_summary = 'energy'
        self.key_plot = 'calculator steps'

        self.plot_label = 'Summary of runs'
        self.plot_xlabel = 'run'
        self.plot_ylabel = self.key_plot
