# creates: citations.png citations.csv
from __future__ import print_function
import os
import datetime

import matplotlib.pyplot as plt


months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def f(filename):
    papers = {}
    lines = open(filename).readlines()
    n = 0
    dois = set()
    while n < len(lines):
        line = lines[n]
        tag = line[:2]
        if tag == 'TI':
            ntitle = n
            y = None
            m = 1
            d = 15
        elif tag == 'SO':
            title = ' '.join(lines[i][3:-1] for i in range(ntitle, n))
        elif tag == 'DI':
            doi = line[3:-1]
        elif tag == 'PY':
            y = int(line.split()[1])
        elif tag == 'PD':
            for w in line.split()[1:]:
                if w[0].isdigit():
                    w = int(w)
                    if w < 100:
                        d = w
                    else:
                        y = w
                else:
                    if '-' in w:
                        w = w.split('-')[-1]
                    m = months.index(w) + 1
        elif tag == '\n':
            date = datetime.date(y, m, d)
            if doi not in dois:
                dois.add(doi)
                papers[doi] = (date, title)
        n += 1

    return papers


plt.figure(figsize=(8, 4))
total = {}
for bib in ['gpaw1', 'tddft', 'lcao', 'gpaw2', 'response']:
    papers = {}
    for line in open(bib + '.txt'):
        date, doi, title = line.split(' ', 2)
        papers[doi] = (datetime.date(*[int(x) for x in date.split('-')]),
                       title.strip())
    if os.path.isfile(bib + '.bib'):
        papers.update(f(bib + '.bib'))
    papers = [(papers[doi][0], doi, papers[doi][1]) for doi in papers]
    papers.sort()
    plt.plot([paper[0] for paper in papers], range(1, len(papers) + 1),
             '-o', label=bib)
    fd = open(bib + '.txt', 'w')
    for date, doi, title in papers:
        fd.write('%d-%02d-%02d %s %s\n' % (date.year, date.month, date.day,
                                           doi, title))
        total[doi] = (date, title)
    fd.close()
    x = dict([(p[1], 0) for p in papers])
    print((bib, len(papers), len(x), len(total)))


allpapers = [(paper[0], doi, paper[1]) for doi, paper in total.items()]
allpapers.sort()
plt.plot([paper[0] for paper in allpapers], range(1, len(allpapers) + 1),
         '-o', label='total')

fd = open('citations.csv', 'w')
n = len(allpapers)
for date, doi, title in allpapers[::-1]:
    fd.write('%d,"`%s <http://dx.doi.org/%s>`__"\n' % (n, title, doi))
    n -= 1
fd.close()

plt.xlabel('date')
plt.ylabel('number of citations')
plt.legend(loc='upper left')
plt.savefig('citations.png')
