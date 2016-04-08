""" Implementaiton of a population for maintaining a GA population and
proposing structures to pair. """
from random import randrange, random
from math import tanh, sqrt

from ase.db.core import now


def count_looks_like(a, all_cand, comp):
    """Utility method for counting occurences."""
    n = 0
    for b in all_cand:
        if a.info['confid'] == b.info['confid']:
            continue
        if comp.looks_like(a, b):
            n += 1
    return n


class Population(object):
    """Population class which maintains the current population
    and proposes which candidates to pair together.

    Parameters:

    data_connection: DataConnection object
        Bla bla bla.
    population_size: int
        The number of candidates in the population.
    comparator: Comparator object
        this will tell if two configurations are equal.
        Default compare atoms objects directly.
    logfile: str
        Text file that contains information about the population
        The format is::

            timestamp: generation(if available): id1,id2,id3...

        Using this file greatly speeds up convergence checks.
        Default None meaning that no file is written.
    """
    def __init__(self, data_connection, population_size,
                 comparator=None, logfile=None):
        self.dc = data_connection
        self.pop_size = population_size
        if comparator is None:
            from ase.ga.standard_comparators import AtomsComparator
            comparator = AtomsComparator()
        self.comparator = comparator
        self.logfile = logfile
        self.pop = []
        self.pairs = None
        self.all_cand = None
        self.__initialize_pop__()

    def __initialize_pop__(self):
        """ Private method that initalizes the population when
            the population is created. """

        # Get all relaxed candidates from the database
        all_cand = self.dc.get_all_relaxed_candidates()
        all_cand.sort(key=lambda x: x.get_raw_score(), reverse=True)
        # all_cand.sort(key=lambda x: x.get_potential_energy())

        # Fill up the population with the self.pop_size most stable
        # unique candidates.
        i = 0
        while i < len(all_cand) and len(self.pop) < self.pop_size:
            c = all_cand[i]
            i += 1
            eq = False
            for a in self.pop:
                if self.comparator.looks_like(a, c):
                    eq = True
                    break
            if not eq:
                self.pop.append(c)

        for a in self.pop:
            a.info['looks_like'] = count_looks_like(a, all_cand,
                                                    self.comparator)

        self.all_cand = all_cand
        self.__calc_participation__()

    def __calc_participation__(self):
        """ Determines, from the database, how many times each
            candidate has been used to generate new candidates. """
        (participation, pairs) = self.dc.get_participation_in_pairing()
        for a in self.pop:
            if a.info['confid'] in participation.keys():
                a.info['n_paired'] = participation[a.info['confid']]
            else:
                a.info['n_paired'] = 0
        self.pairs = pairs

    def update(self):
        """ New candidates can be added to the database
            after the population object has been created.
            This method extracts these new candidates from the
            database and includes them in the population. """

        if len(self.pop) == 0:
            self.__initialize_pop__()

        new_cand = self.dc.get_all_relaxed_candidates(only_new=True)
        for a in new_cand:
            self.__add_candidate__(a)
            self.all_cand.append(a)
        self.__calc_participation__()
        self._write_log()

    def get_current_population(self):
        """ Returns a copy of the current population. """
        self.update()
        return [a.copy() for a in self.pop]

    def get_population_after_generation(self, gen):
        """ Returns a copy of the population as it where
        after generation gen"""
        if self.logfile is not None:
            f = open(self.logfile, 'r')
            gens = {}
            for l in f:
                _, no, popul = l.split(':')
                gens[int(no)] = popul.split(',')
            f.close()
            return [c.copy() for c in self.all_cand[::-1]
                    if c.info['relax_id'] in gens[gen]]

        all_candidates = [c for c in self.all_cand
                          if c.info['key_value_pairs']['generation'] <= gen]
        cands = [all_candidates[0]]
        for b in all_candidates:
            if b not in cands:
                for a in cands:
                    if self.comparator.looks_like(a, b):
                        break
                else:
                    cands.append(b)
        pop = cands[:self.pop_size]
        return [a.copy() for a in pop]

    def __add_candidate__(self, a):
        """ Adds a single candidate to the population. """

        #check if the structure is too low in raw score
        if a.get_raw_score() < self.pop[-1].get_raw_score() \
                and len(self.pop) == self.pop_size:
            return

        # check if the new candidate should
        # replace a similar structure in the population
        for (i, b) in enumerate(self.pop):
            if self.comparator.looks_like(a, b):
                if b.get_raw_score() < a.get_raw_score():
                    del self.pop[i]
                    a.info['looks_like'] = count_looks_like(a,
                                                            self.all_cand,
                                                            self.comparator)
                    self.pop.append(a)
                    self.pop.sort(key=lambda x: x.get_raw_score(),
                                  reverse=True)
                return

        # the new candidate needs to be added, so remove the highest
        # energy one
        if len(self.pop) == self.pop_size:
            del self.pop[-1]

        # add the new candidate
        a.info['looks_like'] = count_looks_like(a,
                                                self.all_cand,
                                                self.comparator)
        self.pop.append(a)
        self.pop.sort(key=lambda x: x.get_raw_score(), reverse=True)

    def __get_fitness__(self, indecies, with_history=True):
        """Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816

        Sign change on the fitness compared to the formulation in the
        abovementioned paper due to maximizing raw_score instead of
        minimizing energy. (Set raw_score=-energy to optimize the energy)
        """

        scores = [x.get_raw_score() for x in self.pop]
        min_s = min(scores)
        max_s = max(scores)
        T = min_s - max_s
        if isinstance(indecies, int):
            indecies = [indecies]

        f = [0.5 * (1. - tanh(2. * (scores[i] - max_s) / T - 1.))
             for i in indecies]
        if with_history:
            M = [float(self.pop[i].info['n_paired']) for i in indecies]
            L = [float(self.pop[i].info['looks_like']) for i in indecies]
            f = [f[i] * 1. / sqrt(1. + M[i]) * 1. / sqrt(1. + L[i])
                 for i in xrange(len(f))]
        return f

    def get_two_candidates(self, with_history=True):
        """ Returns two candidates for pairing employing the
            fitness criteria from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            and the roulete wheel selection scheme described in
            R.L. Johnston Dalton Transactions,
            Vol. 22, No. 22. (2003), pp. 4193-4207
        """

        if len(self.pop) < 2:
            self.update()

        if len(self.pop) < 2:
            return None

        fit = self.__get_fitness__(range(len(self.pop)), with_history)
        fmax = max(fit)
        c1 = self.pop[0]
        c2 = self.pop[0]
        used_before = False
        while c1.info['confid'] == c2.info['confid'] and not used_before:
            nnf = True
            while nnf:
                t = randrange(0, len(self.pop), 1)
                if fit[t] > random() * fmax:
                    c1 = self.pop[t]
                    nnf = False
            nnf = True
            while nnf:
                t = randrange(0, len(self.pop), 1)
                if fit[t] > random() * fmax:
                    c2 = self.pop[t]
                    nnf = False

            c1id = c1.info['confid']
            c2id = c2.info['confid']
            used_before = (min([c1id, c2id]), max([c1id, c2id])) in self.pairs
        return (c1.copy(), c2.copy())

    def _write_log(self):
        """Writes the population to a logfile.

        The format is::

            timestamp: generation(if available): id1,id2,id3..."""
        if self.logfile is not None:
            ids = [str(a.info['relax_id']) for a in self.pop]
            if ids != []:
                try:
                    gen_nums = [c.info['key_value_pairs']['generation']
                                for c in self.all_cand]
                    max_gen = max(gen_nums)
                except KeyError:
                    max_gen = ' '
                f = open(self.logfile, 'a')
                f.write('{time}: {gen}: {pop}\n'.format(time=now(),
                                                        pop=','.join(ids),
                                                        gen=max_gen))
                f.close()
