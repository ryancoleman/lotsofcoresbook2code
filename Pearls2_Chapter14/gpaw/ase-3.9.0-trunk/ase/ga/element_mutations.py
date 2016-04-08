"""Mutation classes, that mutate the elements in the supplied
atoms objects."""
import random
import numpy as np

from ase.data import atomic_numbers
from ase.ga.offspring_creator import OffspringCreator


def chunks(l, n):
    """split a list into smaller chunks"""
    return [l[i:i + n] for i in range(0, len(l), n)]


class ElementMutation(OffspringCreator):
    """The base class for all operators where the elements
    of the atoms objects are mutated"""
    def __init__(self, element_pool, max_diff_elements,
                 min_percentage_elements, verbose):
        OffspringCreator.__init__(self, verbose)
        if not isinstance(element_pool[0], (list, np.ndarray)):
            self.element_pools = [element_pool]
        else:
            self.element_pools = element_pool

        if max_diff_elements is None:
            self.max_diff_elements = [1e6 for _ in self.element_pools]
        elif isinstance(max_diff_elements, int):
            self.max_diff_elements = [max_diff_elements]
        else:
            self.max_diff_elements = max_diff_elements
        assert len(self.max_diff_elements) == len(self.element_pools)

        if min_percentage_elements is None:
            self.min_percentage_elements = [0 for _ in self.element_pools]
        elif isinstance(min_percentage_elements, (int, float)):
            self.min_percentage_elements = [min_percentage_elements]
        else:
            self.min_percentage_elements = min_percentage_elements
        assert len(self.min_percentage_elements) == len(self.element_pools)

    def get_new_individual(self, parents):
        raise NotImplementedError

    def get_mutation_index_list_and_choices(self, atoms):
        """Returns a list of the indices that are going to
        be mutated and a list of possible elements to mutate
        to. The lists obey the criteria set in the initialization.
        """
        itbm_ok = False
        while not itbm_ok:
            itbm = random.choice(range(len(atoms)))  # index to be mutated
            itbm_ok = True
            for i, e in enumerate(self.element_pools):
                if atoms[itbm].symbol in e:
                    elems = e[:]
                    elems_in, indices_in = zip(*[(a.symbol, a.index)
                                                 for a in atoms
                                                 if a.symbol in elems])
                    max_diff_elem = self.max_diff_elements[i]
                    min_percent_elem = self.min_percentage_elements[i]
                    if min_percent_elem == 0:
                        min_percent_elem = 1. / len(elems_in)
                    break
            else:
                itbm_ok = False

        # Check that itbm obeys min/max criteria
        diff_elems_in = len(set(elems_in))
        if diff_elems_in == max_diff_elem:
            # No more different elements allowed -> one element mutation
            ltbm = []  # list to be mutated
            for i in xrange(len(atoms)):
                if atoms[i].symbol == atoms[itbm].symbol:
                    ltbm.append(i)
        else:
            # Fewer or too many different elements already
            if self.verbose:
                print(int(min_percent_elem * len(elems_in)),
                      min_percent_elem, len(elems_in))
            all_chunks = chunks(indices_in,
                                int(min_percent_elem * len(elems_in)))
            itbm_num_of_elems = 0
            for a in atoms:
                if a.index == itbm:
                    break
                if a.symbol in elems:
                    itbm_num_of_elems += 1
            ltbm = all_chunks[itbm_num_of_elems /
                              (int(min_percent_elem * len(elems_in))) - 1]

        elems.remove(atoms[itbm].symbol)

        return ltbm, elems


class RandomElementMutation(ElementMutation):
    """Mutation that exchanges an element with a randomly chosen element from
    the supplied pool of elements
    If the individual consists of different groups of elements the element
    pool can be supplied as a list of lists

    Parameters:

    element_pool: List of elements in the phase space. The elements can be
        grouped if the individual consist of different types of elements.
        The list should then be a list of lists e.g. [[list1], [list2]]

    max_diff_elements: The maximum number of different elements in the
        individual. Default is infinite. If the elements are grouped
        max_diff_elements should be supplied as a list with each input
        corresponding to the elements specified in the same input in
        element_pool.

    min_percentage_elements: The minimum percentage of any element in the
        individual. Default is any number is allowed. If the elements are
        grouped min_percentage_elements should be supplied as a list with
        each input corresponding to the elements specified in the same input
        in element_pool.

    Example: element_pool=[[A,B,C,D],[x,y,z]], max_diff_elements=[3,2],
        min_percentage_elements=[.25, .5]
        An individual could be "D,B,B,C,x,x,x,x,z,z,z,z"
    """
    def __init__(self, element_pool, max_diff_elements=None,
                 min_percentage_elements=None, verbose=False):
        ElementMutation.__init__(self, element_pool, max_diff_elements,
                                 min_percentage_elements, verbose)
        self.descriptor = 'RandomElementMutation'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        ltbm, choices = self.get_mutation_index_list_and_choices(f)

        new_element = random.choice(choices)
        for a in f:
            if a.index in ltbm:
                a.symbol = new_element
            indi.append(a)

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


def mendeleiev_table():
    """Returns the mendeleiev table as a python list of lists.
    Each cell contains either None or a pair (symbol, atomic number),
    or a list of pairs for the cells \* and \**.
    """
    import re
    elems = 'HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKrRb'
    elems += 'SrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm'
    elems += 'YbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMd'
    elems += 'NoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'
    L = [(e, i + 1)
         for (i, e) in enumerate(re.compile('[A-Z][a-z]*').findall(elems))]
    for i, j in ((88, 103), (56, 71)):
        L[i] = L[i:j]
        L[i + 1:] = L[j:]
    for i, j in ((12, 10), (4, 10), (1, 16)):
        L[i:i] = [None] * j
    return [L[18 * i:18 * (i + 1)] for i in range(7)]


def get_row_column(element):
    """Returns the row and column of the element in the periodic table"""
    t = mendeleiev_table()
    for r in range(len(t)):
        en = (element, atomic_numbers[element])
        if en in t[r]:
            return r, t[r].index(en)


class MoveDownMutation(ElementMutation):
    """
    Mutation that exchanges an element with an element one step
    (or more steps if fewer is forbidden) down the same
    column in the periodic table.

    This mutation is introduced and used in:
    P. B. Jensen et al., Phys. Chem. Chem. Phys., 16, 36, 19732-19740 (2014)

    The idea behind is that elements close to each other in the
    periodic table is chemically similar, and therefore exhibit
    similar properties. An individual in the population is
    typically close to fittest possible, exchanging an element
    with a similar element will normally result in a slight
    increase (or decrease) in fitness.

    Parameters:

    element_pool: List of elements in the phase space. The elements can be
        grouped if the individual consist of different types of elements.
        The list should then be a list of lists e.g. [[list1], [list2]]

    max_diff_elements: The maximum number of different elements in the
        individual. Default is infinite. If the elements are grouped
        max_diff_elements should be supplied as a list with each input
        corresponding to the elements specified in the same input in
        element_pool.

    min_percentage_elements: The minimum percentage of any element in the
        individual. Default is any number is allowed. If the elements are
        grouped min_percentage_elements should be supplied as a list with
        each input corresponding to the elements specified in the same input
        in element_pool.

    Example: element_pool=[[A,B,C,D],[x,y,z]], max_diff_elements=[3,2],
        min_percentage_elements=[.25, .5]
        An individual could be "D,B,B,C,x,x,x,x,z,z,z,z"
    """
    def __init__(self, element_pool, max_diff_elements=None,
                 min_percentage_elements=None, verbose=False):
        ElementMutation.__init__(self, element_pool, max_diff_elements,
                                 min_percentage_elements, verbose)
        self.descriptor = 'MoveDownMutation'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        ltbm, choices = self.get_mutation_index_list_and_choices(f)
        ## periodic table row, periodic table column
        ptrow, ptcol = get_row_column(f[ltbm[0]].symbol)

        popped = []
        m = 0
        for j in range(len(choices)):
            e = choices[j - m]
            row, column = get_row_column(e)
            if row <= ptrow or column != ptcol:
                #Throw away if above (lower numbered row)
                #or in a different column in the periodic table
                popped.append(choices.pop(j - m))
                m += 1

        used_descriptor = self.descriptor
        if len(choices) == 0:
            msg = '{0},{2} cannot be mutated by {1}, '
            msg = msg.format(f.info['confid'],
                             self.descriptor,
                             f[ltbm[0]].symbol)
            msg += 'doing random mutation instead'
            if self.verbose:
                print(msg)
            used_descriptor = 'RandomElementMutation_from_{0}'
            used_descriptor = used_descriptor.format(self.descriptor)
            random.shuffle(popped)
            choices = popped
        else:
            # Sorting the element that lie below and in the same column
            # in the periodic table so that the one closest below is first
            choices.sort(key=lambda x: get_row_column(x)[0])
        new_element = choices[0]
        
        for a in f:
            if a.index in ltbm:
                a.symbol = new_element
            indi.append(a)
            
        return (self.finalize_individual(indi),
                used_descriptor + ': {0}'.format(f.info['confid']))
    

class MoveUpMutation(ElementMutation):
    """
    Mutation that exchanges an element with an element one step
    (or more steps if fewer is forbidden) up the same
    column in the periodic table.
    
    This mutation is introduced and used in:
    P. B. Jensen et al., Phys. Chem. Chem. Phys., 16, 36, 19732-19740 (2014)
    
    See MoveDownMutation for the idea behind

    Parameters:

    element_pool: List of elements in the phase space. The elements can be
        grouped if the individual consist of different types of elements.
        The list should then be a list of lists e.g. [[list1], [list2]]
    
    max_diff_elements: The maximum number of different elements in the
        individual. Default is infinite. If the elements are grouped
        max_diff_elements should be supplied as a list with each input
        corresponding to the elements specified in the same input in
        element_pool.

    min_percentage_elements: The minimum percentage of any element in the
        individual. Default is any number is allowed. If the elements are
        grouped min_percentage_elements should be supplied as a list with
        each input corresponding to the elements specified in the same input
        in element_pool.
    
    Example: element_pool=[[A,B,C,D],[x,y,z]], max_diff_elements=[3,2],
        min_percentage_elements=[.25, .5]
        An individual could be "D,B,B,C,x,x,x,x,z,z,z,z"
    """
    def __init__(self, element_pool, max_diff_elements=None,
                 min_percentage_elements=None, verbose=False):
        ElementMutation.__init__(self, element_pool, max_diff_elements,
                                 min_percentage_elements, verbose)
        self.descriptor = 'MoveUpMutation'
        
    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        ltbm, choices = self.get_mutation_index_list_and_choices(f)

        ## periodic table row, periodic table column
        ptrow, ptcol = get_row_column(f[ltbm[0]].symbol)

        popped = []
        m = 0
        for j in range(len(choices)):
            e = choices[j - m]
            row, column = get_row_column(e)
            if row >= ptrow or column != ptcol:
                #Throw away if below (higher numbered row)
                #or in a different column in the periodic table
                popped.append(choices.pop(j - m))
                m += 1

        used_descriptor = self.descriptor
        if len(choices) == 0:
            msg = '{0},{2} cannot be mutated by {1}, '
            msg = msg.format(f.info['confid'],
                             self.descriptor,
                             f[ltbm[0]].symbol)
            msg += 'doing random mutation instead'
            if self.verbose:
                print(msg)
            used_descriptor = 'RandomElementMutation_from_{0}'
            used_descriptor = used_descriptor.format(self.descriptor)
            random.shuffle(popped)
            choices = popped
        else:
            # Sorting the element that lie above and in the same column
            # in the periodic table so that the one closest above is first
            choices.sort(key=lambda x: get_row_column(x)[0], reverse=True)
        new_element = choices[0]

        for a in f:
            if a.index in ltbm:
                a.symbol = new_element
            indi.append(a)
            
        return (self.finalize_individual(indi),
                used_descriptor + ': {0}'.format(f.info['confid']))
        
        
class MoveRightMutation(ElementMutation):
    """
    Mutation that exchanges an element with an element one step
    (or more steps if fewer is forbidden) to the right in the
    same row in the periodic table.
    
    This mutation is introduced and used in:
    P. B. Jensen et al., Phys. Chem. Chem. Phys., 16, 36, 19732-19740 (2014)
    
    See MoveDownMutation for the idea behind

    Parameters:

    element_pool: List of elements in the phase space. The elements can be
        grouped if the individual consist of different types of elements.
        The list should then be a list of lists e.g. [[list1], [list2]]
    
    max_diff_elements: The maximum number of different elements in the
        individual. Default is infinite. If the elements are grouped
        max_diff_elements should be supplied as a list with each input
        corresponding to the elements specified in the same input in
        element_pool.

    min_percentage_elements: The minimum percentage of any element in the
        individual. Default is any number is allowed. If the elements are
        grouped min_percentage_elements should be supplied as a list with
        each input corresponding to the elements specified in the same input
        in element_pool.
    
    Example: element_pool=[[A,B,C,D],[x,y,z]], max_diff_elements=[3,2],
        min_percentage_elements=[.25, .5]
        An individual could be "D,B,B,C,x,x,x,x,z,z,z,z"
    """
    def __init__(self, element_pool, max_diff_elements=None,
                 min_percentage_elements=None, verbose=False):
        ElementMutation.__init__(self, element_pool, max_diff_elements,
                                 min_percentage_elements, verbose)
        self.descriptor = 'MoveRightMutation'
        
    def get_new_individual(self, parents):
        f = parents[0]
        
        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]
        
        ltbm, choices = self.get_mutation_index_list_and_choices(f)
        ## periodic table row, periodic table column
        ptrow, ptcol = get_row_column(f[ltbm[0]].symbol)

        popped = []
        m = 0
        for j in range(len(choices)):
            e = choices[j - m]
            row, column = get_row_column(e)
            if row != ptrow or column <= ptcol:
                #Throw away if to the left (a lower numbered column)
                #or in a different row in the periodic table
                popped.append(choices.pop(j - m))
                m += 1

        used_descriptor = self.descriptor
        if len(choices) == 0:
            msg = '{0},{2} cannot be mutated by {1}, '
            msg = msg.format(f.info['confid'],
                             self.descriptor,
                             f[ltbm[0]].symbol)
            msg += 'doing random mutation instead'
            if self.verbose:
                print(msg)
            used_descriptor = 'RandomElementMutation_from_{0}'
            used_descriptor = used_descriptor.format(self.descriptor)
            random.shuffle(popped)
            choices = popped
        else:
            #Sorting so the element closest to the right is first
            choices.sort(key=lambda x: get_row_column(x)[1])
        new_element = choices[0]
        
        for a in f:
            if a.index in ltbm:
                a.symbol = new_element
            indi.append(a)
            
        return (self.finalize_individual(indi),
                used_descriptor + ': {0}'.format(f.info['confid']))
        
        
class MoveLeftMutation(ElementMutation):
    """
    Mutation that exchanges an element with an element one step
    (or more steps if fewer is forbidden) to the left in the
    same row in the periodic table.
    
    This mutation is introduced and used in:
    P. B. Jensen et al., Phys. Chem. Chem. Phys., 16, 36, 19732-19740 (2014)
    
    See MoveDownMutation for the idea behind

    Parameters:

    element_pool: List of elements in the phase space. The elements can be
        grouped if the individual consist of different types of elements.
        The list should then be a list of lists e.g. [[list1], [list2]]
    
    max_diff_elements: The maximum number of different elements in the
        individual. Default is infinite. If the elements are grouped
        max_diff_elements should be supplied as a list with each input
        corresponding to the elements specified in the same input in
        element_pool.

    min_percentage_elements: The minimum percentage of any element in the
        individual. Default is any number is allowed. If the elements are
        grouped min_percentage_elements should be supplied as a list with
        each input corresponding to the elements specified in the same input
        in element_pool.
    
    Example: element_pool=[[A,B,C,D],[x,y,z]], max_diff_elements=[3,2],
        min_percentage_elements=[.25, .5]
        An individual could be "D,B,B,C,x,x,x,x,z,z,z,z"
    """
    def __init__(self, element_pool, max_diff_elements=None,
                 min_percentage_elements=None, verbose=False):
        ElementMutation.__init__(self, element_pool, max_diff_elements,
                                 min_percentage_elements, verbose)
        self.descriptor = 'MoveLeftMutation'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        ltbm, choices = self.get_mutation_index_list_and_choices(f)
        ## periodic table row, periodic table column
        ptrow, ptcol = get_row_column(f[ltbm[0]].symbol)

        popped = []
        m = 0
        for j in range(len(choices)):
            e = choices[j - m]
            row, column = get_row_column(e)
            if row != ptrow or column >= ptcol:
                #Throw away if to the right (a higher numbered column)
                #or in a different row in the periodic table
                popped.append(choices.pop(j - m))
                m += 1

        used_descriptor = self.descriptor
        if len(choices) == 0:
            msg = '{0},{2} cannot be mutated by {1}, '
            msg = msg.format(f.info['confid'],
                             self.descriptor,
                             f[ltbm[0]].symbol)
            msg += 'doing random mutation instead'
            if self.verbose:
                print(msg)
            used_descriptor = 'RandomElementMutation_from_{0}'
            used_descriptor = used_descriptor.format(self.descriptor)
            random.shuffle(popped)
            choices = popped
        else:
            #Sorting so the element closest to the left is first
            choices.sort(key=lambda x: get_row_column(x)[1], reverse=True)
        new_element = choices[0]
        
        for a in f:
            if a.index in ltbm:
                a.symbol = new_element
            indi.append(a)
            
        return (self.finalize_individual(indi),
                used_descriptor + ': {0}'.format(f.info['confid']))
