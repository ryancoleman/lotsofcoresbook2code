import random
import numpy as np
from operator import itemgetter

from ase.ga.offspring_creator import OffspringCreator
from ase.ga.utilities import get_distance_matrix, get_nndist
from ase import Atoms


class Mutation(OffspringCreator):
    """Base class for all particle mutation type operators.
    Do not call this class directly."""
    def __init__(self, num_muts):
        OffspringCreator.__init__(self)
        self.descriptor = 'Mutation'
        self.num_muts = num_muts

    @classmethod
    def interchange2(cls, atoms, i1, i2):
        """Switches identity of the atoms on index i1 and i2 in
        the supplied atoms object."""
        p1 = atoms[int(i1)].position.copy()
        atoms[int(i1)].position = atoms[int(i2)].position.copy()
        atoms[int(i2)].position = p1

    @classmethod
    def get_atomic_configuration(cls, atoms, elements=None, eps=4e-2):
        """Returns the atomic configuration of the particle as a list of
        lists. Each list contain the indices of the atoms sitting at the
        same distance from the geometrical center of the particle. Highly
        symmetrical particles will often have many atoms in each shell.

        For further elaboration see:
        J. Montejano-Carrizales and J. Moran-Lopez, Geometrical
        characteristics of compact nanoclusters, Nanostruct. Mater., 1,
        5, 397-409 (1992)

        Parameters:

        elements: Only take into account the elements specified in this
            list. Default is to take all elements into account.

        eps: The distance allowed to separate elements within each shell."""
        atoms = atoms.copy()
        if elements is None:
            e = list(set(atoms.get_chemical_symbols()))
        else:
            e = elements
        atoms.set_constraint()
        atoms.center()
        geo_mid = np.array([(atoms.get_cell() / 2.)[i][i] for i in range(3)])
        dists = [(np.linalg.norm(geo_mid - atoms[i].position), i)
                 for i in range(len(atoms))]
        dists.sort(key=itemgetter(0))
        atomic_conf = []
        old_dist = -10.
        for dist, i in dists:
            if abs(dist - old_dist) > eps:
                atomic_conf.append([i])
            else:
                atomic_conf[-1].append(i)
            old_dist = dist
        sorted_elems = sorted(set(atoms.get_chemical_symbols()))
        if e is not None and sorted(e) != sorted_elems:
            for shell in atomic_conf:
                torem = []
                for i in shell:
                    if atoms[i].symbol not in e:
                        torem.append(i)
                for i in torem:
                    shell.remove(i)
        return atomic_conf

    @classmethod
    def get_list_of_possible_permutations(cls, atoms, l1, l2):
        """Returns a list of available permutations from the two
        lists of indices, l1 and l2. Checking that identical elements
        are not permuted."""
        possible_permutations = []
        for i in l1:
            for j in l2:
                if atoms[int(i)].symbol != atoms[int(j)].symbol:
                    possible_permutations.append((i, j))
        return possible_permutations


class RandomMutation(Mutation):
    """Moves a random atom the supplied length in a random direction."""
    def __init__(self, length=2., num_muts=1):
        Mutation.__init__(self, num_muts)
        self.descriptor = 'RandomMutation'
        self.length = length

    def mutate(self, atoms):
        """ Does the actual mutation. """
        tbm = random.choice(range(len(atoms)))

        indi = Atoms()
        for a in atoms:
            if a.index == tbm:
                a.position += self.random_vector(self.length)
            indi.append(a)
        return indi

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.mutate(f)
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))

    @classmethod
    def random_vector(cls, l):
        """return random vector of length l"""
        vec = np.array([random.random() * 2 - 1 for i in range(3)])
        vl = np.linalg.norm(vec)
        return np.array([v * l / vl for v in vec])


class RandomPermutation(Mutation):
    """Permutes two random atoms.

    Parameters:

    num_muts: the number of times to perform this operation."""
    def __init__(self, num_muts=1):
        Mutation.__init__(self, num_muts)
        self.descriptor = 'RandomPermutation'

    def get_new_individual(self, parents):
        f = parents[0].copy()

        diffatoms = len(set(f.numbers))
        assert diffatoms > 1, 'Permutations with one atomic type is not valid'

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]
        
        for _ in xrange(self.num_muts):
            i1 = random.randrange(len(f))
            i2 = random.randrange(len(f))
            while f[i1].symbol == f[i2].symbol:
                i2 = random.randrange(len(f))
            self.interchange2(f, i1, i2)
            
        for atom in f:
            indi.append(atom)
            
        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


class COM2surfPermutation(Mutation):
    """The Center Of Mass to surface (COM2surf) permutation operator
    described in
    S. Lysgaard et al., Top. Catal., 2014, 57 (1-4), pp 33-39

    Parameters:
    
    elements: which elements should be included in this permutation,
        for example: include all metals and exclude all adsorbates
    
    min_ratio: minimum ratio of each element in the core or surface region.
        If less than minimum ratio is present in the core, the region defining
        the core will be extended untill the minimum ratio is met, and vice
        versa for the surface region. It has the potential reach the
        recursive limit if an element has a smaller total ratio in the
        complete particle. In that case remember to decrease this min_ratio.

    num_muts: the number of times to perform this operation.
    """
    def __init__(self, elements=None, min_ratio=0.25, num_muts=1):
        Mutation.__init__(self, num_muts)
        self.descriptor = 'COM2surfPermutation'
        self.min_ratio = min_ratio
        self.elements = elements
        
    def get_new_individual(self, parents):
        f = parents[0].copy()
    
        diffatoms = len(set(f.numbers))
        assert diffatoms > 1, 'Permutations with one atomic type is not valid'

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        for _ in xrange(self.num_muts):
            elems = self.elements
            atomic_conf = self.get_atomic_configuration(f, elements=elems)
            core = self.get_core_indices(f, atomic_conf, self.min_ratio)
            shell = self.get_shell_indices(f, atomic_conf, self.min_ratio)
            permuts = self.get_list_of_possible_permutations(f, core, shell)
            swap = random.choice(permuts)
            self.interchange2(f, *swap)

        for atom in f:
            indi.append(atom)
            
        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))
        
    def get_core_indices(self, atoms, atomic_conf, min_ratio, recurs=0):
        """Recursive function that returns the indices in the core subject to
        the min_ratio constraint. The indices are found from the supplied
        atomic configuration."""
        elements = list(set([atoms[i].symbol
                             for subl in atomic_conf for i in subl]))
        
        core = [i for subl in atomic_conf[:1 + recurs] for i in subl]
        while len(core) < 1:
            recurs += 1
            core = [i for subl in atomic_conf[:1 + recurs] for i in subl]
            
        for elem in elements:
            ratio = len([i for i in core
                         if atoms[i].symbol == elem]) / float(len(core))
            if ratio < min_ratio:
                return self.get_core_indices(atoms, atomic_conf,
                                             min_ratio, recurs + 1)
        return core
                
    def get_shell_indices(self, atoms, atomic_conf, min_ratio, recurs=0):
        """Recursive function that returns the indices in the surface
        subject to the min_ratio constraint. The indices are found from
        the supplied atomic configuration."""
        elements = list(set([atoms[i].symbol
                             for subl in atomic_conf for i in subl]))
        
        shell = [i for subl in atomic_conf[-1 - recurs:] for i in subl]
        while len(shell) < 1:
            recurs += 1
            shell = [i for subl in atomic_conf[-1 - recurs:] for i in subl]
            
        for elem in elements:
            ratio = len([i for i in shell
                         if atoms[i].symbol == elem]) / float(len(shell))
            if ratio < min_ratio:
                return self.get_shell_indices(atoms, atomic_conf,
                                              min_ratio, recurs + 1)
        return shell


class _NeighborhoodPermutation(Mutation):
    """Helper class that holds common functions to all permutations
    that look at the neighborhoods of each atoms."""
    def get_possible_poor2rich_permutations(self, atoms, inverse=False,
                                            recurs=0, distance_matrix=None):
        dm = distance_matrix
        if dm is None:
            dm = get_distance_matrix(atoms)
        # Adding a small value (0.2) to overcome slight variations
        # in the average bond length
        nndist = get_nndist(atoms, dm) + 0.2
        same_neighbors = {}
        
        def f(x):
            return x[1]
        for i, atom in enumerate(atoms):
            same_neighbors[i] = 0
            neighbors = [j for j in xrange(len(dm[i])) if dm[i][j] < nndist]
            for n in neighbors:
                if atoms[n].symbol == atom.symbol:
                    same_neighbors[i] += 1
        sorted_same = sorted(same_neighbors.iteritems(), key=f)
        if inverse:
            sorted_same.reverse()
        poor_indices = [j[0] for j in sorted_same
                        if abs(j[1] - sorted_same[0][1]) <= recurs]
        rich_indices = [j[0] for j in sorted_same
                        if abs(j[1] - sorted_same[-1][1]) <= recurs]
        permuts = self.get_list_of_possible_permutations(atoms, poor_indices,
                                                         rich_indices)
        
        if len(permuts) == 0:
            return self.get_possible_poor2rich_permutations(atoms, inverse,
                                                            recurs + 1, dm)
        return permuts
            
        
class Poor2richPermutation(_NeighborhoodPermutation):
    """The poor to rich (Poor2rich) permutation operator described in
    S. Lysgaard et al., Top. Catal., 2014, 57 (1-4), pp 33-39
    
    Permutes two atoms from regions short of the same elements, to
    regions rich in the same elements.
    (Inverse of Rich2poorPermutation)

    Parameters:

    elements: Which elements to take into account in this permutation
    """
    def __init__(self, elements=[], num_muts=1):
        _NeighborhoodPermutation.__init__(self, num_muts)
        self.descriptor = 'Poor2richPermutation'
        self.elements = elements

    def get_new_individual(self, parents):
        f = parents[0].copy()
    
        diffatoms = len(set(f.numbers))
        assert diffatoms > 1, 'Permutations with one atomic type is not valid'
        
        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]
        
        for _ in xrange(self.num_muts):
            atoms = f.copy()
            del atoms[[atom.index for atom in atoms
                       if atom.symbol not in self.elements]]
            permuts = self.get_possible_poor2rich_permutations(atoms)
            swap = random.choice(permuts)
            self.interchange2(f, *swap)
            
        for atom in f:
            indi.append(atom)
            
        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))

        
class Rich2poorPermutation(_NeighborhoodPermutation):
    """
    The rich to poor (Rich2poor) permutation operator described in
    S. Lysgaard et al., Top. Catal., 2014, 57 (1-4), pp 33-39

    Permutes two atoms from regions rich in the same elements, to
    regions short of the same elements.
    (Inverse of Poor2richPermutation)

    Parameters:

    elements: Which elements to take into account in this permutation
    """
    def __init__(self, elements=None, num_muts=1):
        _NeighborhoodPermutation.__init__(self, num_muts)
        self.descriptor = 'Rich2poorPermutation'
        self.elements = elements
        self.num_muts = num_muts

    def get_new_individual(self, parents):
        f = parents[0].copy()

        diffatoms = len(set(f.numbers))
        assert diffatoms > 1, 'Permutations with one atomic type is not valid'

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        if self.elements is None:
            elems = list(set(f.get_chemical_symbols()))
        else:
            elems = self.elements
        for _ in xrange(self.num_muts):
            atoms = f.copy()
            del atoms[[atom.index for atom in atoms
                       if atom.symbol not in elems]]
            permuts = self.get_possible_poor2rich_permutations(atoms,
                                                               inverse=True)
            swap = random.choice(permuts)
            self.interchange2(f, *swap)

        for atom in f:
            indi.append(atom)

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


class SymmetricSubstitute(Mutation):
    """Permute all atoms within a subshell of the symmetric particle.
    The atoms within a subshell all have the same distance to the center,
    these are all equivalent under the particle point group symmetry.

    """
    def __init__(self, elements=None, num_muts=1):
        Mutation.__init__(self, num_muts)
        self.descriptor = 'SymmetricSubstitute'
        self.elements = elements

    def substitute(self, atoms):
        """Does the actual substitution"""
        atoms = atoms.copy()
        aconf = self.get_atomic_configuration(atoms,
                                              elements=self.elements)
        itbm = random.randint(0, len(aconf) - 1)
        to_element = random.choice(self.elements)

        for i in aconf[itbm]:
            atoms[i].symbol = to_element

        return atoms

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.substitute(f)
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


class RandomSubstitute(Mutation):
    """Substitutes one atom with another atom type. The possible atom types
    are supplied in the parameter elements"""
    def __init__(self, elements=None, num_muts=1):
        Mutation.__init__(self, num_muts)
        self.descriptor = 'RandomSubstitute'
        self.elements = elements

    def substitute(self, atoms):
        """Does the actual substitution"""
        atoms = atoms.copy()
        if self.elements is None:
            elems = list(set(atoms.get_chemical_symbols()))
        else:
            elems = self.elements[:]
        possible_indices = [a.index for a in atoms
                            if a.symbol in elems]
        itbm = random.choice(possible_indices)
        elems.remove(atoms[itbm].symbol)
        new_symbol = random.choice(elems)
        atoms[itbm].symbol = new_symbol

        return atoms

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.substitute(f)
        indi = self.initialize_individual(f, indi)
        indi.info['data']['parents'] = [f.info['confid']]

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))
