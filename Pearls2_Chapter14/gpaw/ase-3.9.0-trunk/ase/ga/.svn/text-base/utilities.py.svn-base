""" Various utility methods used troughout the GA. """
from ase.data import covalent_radii
import itertools
import numpy as np
from ase.io import write, read
import os
import time
import math


def closest_distances_generator(atom_numbers, ratio_of_covalent_radii):
    """ Generates the blmin dict used across the GA.
        The distances are based on the covalent radii of the atoms.
    """
    cr = covalent_radii
    ratio = ratio_of_covalent_radii

    blmin = dict()
    for i in atom_numbers:
        blmin[(i, i)] = cr[i] * 2 * ratio
        for j in atom_numbers:
            if i == j:
                continue
            if (i, j) in blmin.keys():
                continue
            blmin[(i, j)] = blmin[(j, i)] = ratio * (cr[i] + cr[j])
    return blmin


def get_mic_distance(p1, p2, cell, pbc):
    """ This method calculates the shortest distance between p1 and p2
         through the cell boundaries defined by cell and pbc.
         This method works for reasonable unit cells, but not for extremely
         elongated ones.
    """
    ct = cell.T
    pos = np.mat((p1, p2))
    scaled = np.linalg.solve(ct, pos.T).T
    for i in xrange(3):
        if pbc[i]:
            scaled[:, i] %= 1.0
            scaled[:, i] %= 1.0
    P = np.dot(scaled, cell)

    pbc_directions = [[-1, 1] * int(direction) + [0] for direction in pbc]
    translations = np.mat(list(itertools.product(*pbc_directions))).T
    p0r = np.tile(np.reshape(P[0, :], (3, 1)), (1, translations.shape[1]))
    p1r = np.tile(np.reshape(P[1, :], (3, 1)), (1, translations.shape[1]))
    dp_vec = p0r + ct * translations
    d = np.min(np.power(p1r - dp_vec, 2).sum(axis=0))**0.5
    return d


def db_call_with_error_tol(db_cursor, expression, args=[]):
    """ In case the GA is used on older versions of networking
         filesystems there might be some delays. For this reason
         some extra error tolerance when calling the SQLite db is
         employed.
    """
    import sqlite3
    i = 0
    while i < 10:
        try:
            db_cursor.execute(expression, args)
            return
        except sqlite3.OperationalError, e:
            print(e)
            time.sleep(2.)
        i += 1
    raise sqlite3.OperationalError(
        'Database still locked after 10 attempts (20 s)')


def save_trajectory(confid, trajectory, folder):
    """ Saves traj files to the database folder.
         This method should never be used directly,
         but only through the DataConnection object.
    """
    fname = os.path.join(folder, 'traj%05d.traj' % confid)
    write(fname, trajectory)
    return fname


def get_trajectory(fname):
    """ Extra error tolerance when loading traj files. """
    fname = str(fname)
    try:
        t = read(fname)
    except IOError, e:
        print('get_trajectory error ' + e)
    return t


def atoms_too_close(a, bl):
    """ Checks if any atoms in a are too close, as defined by
        the distances in the bl dictionary. """
    num = a.numbers
    for i in xrange(len(a)):
        for j in xrange(i + 1, len(a)):
            if a.get_distance(i, j, True) < bl[(num[i], num[j])]:
                return True
    return False


def atoms_too_close_two_sets(a, b, bl):
    """ Checks if any atoms in a are too close to an atom in b,
        as defined by the bl dictionary. """
    tot = a + b
    num = tot.numbers
    for i in xrange(len(a)):
        for j in xrange(len(a), len(tot)):
            if tot.get_distance(i, j, True) < bl[(num[i], num[j])]:
                return True
    return False


def get_all_atom_types(slab, atom_numbers_to_optimize):
    """ Utility method used to extract all unique atom types
        from the atoms object slab and the list of atomic numbers
        atom_numbers_to_optimize. """
    from_slab = list(set(slab.numbers))
    from_top = list(set(atom_numbers_to_optimize))
    from_slab.extend(from_top)
    return list(set(from_slab))


def get_distance_matrix(atoms, self_distance=1000):
    """ Returns a numpy matrix with the distances between the atoms
        in the supplied atoms object, with the indices of the matrix
        corresponding to the indices in the atoms object.
        The parameter self_distance will be put in the diagonal
        elements ([i][i])
    """
    dm = np.zeros([len(atoms), len(atoms)])
    for i in xrange(len(atoms)):
        dm[i][i] = self_distance
        for j in xrange(i + 1, len(atoms)):
            rij = atoms.get_distance(i, j)
            dm[i][j] = rij
            dm[j][i] = rij
    return dm


def get_rdf(atoms, rmax, nbins, distance_matrix=None):
    """
    Returns two numpy arrays; the radial distribution function
    and the corresponding distances of the supplied atoms object
    """
    dm = distance_matrix
    if dm is None:
        dm = get_distance_matrix(atoms)
    rdf = np.zeros(nbins + 1)
    dr = float(rmax / nbins)
    for i in xrange(len(atoms)):
        for j in xrange(i + 1, len(atoms)):
            rij = dm[i][j]
            index = int(math.ceil(rij / dr))
            if index <= nbins:
                rdf[index] += 1

    # Normalize
    phi = len(atoms) / atoms.get_volume()
    norm = 2.0 * math.pi * dr * phi * len(atoms)

    dists = [0]
    for i in xrange(1, nbins + 1):
        rrr = (i - 0.5) * dr
        dists.append(rrr)
        rdf[i] /= (norm * ((rrr**2) + (dr**2) / 12.))

    return rdf, np.array(dists)


def get_nndist(atoms, distance_matrix):
    """
    Returns an estimate of the nearest neighbor bond distance
    in the supplied atoms object given the supplied distance_matrix.
    The estimate comes from the first peak in the radial distribution
    function.
    """
    rmax = np.sqrt(sum([sum(c**2) for c in atoms.cell])) / 2.
    nbins = 400
    rdf, dists = get_rdf(atoms, rmax, nbins, distance_matrix)
    i = 0
    while np.gradient(rdf)[i] >= 0:
        i += 1
    return dists[i]


def get_nnmat(atoms):
    """
    Calculate the nearest neighbor matrix as specified in
    S. Lysgaard et al., Top. Catal., 2014, 57 (1-4), pp 33-39

    Returns an array of average numbers of nearest neighbors
    the order is determined by self.elements.
    Example: self.elements = ["Cu", "Ni"]
    get_nnmat returns a single list [Cu-Cu bonds/N(Cu),
    Cu-Ni bonds/N(Cu), Ni-Cu bonds/N(Ni), Ni-Ni bonds/N(Ni)]
    where N(element) is the number of atoms of the type element
    in the atoms object.

    The distance matrix can be quite costly to calculate every
    time nnmat is required (and disk intensive if saved), thus
    it makes sense to calculate nnmat along with e.g. the
    potential energy and save it in atoms.info['data']['nnmat'].
    """
    if 'nnmat' in atoms.info['data']:
        return atoms.info['data']['nnmat']
    elements = sorted(set(atoms.get_chemical_symbols()))
    nnmat = np.zeros((len(elements), len(elements)))
    dm = get_distance_matrix(atoms)
    nndist = get_nndist(atoms, dm) + 0.2
    for i in xrange(len(atoms)):
        row = [j for j in xrange(len(elements))
               if atoms[i].symbol == elements[j]][0]
        neighbors = [j for j in xrange(len(dm[i])) if dm[i][j] < nndist]
        for n in neighbors:
            column = [j for j in xrange(len(elements))
                      if atoms[n].symbol == elements[j]][0]
            nnmat[row][column] += 1
    # divide by the number of that type of atoms in the structure
    for i, el in enumerate(elements):
        nnmat[i] /= len([j for j in range(len(atoms))
                         if atoms[int(j)].symbol == el])
    # makes a single list out of a list of lists
    nnlist = np.reshape(nnmat, (len(nnmat)**2))
    return nnlist
