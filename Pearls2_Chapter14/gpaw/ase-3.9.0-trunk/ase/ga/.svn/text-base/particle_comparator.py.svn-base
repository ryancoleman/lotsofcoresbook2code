"""Comparators originally meant to be used with particles"""
import numpy as np
from ase.ga.utilities import get_nnmat


class NNMatComparator(object):
    """Use the nearest neighbor matrix to determine differences
    in the distribution (and to a slighter degree structure)
    of atoms. As specified in
    S. Lysgaard et al., Top. Catal., 57 (1-4), pp 33-39, (2014)"""
    def __init__(self, d=0.2, elements=[]):
        self.d = d
        self.elements = elements

    def looks_like(self, a1, a2):
        """ Return if structure a1 or a2 are similar or not. """
        elements = self.elements
        if elements == []:
            elements = sorted(set(a1.get_chemical_symbols()))
        a1, a2 = a1.copy(), a2.copy()
        del a1[[a.index for a in a1 if a.symbol not in elements]]
        del a2[[a.index for a in a2 if a.symbol not in elements]]

        nnmat_a1 = get_nnmat(a1)
        nnmat_a2 = get_nnmat(a2)

        diff = np.linalg.norm(nnmat_a1 - nnmat_a2)

        if diff < self.d:
            return True
        else:
            return False


class AdsorptionSitesComparator(object):
    """Compares the metal atoms in the adsorption sites and returns True
    if less than min_diff_adsorption_sites of the sites with adsorbates
    consist of different atoms.

    Ex:
    a1.info['data']['adsorbates_site_atoms'] =
    [('Cu','Ni'),('Cu','Ni'),('Ni'),('Ni')]

    a2.info['data']['adsorbates_site_atoms'] =
    [('Cu','Ni'),('Ni','Ni', 'Ni'),('Ni'),('Ni')]

    will have a difference of 2:
    (2*('Cu','Ni')-1*('Cu','Ni')=1, 1*('Ni','Ni','Ni')=1, 2*('Ni')-2*('Ni')=0)

    """
    def __init__(self, min_diff_adsorption_sites=2):
        self.min_diff_adsorption_sites = min_diff_adsorption_sites

    def looks_like(self, a1, a2):
        s = 'adsorbates_site_atoms'
        if not all([(s in a.info['data'] and
                     a.info['data'][s] != [])
                    for a in [a1, a2]]):
            return False

        counter = {}
        for asa in a1.info['data'][s]:
            t_asa = tuple(sorted(asa))
            if t_asa not in counter.keys():
                counter[t_asa] = 1
            else:
                counter[t_asa] += 1

        for asa in a2.info['data'][s]:
            t_asa = tuple(sorted(asa))
            if t_asa not in counter.keys():
                counter[t_asa] = -1
            else:
                counter[t_asa] -= 1

        diffs = len([k for k, v in counter.iteritems() if v != 0])

        if diffs < self.min_diff_adsorption_sites:
            return True

        return False


class AdsorptionMetalsComparator(object):
    """Compares the number of adsorbate-metal bonds and returns True if the
    number for a1 and a2 differs by less than the supplied parameter
    ``same_adsorption_number``

    Ex:
    a1.info['data']['adsorbates_bound_to'] = {'Cu':1, 'Ni':3}
    a2.info['data']['adsorbates_bound_to'] = {'Cu':.5, 'Ni':3.5}
    will have a difference of .5 in both elements:
    """
    def __init__(self, same_adsorption_number):
        self.same_adsorption_number = same_adsorption_number

    def looks_like(self, a1, a2):
        s = 'adsorbates_bound_to'
        if not all([(s in a.info['data'] and
                     any(a.info['data'][s].values()))
                    for a in [a1, a2]]):
            return False

        diffs = [a1.info['data'][s][k] - a2.info['data'][s][k]
                 for k in a1.info['data'][s].keys()]
        for d in diffs:
            if abs(d) < self.same_adsorption_number:
                return True
        return False
