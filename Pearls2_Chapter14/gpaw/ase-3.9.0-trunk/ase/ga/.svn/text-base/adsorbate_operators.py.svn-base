"""Adsorbate operators that adds an adsorbate to the surface
of a particle or given structure, using a supplied list of sites."""
import numpy as np
import random
from ase import Atoms, Atom
from ase.structure import molecule
from ase.ga.offspring_creator import OffspringCreator
from ase.calculators.neighborlist import NeighborList as aseNeighborList


class AdsorbateOperator(OffspringCreator):
    """Base class for all operators that add, move or remove adsorbates.

    Don't use this operator directly!"""

    def __init__(self, adsorbate):
        OffspringCreator.__init__(self)
        self.adsorbate = self.convert_adsorbate(adsorbate)
        self.descriptor = 'AdsorbateOperator'

    def get_new_individual(self, parents):
        raise NotImplementedError

    def add_adsorbate(self, atoms, sites_list, min_adsorbate_distance=1.5):
        """Adds the adsorbate in self.adsorbate to the supplied atoms
        object at the first free site in the specified sites_list. A site
        is free if no other adsorbates can be found in a sphere of radius
        min_adsorbate_distance around the chosen site.

        Parameters:

        atoms: Atoms object
            the atoms object that the adsorbate will be added to

        sites_list: list
            a list of dictionaries, each dictionary should be of the
            following form:
            {'height': h, 'normal': n, 'adsorbate_position': ap,
            'site': si, 'surface': su}

        min_adsorbate_distance: float
            the radius of the sphere inside which no other
            adsorbates should be found
        """
        i = 0
        while self.is_site_occupied(atoms, sites_list[i],
                                    min_adsorbate_distance):
            i += 1
            if i >= len(sites_list):
                return False
        site = sites_list[i]

        ## Make the correct position
        height = site['height']
        normal = np.array(site['normal'])
        pos = np.array(site['adsorbate_position']) + normal * height

        ## Rotate the adsorbate according to the normal
        ads = self.adsorbate.copy()
        if len(ads) > 1:
            avg_pos = np.average(ads[1:].positions, 0)
            ads.rotate(avg_pos - ads[0].position, normal)
        ads.translate(pos - ads[0].position)

        atoms.extend(ads)

        return True

    def remove_adsorbate(self, atoms, sites_list, for_move=False):
        """Removes an adsorbate from the atoms object at the first occupied
        site in sites_list. If no adsorbates can be found, one will be
        added instead.
        """
        i = 0
        while not self.is_site_occupied(atoms, sites_list[i],
                                        min_adsorbate_distance=0.2):
            ## very small min_adsorbate_distance used for testing
            i += 1
            if i >= len(sites_list):
                if for_move:
                    return False
                print('removal not possible will add instead')
                return self.add_adsorbate(atoms, sites_list)
        sites_list[i]['occupied'] = False
        site = sites_list[i]

        ## Make the correct position
        height = site['height']
        normal = np.array(site['normal'])
        pos = np.array(site['adsorbate_position']) + normal * height

        ads_ind = self.get_adsorbate_indices(atoms, pos)
        ads_ind.sort(reverse=True)

        # print('removing', ads_ind, [atoms[j].symbol for j in ads_ind], pos)
        for k in ads_ind:
            atoms.pop(k)

        return True

    def get_adsorbate_indices(self, atoms, position):
        """Returns the indices of the adsorbate at the supplied position"""
        dmin = 1000.
        for a in atoms:
            d = np.linalg.norm(a.position - position)
            if d < dmin:
                dmin = d
                ind = a.index

        mbl = 1.5  # max_bond_length
        nl = aseNeighborList([mbl / 2. for i in atoms],
                             skin=0.0, self_interaction=False)
        nl.update(atoms)

        return list(set(self._get_indices_in_adsorbate(atoms, nl, ind, [])))

    def _get_indices_in_adsorbate(self, atoms, neighborlist,
                                  index, molecule_indices=[]):
        """Internal recursive function that help
        determine adsorbate indices"""
        mi = molecule_indices
        nl = neighborlist
        mi.append(index)
        neighbors, _ = nl.get_neighbors(index)
        ads_syms = self.adsorbate.copy().get_chemical_symbols()
        for n in neighbors:
            if int(n) not in mi:
                if atoms[int(n)].symbol in ads_syms:
                    mi = self._get_indices_in_adsorbate(atoms, nl, n, mi)
        return mi

    def is_site_occupied(self, atoms, site, min_adsorbate_distance):
        """Returns True if the site on the atoms object is occupied by
        creating a sphere of radius min_adsorbate_distance and checking
        that no other adsorbate is inside the sphere."""
        if site['occupied'] == 1:
            return True
        ads = self.adsorbate.get_chemical_symbols()
        height = site['height']
        normal = np.array(site['normal'])
        pos = np.array(site['adsorbate_position']) + normal * height
        dists = [np.linalg.norm(pos - a.position)
                 for a in atoms if a.symbol in ads]
        for d in dists:
            if d < min_adsorbate_distance:
                # print('under min d', d, pos)
                site['occupied'] = 1
                return True
        return False

    @classmethod
    def convert_adsorbate(cls, adsorbate):
        """Converts the adsorbate to an Atoms object"""
        if isinstance(adsorbate, Atoms):
            ads = adsorbate
        elif isinstance(adsorbate, Atom):
            ads = Atoms([adsorbate])
        else:
            # Hope it is a useful string or something like that
            if adsorbate == 'CO':
                # CO otherwise comes out as OC - very inconvenient
                ads = molecule(adsorbate, symbols=adsorbate)
            else:
                ads = molecule(adsorbate)
        ads.translate(-ads[0].position)
        return ads


class AddAdsorbate(AdsorbateOperator):
    """
    Use this operator to add adsorbates to the surface.

    Supply a list of adsorption_sites of the form:
    [{'adsorbate_position':[.., .., ..], 'normal':surface_normal_vector,
    'height':height, 'site':site, 'surface':surface}, {...}, ...]

    The adsorbate will be positioned at:
    adsorbate_position + surface_normal_vector * height
    The site and surface parameters are supplied to be able to keep
    track of which sites and surfaces are filled - useful to determine
    beforehand or know afterwards.

    If the surface is allowed to change during the algorithm run,
    a list of adsorbate sites should not be supplied. It would instead
    be generated for every case, however this has not been implemented
    here yet.

    Site and surface preference can be supplied. If both are supplied site
    will be considered first.
    """
    def __init__(self, adsorbate,
                 minimum_adsorbate_distance=2.,
                 adsorption_sites=None,
                 site_preference=None,
                 surface_preference=None):
        AdsorbateOperator.__init__(self, adsorbate)
        self.descriptor = 'AddAdsorbate'

        self.min_adsorbate_distance = minimum_adsorbate_distance

        if adsorption_sites is None:
            raise NotImplementedError
        ## Adding 0's to the end of all sites to specify not filled
        for s in adsorption_sites:
            s.update({'occupied': 0})
        self.adsorption_sites = adsorption_sites
        self.site_preference = site_preference
        self.surface_preference = surface_preference

        self.num_muts = 1

    def get_new_individual(self, parents):
        """Returns the new individual as an atoms object"""
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        for atom in f:
            indi.append(atom)

        for _ in xrange(self.num_muts):
            random.shuffle(self.adsorption_sites)

            if self.surface_preference is not None:
                def func(x):
                    return x['surface'] == self.surface_preference
                self.adsorption_sites.sort(key=func, reverse=True)

            if self.site_preference is not None:
                def func(x):
                    return x['site'] == self.site_preference
                self.adsorption_sites.sort(key=func, reverse=True)

            added = self.add_adsorbate(indi, self.adsorption_sites,
                                       self.min_adsorbate_distance)
            if not added:
                break

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


class RemoveAdsorbate(AdsorbateOperator):
    """This operator removes an adsorbate from the surface. It works
    exactly (but doing the opposite) as the AddAdsorbate operator."""
    def __init__(self, adsorbate,
                 adsorption_sites=None,
                 site_preference=None,
                 surface_preference=None):
        AdsorbateOperator.__init__(self, adsorbate)
        self.descriptor = 'RemoveAdsorbate'

        if adsorption_sites is None:
            raise NotImplementedError
        # Adding 0's to the end of all sites to specify not filled
        for s in adsorption_sites:
            s.update({'occupied': 0})
        self.adsorption_sites = adsorption_sites
        self.site_preference = site_preference
        self.surface_preference = surface_preference

        self.num_muts = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        for atom in f:
            indi.append(atom)

        for _ in xrange(self.num_muts):
            random.shuffle(self.adsorption_sites)

            if self.surface_preference is not None:
                def func(x):
                    return x['surface'] == self.surface_preference
                self.adsorption_sites.sort(key=func, reverse=True)

            if self.site_preference is not None:
                def func(x):
                    return x['site'] == self.site_preference
                self.adsorption_sites.sort(key=func, reverse=True)

            removed = self.remove_adsorbate(indi, self.adsorption_sites)

            if not removed:
                break

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))


class MoveAdsorbate(AdsorbateOperator):
    """This operator removes an adsorbate from the surface and adds it
    again at a different position, i.e. effectively moving the adsorbate."""
    def __init__(self, adsorbate,
                 minimum_adsorbate_distance=2.,
                 adsorption_sites=None,
                 site_preference_from=None,
                 surface_preference_from=None,
                 site_preference_to=None,
                 surface_preference_to=None):
        AdsorbateOperator.__init__(self, adsorbate)
        self.descriptor = 'MoveAdsorbate'

        self.min_adsorbate_distance = minimum_adsorbate_distance

        if adsorption_sites is None:
            raise NotImplementedError
        ## Adding 0's to the end of all sites to specify not filled
        for s in adsorption_sites:
            s.update({'occupied': 0})
        self.adsorption_sites = adsorption_sites
        self.site_preference_from = site_preference_from
        self.surface_preference_from = surface_preference_from
        self.site_preference_to = site_preference_to
        self.surface_preference_to = surface_preference_to

        self.num_muts = 1

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f)
        indi.info['data']['parents'] = [f.info['confid']]

        for atom in f:
            indi.append(atom)

        for _ in xrange(self.num_muts):
            random.shuffle(self.adsorption_sites)
            if self.surface_preference_from is not None:
                def func(x):
                    return x['surface'] == self.surface_preference_from
                self.adsorption_sites.sort(key=func, reverse=True)

            if self.site_preference_from is not None:
                def func(x):
                    return x['site'] == self.site_preference_from
                self.adsorption_sites.sort(key=func, reverse=True)

            removed = self.remove_adsorbate(indi, self.adsorption_sites,
                                            for_move=True)

            random.shuffle(self.adsorption_sites)
            if self.surface_preference_to is not None:
                def func(x):
                    return x['surface'] == self.surface_preference_to
                self.adsorption_sites.sort(key=func, reverse=True)

            if self.site_preference_to is not None:
                def func(x):
                    return x['site'] == self.site_preference_to
                self.adsorption_sites.sort(key=func, reverse=True)

            added = self.add_adsorbate(indi, self.adsorption_sites,
                                       self.min_adsorbate_distance)

            if (not removed) or (not added):
                break

        return (self.finalize_individual(indi),
                self.descriptor + ': {0}'.format(f.info['confid']))
