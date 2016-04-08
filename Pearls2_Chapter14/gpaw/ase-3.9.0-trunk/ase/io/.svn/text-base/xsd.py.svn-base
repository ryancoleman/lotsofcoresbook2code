import numpy as np
import xml.etree.ElementTree as ET

from ase import Atoms


def read_xsd(filename):
    tree = ET.parse(filename)
    root = tree.getroot()

    atomtreeroot = root.find('AtomisticTreeRoot')
    symmetrysystem = atomtreeroot.find('SymmetrySystem')
    mappingset = symmetrysystem.find('MappingSet')
    mappingfamily = mappingset.find('MappingFamily')
    system = mappingfamily.find('IdentityMapping')

    coords = list()
    cell = list()
    formula = str()

    for atom in system:
        if atom.tag == 'Atom3d':
            symbol = atom.get('Components')
#            symbols.append(symbol)
            formula += symbol
        
            xyz = atom.get('XYZ')
            coord = [float(coord) for coord in xyz.split(',')]
            coords.append(coord)
        elif atom.tag == 'SpaceGroup':
            avec = [float(vec) for vec in atom.get('AVector').split(',')]
            bvec = [float(vec) for vec in atom.get('BVector').split(',')]
            cvec = [float(vec) for vec in atom.get('CVector').split(',')]
        
            cell.append(avec)
            cell.append(bvec)
            cell.append(cvec)

#    print formula
#    print np.array(coords)
#    print np.array(cell)

    atoms = Atoms(formula, cell=cell)
    atoms.set_scaled_positions(coords)
    return atoms
