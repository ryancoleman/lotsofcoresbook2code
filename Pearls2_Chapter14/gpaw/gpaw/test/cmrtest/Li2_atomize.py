from __future__ import print_function
import os

from ase.structure import molecule
from ase.io import read, write
from ase.parallel import rank

from gpaw import GPAW, restart

import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

import cmr
#from cmr.tools.log import Log
#cmr.logger.set_message_selection(Log.MSG_TYPE_ALL)

calculate = True
recalculate = True
analyse_from_dir = True # analyse local cmr files

upload_to_db = False  # upload cmr files to the database
analyse_from_db = False # analyse database

create_group = True # group calculations beloging to a given reaction

clean = False

if create_group: assert analyse_from_dir or analyse_from_db

if analyse_from_db: assert upload_to_db

symbol = 'Li'

# define the project in order to find it in the database!
project_id = 'my first project: atomize'

vacuum = 3.5

# calculator parameters
xc = 'LDA'
mode = 'lcao'
h = 0.20

cmr_params_template = {
    'db_keywords': [project_id],
    # add project_id also as a field to support search across projects
    'project_id': project_id,
    # user's tags
    'U_vacuum': vacuum,
    'U_xc': xc,
    'U_mode': mode,
    'U_h': h,
    }

if calculate:

    # molecule
    formula = symbol + '2'
    # set formula name to be written into the cmr file
    cmr_params = cmr_params_template.copy()
    cmr_params['U_formula'] = formula

    cmrfile = formula + '.cmr'

    system = molecule(formula)
    system.center(vacuum=vacuum)
    # Note: Molecules do not need broken cell symmetry!
    if 0:
        system.cell[1, 1] += 0.01
        system.cell[2, 2] += 0.02

    # Hund rule (for atoms)
    hund = (len(system) == 1)
    cmr_params['U_hund'] = hund

    # first calculation: LDA lcao
    calc = GPAW(mode=mode, xc=xc, h=h, hund=hund, txt=formula + '.txt')
    system.set_calculator(calc)
    e = system.get_potential_energy()
    # write gpw file
    calc.write(formula)
    # add total energy to users tags
    cmr_params['U_potential_energy'] = e
    # write the information 'as in' corresponding trajectory file
    # plus cmr_params into cmr file
    write(cmrfile, system, cmr_params=cmr_params)

    del calc

    # atom
    formula = symbol
    # set formula name to be written into the cmr file
    cmr_params = cmr_params_template.copy()
    cmr_params['U_formula'] = formula

    cmrfile = formula + '.cmr'

    system = molecule(formula)
    system.center(vacuum=vacuum)
    # Note: Li does not need broken cell symmetry! Many other atoms do!
    if 0:
        system.cell[1, 1] += 0.01
        system.cell[2, 2] += 0.02

    # Hund rule (for atoms)
    hund = (len(system) == 1)
    cmr_params['U_hund'] = hund

    # first calculation: LDA lcao
    calc = GPAW(mode=mode, xc=xc, h=h, hund=hund, txt=formula + '.txt')
    system.set_calculator(calc)
    e = system.get_potential_energy()
    # write gpw file
    calc.write(formula)
    # add total energy to users tags
    cmr_params['U_potential_energy'] = e
    # write the information 'as in' corresponding trajectory file
    # plus cmr_params into cmr file
    write(cmrfile, system, cmr_params=cmr_params)

    del calc

if recalculate:

    # now calculate PBE energies on LDA orbitals

    # molecule
    formula = symbol + '2'
    system, calc = restart(formula, txt=None)

    ediff = calc.get_xc_difference('PBE')

    cmrfile = formula + '.cmr'

    # add new results to the cmrfile
    data = cmr.read(cmrfile)
    data.set_user_variable('U_potential_energy_PBE', data['U_potential_energy'] + ediff)
    data.write(cmrfile)

    del calc

    # atom
    formula = symbol
    system, calc = restart(formula, txt=None)

    ediff = calc.get_xc_difference('PBE')

    cmrfile = formula + '.cmr'

    # add new results to the cmrfile
    data = cmr.read(cmrfile)
    data.set_user_variable('U_potential_energy_PBE', data['U_potential_energy'] + ediff)
    data.write(cmrfile)

    del calc

if analyse_from_dir:

    # analyze the results from cmr files in the local directory
    from cmr.ui import DirectoryReader

    # read all compounds in the project with lcao and LDA orbitals
    reader = DirectoryReader(directory='.', ext='.cmr')
    all = reader.find(name_value_list=[('U_mode', 'lcao'), ('U_xc', 'LDA')],
                      keyword_list=[project_id])
    if rank == 0:
        print('results from cmr files in the local directory')
    # print requested results
    # column_length=0 aligns data in the table (-1 : data unaligned is default)
    all.print_table(column_length=0,
                    columns=['U_formula', 'U_vacuum',
                             'U_xc', 'U_h', 'U_hund',
                             'U_potential_energy', 'U_potential_energy_PBE',
                             'ase_temperature'])

    # access the results directly and calculate atomization energies
    f2 = symbol + '2'
    f1 = symbol

    if rank == 0:

        # results are accesible only on master rank

        r2 = all.get('U_formula', f2)
        r1 = all.get('U_formula', f1)

        # calculate atomization energies (ea)
        ea_LDA = 2 * r1['U_potential_energy'] - r2['U_potential_energy']
        ea_PBE = 2 * r1['U_potential_energy_PBE'] - r2['U_potential_energy_PBE']
        print('atomization energy [eV] ' + xc + ' = ' + str(ea_LDA))
        print('atomization energy [eV] PBE = ' + str(ea_PBE))

        if create_group:
            # ea_LDA and ea_PBE define a group
            group = cmr.create_group();
            group.add(r1['db_hash']);
            group.add(r2['db_hash']);
            group.set_user_variable('U_ea_LDA', ea_LDA)
            group.set_user_variable('U_ea_PBE', ea_PBE)
            group.set_user_variable('U_description', 'atomization energy [eV]')
            group.set_user_variable('U_reaction', '2 * ' + symbol + ' - ' + symbol + '2')
            group.set_user_variable('db_keywords', [project_id])
            group.set_user_variable('project_id', project_id)
            group.write(symbol + '2_atomize_from_dir.cmr');

    if True:

        all = reader.find(keyword_list=[project_id])
        
        if rank == 0:
            print('contents of the cmr files present in the local directory')
        # print requested results
        # column_length=0 aligns data in the table (-1 : data unaligned is default)
        all.print_table(column_length=0,
                        columns=['U_formula', 'U_vacuum',
                                 'U_xc', 'U_h', 'U_hund',
                                 'U_potential_energy', 'U_potential_energy_PBE',
                                 'ase_temperature', 'U_reaction', 'U_ea_LDA', 'U_ea_PBE', 'U_description'])

if upload_to_db:

    # upload cmr files to the database

    if rank == 0:
        os.system('cmr --commit ' + symbol + '*.cmr')

if analyse_from_db:

    # analyze the results from the database
    # analysis can only be performed on rank 0!!
    from cmr.ui import DBReader
    reader = DBReader()
    all = reader.find(name_value_list=[('U_mode', 'lcao'),
                                       ('U_xc', 'LDA'),
                                       #('db_user', '')
                                       ],
                      keyword_list=[project_id])
    
    if rank == 0:
        print('results from the database')
    # print requested results
    # column_length=0 aligns data in the table (-1 : data unaligned is default)
    all.print_table(column_length=0,
                    columns=['U_formula', 'U_vacuum',
                             'U_xc', 'U_h', 'U_hund',
                             'U_potential_energy', 'U_potential_energy_PBE',
                             'ase_temperature'])

    # access the results directly and calculate atomization energies
    f2 = symbol + '2'
    f1 = symbol

    # results are accesible only on master rank
    r1 = all.get('U_formula', f1)
    r2 = all.get('U_formula', f2)

    # check if results were successfully retrieved, otherwise we have to wait
    if r1 is None or r2 is None:
        print("Results are not yet in the database. Wait, and try again.")
    else:
        # calculate atomization energies (ea)
        ea_LDA = 2 * r1['U_potential_energy'] - r2['U_potential_energy']
        ea_PBE = 2 * r1['U_potential_energy_PBE'] - r2['U_potential_energy_PBE']
        if rank == 0:
            print('atomization energy [eV] ' + xc + ' = ' + str(ea_LDA))
            print('atomization energy [eV] PBE = ' + str(ea_PBE))

        if create_group:
            # ea_LDA and ea_PBE define a group
            group = cmr.create_group();
            group.add(r1['db_hash']);
            group.add(r2['db_hash']);
            group.set_user_variable('U_ea_LDA', ea_LDA)
            group.set_user_variable('U_ea_PBE', ea_PBE)
            group.set_user_variable('U_description', 'atomization energy [eV] (from database)')
            group.set_user_variable('U_reaction', '2 * ' + symbol + ' - ' + symbol + '2')
            group.set_user_variable('db_keywords', [project_id])
            group.set_user_variable('project_id', project_id)
            group.write(symbol + '2_atomize_from_db.cmr');
            group.write(".cmr");

    if True:

        all = reader.find(keyword_list=[project_id])
        
        if rank == 0:
            print('contents of the database')
        # print requested results
        # column_length=0 aligns data in the table (-1 : data unaligned is default)
        all.print_table(column_length=0,
                        columns=['U_formula', 'U_vacuum',
                                 'U_xc', 'U_h', 'U_hund',
                                 'U_potential_energy', 'U_potential_energy_PBE',
                                 'ase_temperature', 'U_reaction', 'U_ea_LDA', 'U_ea_PBE', 'U_description'])

if clean:

    if rank == 0:
        for file in [symbol + '.cmr', symbol + '.gpw', symbol + '.txt',
                     symbol + '2.cmr', symbol + '2.gpw', symbol + '2.txt',
                     symbol + '2_atomize_from_dir.cmr',
                     symbol + '2_atomize_from_db.cmr']:
            if os.path.exists(file): os.unlink(file)


