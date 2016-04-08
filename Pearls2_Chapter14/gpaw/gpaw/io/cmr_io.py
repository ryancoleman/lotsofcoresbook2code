# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import os
import time
import xml.sax

import numpy as np

import ase
from ase.version import version as ase_version
import gpaw
from gpaw.version import version_base as gpaw_version

try:
    #new style cmr io
    import cmr
    from cmr.base.converter import Converter
    from cmr.io import Flags
    from cmr.tools.functions import create_db_filename as cdbfn
    from cmr.definitions import CALCULATOR_GPAW

    def get_reader(name):
        reader = cmr.read(name, mode=Flags.READ_MODE_ORIGINAL_TYPE)
        return reader
        
    def get_writer():
        return Converter.get_xml_writer(CALCULATOR_GPAW, calculator_version=gpaw_version)
    
    def create_db_filename(param, ext=".db"):
        return cdbfn(param, ext=ext)
    
except:    
    #old style cmr io
    import cmr
    from cmr import create_db_filename as cdbfn
    from cmr.io import XMLData
    from cmr.io import READ_DATA
    from cmr.io import EVALUATE
    from cmr.io import WRITE_CONVERT
    from cmr.io import CONVERTION_ORIGINAL
    from cmr.static import CALCULATOR_GPAW
    
    def create_db_filename(param, ext="ignored"):
        return cdbfn()

    def get_reader(name):
        reader = cmr.read(name,
                          read_mode=READ_DATA,
                          evaluation_mode=EVALUATE,
                          convertion_mode=CONVERTION_ORIGINAL)
        return reader
        
    def get_writer():
        data = XMLData()
        data.set_calculator_name(CALCULATOR_GPAW)
        data.set_write_mode(WRITE_CONVERT)
        return data

class Writer:
    """ This class is a wrapper to the db output writer
    and intended to be used with gpaw
    """
    def __init__(self, filename, comm=None):
        self.comm = comm # for possible future use
        self.verbose = False
        self.data = get_writer()
        self.split_array = None #used when array is not filled at once
        self.dimensions = {}
        self.filename = filename
        self.cmr_params = {}

        uname = os.uname()
        self.data['user']=os.getenv('USER', '???')
        self.data['date']=time.asctime()

        self.data['architecture']=uname[4]
        self.data['ase_dir']=os.path.dirname(ase.__file__)
        self.data['ase_version']=ase_version
        self.data['numpy_dir']=os.path.dirname(np.__file__)
        self.data['gpaw_dir']=os.path.dirname(gpaw.__file__)
        self.data["db_calculator_version"] = gpaw_version
        self.data['calculator']="gpaw"
        self.data['location']=uname[1]

        
    def dimension(self, name, value):
        self.dimensions[name]=value
        self.data[name]=value
        if self.verbose:
            print("dimension: ", name, value)

    def __setitem__(self, name, value):
        """ sets the value of a variable in the db-file. Note that only
        values that were defined in the cmr-schema are written.
        (User defined values have to be added with set_user_variable(name, value)
        IMPORTANT: CMR does not support None values for numbers, therefore all variables
        with value None are ignored."""
        if self.verbose:
            print("name value:", name, value)
        if name == "GridSpacing" and value == "None":
            return
        if not value is None:
            self.data[name]=value
        
    def _get_dimension(self, array):
        """retrieves the dimension of a multidimensional array
        by using then len() function until fail
        """
        indent = ""
        measured_dimensions = []
        while True:
            try:
                measured_dimensions.append(len(array))
                if self.verbose:
                    print(indent+"Length:", len(array))
                    indent += " "
                array = array[0]
            except IndexError:
                break
            except TypeError:
                break
        return measured_dimensions

    def _close_array(self):
        """if an array is filled with fill then we don't know
        the exact end of the array. Therefore we check before
        adding a new variable if the end was reached."""
        if self.split_array is None:
            return
        
        measured_dimensions = self._get_dimension(self.split_array[2])
        if self.verbose:
            print("Dimensions:         ", self.split_array[1])
            print("Mesured Dimensions: ", measured_dimensions)
            
        #make the array fit (only fixes the 1st dimension)
        res = np.array(self.split_array[2]).reshape(self.split_array[1])
        self.data[self.split_array[0]] = res

    def add(self, name, shape, array=None, dtype=None, units=None,
            parallel=False, write=True):
        self._close_array()
        if self.verbose:
            print("add:", name, shape, array, dtype, units)
        if array is None:
            dimension = []
            for a in shape:
                dimension.append(self.dimensions[a])
            self.split_array = (name, dimension, [])
        else:
            self.data[name]=array

    def fill(self, array, *indices, **kwargs):
        if self.verbose:
            print("fill (", len(array),"):", array)
        self.split_array[2].append(array)

    def set_db_copy_settings(self, make_db_copy, private):
        """deprecated: This method is going to be removed"""
        print("Warning: gpaw.cmr.readwriter.set_db_copy_settings is deprecated.")
        print("Please update gpaw and CMR")
        pass
        #if self.verbose:
        #    print "set_db_copy_settings", make_db_copy, private
        #self.db_copy=make_db_copy
        #self.private=private
        #self.data.set_db_copy_settings(make_db_copy, private)

    def write(self, string, db, private, **kwargs):
        print("Warning: gpaw.cmr.readwriter.write is deprecated.")
        print("Please update gpaw and CMR")
        pass
        #if self.verbose:
        #    print "write():", string
        #self.data.write(string, db, private, **kwargs)

    def write_additional_db_params(self, cmr_params):
        """writes the user variables and also sets the write attributes for
        the output file"""
        self.cmr_params = cmr_params.copy()
        #cmr.set_params_to_xml_data(self.data, cmr_params)
        
    def close(self):
        if self.verbose:
            print("close()")
        self._close_array()
        if self.cmr_params.has_key("ase_atoms_var"):
            ase_vars = self.cmr_params["ase_atoms_var"]
            for key in ase_vars:
                self.data.set_user_variable(key, ase_vars[key])
            self.cmr_params.pop("ase_atoms_var")
        if self.filename==".db" or self.filename==".cmr":
            # Note: 
            #      .cmr files can currently not be uploaded to the database therefore 
            #      it defaults to .db until supported
            self.cmr_params["output"]=create_db_filename(self.data, ext=".db")
        else:
            self.cmr_params["output"]=self.filename
        try:
            cmr.runtime.pause_ase_barriers(True)
            self.data.write(self.cmr_params, ase_barrier=False)
        except TypeError:
            # for compatibility with older CMR versions:
            self.data.write(self.cmr_params)
        cmr.runtime.pause_ase_barriers(False)
        return [self.data.get_hash()]


class Reader:
    """ This class allows gpaw to access
    to read a db-file
    """
    def __init__(self, name, comm):
        self.verbose = False
        self.reader = self.parameters = get_reader(name)

    def dimension(self, name):
        return self.reader[name]
    
    def __getitem__(self, name):
        if name=='version' and not self.reader.has_key('version') \
            and self.reader.has_key('db_calculator_version'):
                return self.reader['db_calculator_version']
        return self.reader[name]

    def has_array(self, name):
        return self.reader.has_key(name)
    
    def get(self, name, *indices, **kwargs):
        if self.verbose:
            print("incides", indices)
        result = self.reader[name]
        if indices!=():
            for a in indices:
                result = result[a]
            return result
        
        #gpaw wants expressions evaluated
        if type(result)==str or type(result)==unicode:
            try:
                if self.verbose:
                    print("Converting ", result)
                result = eval(result, {})
            except (SyntaxError, NameError):
                pass
        return result
    
    def get_reference(self, name, indices, length=None):
        result = self.reader[name]
        if indices!=():
            for a in indices:
                result = result[a]
            return result
    
    def get_file_object(self, name, indices):
        result = self.reader.retrieve_file(name)
        if indices!=():
            for a in indices:
                result = result[a]
        return result
        

    def get_parameters(self):
        return self.reader.keys()

    def close(self):
        pass


