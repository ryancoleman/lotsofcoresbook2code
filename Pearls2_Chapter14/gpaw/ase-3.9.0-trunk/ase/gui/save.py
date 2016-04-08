# encoding: utf-8
"""save.py - Window for saving one or more configurations.
"""

import gtk
from gettext import gettext as _
from copy import copy
from ase.gui.widgets import pack, cancel_apply_ok, oops, help
from ase.io import string2index, write
import numpy as np
import sys
import os

class SaveWindow(gtk.Window):
    # List of valid file types - translation occurs when *using* this list!
    # Name, suffix, is_graphics, multiimage
    filetypes = [('Automatic', None, None, None),
                 ('XYZ file', 'xyz', False, True),
                 ('ASE trajectory', 'traj', False, True),
                 ('PDB file', 'pdb', False, True),
                 ('Gaussian cube file', 'cube', False, False),
                 ('Python script', 'py', False, True),
                 ('VNL file', 'vnl', False, False),
                 ('Portable Network Graphics', 'png', True, False),
                 ('Persistence of Vision', 'pov', True, False),
                 ('Encapsulated PostScript', 'eps', True, False),
                 ('FHI-aims geometry input', 'in', False, False),
                 ('CASTEP geom file', 'cell', False, True),
                 ('VASP geometry input', 'POSCAR', False, False),
                 ('ASE bundle trajectory', 'bundle', False, True),
                 ('cif file', 'cif', False, True),
                 ]
    
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui
        self.setup()
        self.show()

    def setup(self):
        self.set_title(_('Save ...'))
        vbox = gtk.VBox()
        self.add(vbox)
        confframe = gtk.Frame()  # For selecting which configs to save
        pack(vbox, confframe)
        vbox2 = gtk.VBox()
        confframe.add(vbox2)
        vbox2.show()
        self.radio_thisconf = gtk.RadioButton(
            None, _('Save this configuration only (#%s).') % self.frame)
        self.radio_allconf = gtk.RadioButton(self.radio_thisconf, 
                                             _('Save all configurations.'))
        self.radio_someconf = gtk.RadioButton(self.radio_thisconf, 
                                             _('Save some configurations: '))
        self.whichconf = gtk.Entry(max=20)
        pack(vbox2, [self.radio_thisconf])    
        pack(vbox2, [self.radio_allconf])
        pack(vbox2, [self.radio_someconf, self.whichconf])
        if self.gui.images.nimages <= 1:
            self.radio_thisconf.set_active(True)
            self.radio_allconf.set_sensitive(False)
            self.radio_someconf.set_sensitive(False)
            self.whichconf.set_sensitive(False)
        else:
            self.radio_allconf.set_active(True)
        self.chooser = gtk.FileChooserWidget(gtk.FILE_CHOOSER_ACTION_SAVE)
        try:
            fname = sys.argv[1]
        except IndexError:
            fname = ""
            self.chooser.set_current_folder(os.getcwd())
        self.chooser.set_current_name(fname)
        self.add_filters(self.chooser)
        self.old_chooser_name = self.chooser.get_filter().get_name()
        self.chooser.connect('notify::filter', self.filter_changed)
        
        pack(vbox, self.chooser, expand=True)
        
        # Add buttons
        self.buttons = gtk.HButtonBox()
        savebut = gtk.Button(stock=gtk.STOCK_SAVE)
        savebut.connect('clicked', self.save)
        cancelbut = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancelbut.connect('clicked', lambda x: self.destroy())
        for w in (savebut, cancelbut):
            self.buttons.pack_start(w, 0, 0)
            w.show()
        pack(vbox, [self.buttons], end=True, bottom=True)
        vbox.show()

    def filter_changed(self, *args):
        #print "notify::filter called: " + str(args)
        newname = self.chooser.get_filter().get_name()
        if newname == self.old_chooser_name:
            return  # Nothing has happened.
        oldsuffix = self.name_to_suffix[self.old_chooser_name]
        self.old_chooser_name = newname  # Remember it.
        newsuffix = self.name_to_suffix[newname]
        if newsuffix is None:
            # Change to Automatic - do nothing
            return
        filename = self.chooser.get_filename()
        fileprefix, filesuffix = os.path.splitext(filename)
        if oldsuffix is None:
            # Change away from Automatic, any valid suffix will be changed.
            if filesuffix == '':
                # No old suffix, append new suffix
                newfilename = filename + '.' + newsuffix
            elif filesuffix[1:] in self.name_to_suffix.values():
                # Old suffix is valid, replace with new.
                newfilename = fileprefix + '.' + newsuffix
            else:
                # Old suffix is weird, perhaps a . in filename.  Append the right suffix
                newfilename = filename + '.' + newsuffix
        elif oldsuffix == filesuffix[1:]:
            # Change away from valid suffix, replace old with new.
            newfilename = fileprefix + '.' + newsuffix
        else:
            # Old suffix does not match old file type - do nothing, perhaps user knows 
            # what he is doing...
            return
        # Change the filename
        self.chooser.set_current_name(newfilename)

    def save(self, dummy):
        "The user has pressed the SAVE button."
        filename = self.chooser.get_filename()
        if not filename or filename == "<<filename>>":
            oops("Please specify a file name")
            return

        # Check file type
        suffix = os.path.splitext(filename)[1][1:]
        if 'POSCAR' in filename or 'CONTCAR' in filename:
            suffix = 'POSCAR'
        if suffix == '':
            # No suffix chosen
            filt = self.chooser.get_filter().get_name()
            suffix = self.name_to_suffix[filt]
            if suffix is None:
                oops("Specify file type by giving a suffix or selecting a file type.")
                return
            else:
                filename = filename + '.' + suffix
        else:
            # Suffix given - check that it is not in conflict with selected file type.
            filt = self.chooser.get_filter().get_name()
            suffix2 = self.name_to_suffix[filt]
            if suffix2 is not None and suffix != suffix2:
                oops("Your filename suffix conflicts with the file type you have selected.")
                return
        if suffix not in self.name_to_suffix.values():
            oops("Unknown file suffix "+suffix)
            return
        
        # We now have a consistent file name with an allowed suffix.
        # Find out which images we want to save.
        if self.radio_thisconf.get_active():
            indices = [self.gui.frame]
        elif self.radio_allconf.get_active():
            indices = range(self.gui.images.nimages)
        elif self.radio_someconf.get_active():
            txt = self.whichconf.get_text()
            if not txt:
                oops("Please specify which images to save.")
                return
            try:
                slice = string2index(txt)
            except ValueError:
                oops("ERROR: Failed to parse image specification '%s'" % (txt,))
                return
            indices = range(self.gui.images.nimages)[slice]
            if isinstance(indices, int):
                indices = [indices]
        else:
            raise RuntimeError("No radio button selected - should not be possible!")

        # Now we are ready to write the file!
        extra = {}
        remove_hidden = False
        if self.is_graphics[suffix]:
            bbox = np.empty(4)
            size = np.array([self.gui.width, self.gui.height]) / self.gui.scale
            bbox[0:2] = np.dot(self.gui.center, self.gui.axes[:, :2]) - size / 2
            bbox[2:] = bbox[:2] + size
            extra['rotation'] = self.gui.axes
            extra['show_unit_cell'] = self.gui.ui.get_widget('/MenuBar/ViewMenu/ShowUnitCell').get_active()
            extra['bbox'] = bbox
            extra['colors'] = self.gui.get_colors(rgb=True)[self.gui.images.visible]
            remove_hidden = True
        if len(indices) == 1:
            # Saving a single configuration is always possible.
            write(filename, self.gui.images.get_atoms(indices[0],
                                                      remove_hidden=remove_hidden),
                  **extra)
        elif self.support_multi[suffix]:
            images = [self.gui.images.get_atoms(i, remove_hidden=remove_hidden) 
                      for i in indices]
            write(filename, images, **extra)
        else:
            # We want to write multiple images, but the file format does not support it.
            # The solution is to write multiple files, inserting a number in the file name
            # before the suffix.
            filename = filename.replace('%', '%%') # Preserve % in filenames.
            suffixpos = filename.rfind('.')
            filename = filename[:suffixpos] + '%05d' + filename[suffixpos:]
            for i, idx in enumerate(indices):
                write(filename % (i,), 
                      self.gui.images.get_atoms(idx, remove_hidden=remove_hidden),
                      **extra)
            oops("Wrote %d files" % (len(indices),),
                 (filename % (0,)) + ' .. ' + (filename % (len(indices)-1,)))
        self.destroy()

    def add_filters(self, chooser):
        # Add file type filters
        self.name_to_suffix = {}
        self.is_graphics = {} 
        self.support_multi = {}
        for name, suffix, graphics, multi in self.filetypes:
            if suffix is None:
                name = _(name)
            else:
                name = '%s (%s)' % (_(name), suffix)
            filt = gtk.FileFilter()
            filt.set_name(name)
            if suffix is None:
                filt.add_pattern('*')
            elif suffix == 'POSCAR':
                filt.add_pattern('*POSCAR*')
            else:
                filt.add_pattern('*.'+suffix)
            self.name_to_suffix[name] = suffix
            if suffix is not None:
                self.is_graphics[suffix] = graphics
                self.support_multi[suffix] = multi
            chooser.add_filter(filt)

        
