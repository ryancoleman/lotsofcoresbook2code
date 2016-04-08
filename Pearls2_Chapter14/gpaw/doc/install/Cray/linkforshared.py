    def link_shared_object (self,
                         objects,
                         output_filename,
                         output_dir=None,
                         libraries=None,
                         library_dirs=None,
                         runtime_library_dirs=None,
                         export_symbols=None,
                         debug=0,
                         extra_preargs=None,
                         extra_postargs=None,
                         build_temp=None,
                         target_lang=None):

        if output_dir is None:
            (output_dir, output_filename) = os.path.split(output_filename)
        output_fullname = os.path.join(output_dir, output_filename)
        output_fullname = os.path.abspath(output_fullname)
        linkline = "%s %s" % (output_filename[:-2], output_fullname)
        for l in library_dirs:
            linkline += " -L" + l
        for l in libraries:
            linkline += " -l" + l
        old_fmt = self.static_lib_format
        self.static_lib_format = "%s%.0s"
        self.create_static_lib(objects,
                               output_filename,
                               output_dir,
                               debug,
                               target_lang)

        self.static_lib_format = old_fmt
        print "Append to Setup: ", linkline
