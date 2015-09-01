#!/bin/tcsh

(cd scalar; do_config ; make ; make install)
(cd scalar-double ; do_config ; make ; make install)

(cd parscalar-gigE ; do_config ; make ; make install)
(cd parscalar-gigE-double ; do_config ; make ; make install)
(cd parscalar-single ; do_config ; make ; make install)
