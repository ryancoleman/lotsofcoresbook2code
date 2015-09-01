#!/bin/tcsh

(cd scalar; do_config ; make ; make)
(cd parscalar-gm ; do_config ; make ; make)
(cd scalar-double ; do_config ; make ; make)
(cd parscalar-gm-double ; do_config ; make ; make)
