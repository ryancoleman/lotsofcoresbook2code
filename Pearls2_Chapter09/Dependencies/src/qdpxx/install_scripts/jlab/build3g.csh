#!/bin/tcsh

(cd parscalar-gigE ; do_config ; make ; make)
(cd parscalar-gigE-double ; do_config ; make ; make)
