QMMM
====

:Who:
    Elvar

Have a functioning QM/MM code - which will soon be available in the
development version. It produces reasonable radial distribution
functions btw. classical water and quantum water, and quantum ions,
when compared to experiment and other similar codes.

As the code is written now it assumes that the quantum region is
completely screened by the classical environment (no pbc), and that
there is enough classical stuff (which is cheap) such that the
classical interaction btw. supercells can be treated with the minimum
image convention.

Currently working on:

A) a sensible way to optimize Lennard-Jones
   parameters for the QM system such that QM-MM interactions are as close
   as possible to QM-QM interactions (seems to be the norm in literature)

B) an algorithm to fit classical inter- and intrapotentials to QM
   results (i.e. MM parameters on the run)

C) Quasi-Newton for QM/MM systems.
