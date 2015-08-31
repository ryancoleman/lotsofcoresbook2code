Just include "openvec.h" file in your C/C++ code and you are ready to go.

To compile code examples.
make arch={intel,neon,knc,gccavx}
If no arch is provided, gcc with SSE will be used.

To run a self tests:
make run_selftest

The stencil* samples will allocate 3.9GB,
if you want allocate less memory edit stencil_common.h
and change problem size.
Blocking definitions are located on stencil_common.h too.



Specific information for KNC:

For compiling:
make arch=knc

And for running the tests (if your KNC's name is mic0):
make arch=knc knc_target=mic0 run_selftest


