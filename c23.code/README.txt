DL_MESO_LBE: README
===================

DL_MESO_LBE is a mesoscale modelling code based on the Lattice Boltzmann 
Equation (LBE) method. It forms part of the DL_MESO general purpose mesoscale 
simulation package, developed by Michael Seaton at STFC Daresbury Laboratory 
for CCP5 (the UK Collaborative Computational Project for computer simulation 
of condensed phases) under a grant provided by EPSRC. 

The extract provided here is used with permission of the copyright holders 
and is limited to the serial version of DL_MESO_LBE (with OpenMP 
multithreading) to run the provided two LBE test cases. The full DL_MESO 
package is supplied to individuals under an academic licence, which is free 
of cost to academic scientists pursuing scientific research of a 
non-commercial nature. Please visit www.ccp5.ac.uk/DL_MESO to register. 
Commercial organisations interested in acquiring the package should approach 
Dr M. A. Seaton at Daresbury Laboratory in the first instance. STFC Daresbury 
Laboratory is the sole centre for distribution of the full package.

Build and run on host
---------------------
change to the source directory
	cd 4

To build an executable for DL_MESO_LBE on the host if running in Linux, type 
‘make’ at the command line to generate the executable ‘slbe.out’. To build an 
executable in Windows called ‘slbe.exe’, type ’nmake /f Makefile.nmake’. By 
default, these makefiles will generate the executable using optimization level 
2 and produce a full optimization report (including details in vectorization). 

Additional options for compilation, VTune analysis, Advisor analysis and 
command line execution are available in the makefiles and can be invoked by 
following the above command with one of the following words (e.g. ‘make asm’):

asm		Generate assembly code in the file slbe.s
dump		Disassemble executable and create text file (slbe.dump.txt) 
		with assembler mnemonics
advanced	Carry out VTune Analysis: examine advanced hotspots
ge		Carry out VTune Analysis: apply general exploration
con		Carry out VTune Analysis: examine concurrency
uops		Carry out VTune Analysis: examine micro-operation flows
port		Carry out VTune Analysis: examine port saturation
mem		Carry out VTune Analysis: identify memory access issues
survey		Carry out Advisor Analysis: collect exploration data
trip		Carry out Advisor Analysis: collect data on trip counts
serial		Run DL_MESO_LBE on the host with a single OpenMP thread
run		Run DL_MESO_LBE on the host using all available threads
gui		Run Amplifier GUI
clean		Delete executable file
nuke		Delete executable file, reports, outputs and collected data

DL_MESO_LBE requires two input files in the same directory as the executable 
to run: lbin.sys for system information and lbin.spa for boundary information. 
The supplied test cases are as follows: 

lbin_4fluid	Models separation of 4 immiscible fluids
lbin_diffuse	Models diffusion of 7 solutes in a fluid with thermal effects

One of these files needs to be renamed as lbin.sys to run with DL_MESO_LBE. 
Both test cases use periodic boundary conditions and use the same (empty) 
lbin.spa file. DL_MESO_LBE can be run either through the Makefile or directly 
at the command line with the command ‘slbe.exe’ or ‘./slbe.out’ (depending on 
operating system).

Build and run on MIC
--------------------

To build an executable for DL_MESO_LBE on MIC, either type ‘make mic’ or 
’nmake /f Makefile.nmake mic’. The makefiles will generate an executable 
(‘slbe.out’ or ‘slbe.exe’) using optimisation level 2 and produce a full 
optimization report.

Additional options for MIC compilation are available in the makefiles with 
‘make’ or ‘nmake' followed by one of the following words:

asm-mic		Generate assembly code for MIC in the file slbe.s
dump		Disassemble executable and create text file (slbe.dump.txt) 
		with assembler mnemonics
clean		Delete executable file
nuke		Delete executable file, reports, outputs and collected data

As for the host, DL_MESO_LBE on the MIC will require lbin.sys and lbin.spa 
input files to run each of the available test cases. The executable should be 
run natively on the MIC: either ssh onto the MIC (‘ssh mic0’) and run with 
’slbe.exe’ or ‘./slbe.out’ (after ensuring LD_LIBRARY_PATH is set correctly to 
point to the directory with libiomp5.so), or launch from the host using mpirun 
or mpiexec to start a one-core job and import all environment variables, e.g.

mpirun -envall -n 1 -host mic0 ./mic

The contents of the folder (including executable and input files) may need to 
be copied over to the MIC prior to running DL_MESO_LBE. If this is required,  
the following command can be used:

scp -r * mic0:

Outputs
-------

DL_MESO_LBE writes out information about the system at regular intervals (every 
100 time steps for the test cases): the total mass of the system, the masses of 
each fluid and the total system momentum. These values can be used to verify 
the code is calculating correctly: the fluid masses for both test cases should 
not change over 300 time steps, while the momentum for the 4 fluid system will 
increase to a constant value as the fluids start separating from each other.

A number of structured grid VTK files at the time steps reported above will be 
created: these can be plotted using e.g. Paraview to visualize fluid densities, 
solute concentrations and temperatures.

Queries and copyright message
-----------------------------

Please send any queries about DL_MESO to michael.seaton@stfc.ac.uk

Copyright (c) 2015 STFC Daresbury Laboratory

Dr Michael Seaton,
Scientific Computing Department,
STFC Daresbury Laboratory,
Warrington,
WA4 4AD,
UNITED KINGDOM
