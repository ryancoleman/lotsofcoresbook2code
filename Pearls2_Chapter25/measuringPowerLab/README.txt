type make to build pwr-read on the xeon processor

cd ../NPB and follow the READMEs to download, build and run the 
parallel benchmarks
- Readme.NPB.Xeon.sh #to build and run on Xeon
- Readme.NPB.MIC.sh #to build and run on Xeon



executionDialog_rx.docx (or .pdf) contains the actual tutorial/run-through.
That file should be where you start. 


This document contains Taylor's notes and, hopefully helpful, hints.

Consider using nohup as many of your runs may take hours depending upon
working set size (e.g. C or D) you use. Many systems, such as the one 
I'm using, disconnect your ssh session if it has been inactive more than 
a certain number of minutes. 

To parse and convert my data into a csv form, I used a cygwin shell on
my PC.

Remember to set your LD_LIBRARY_PATH correctly so that your benchmark
can find libiomp5.so, otherwise you'll get the annoying message, "cannot 
open shared object file: No such file or directory."

The script only performs 4 experiments. To have statistically significance
assuming a Gaussian distribution of errors, you should perfrom 16 but that's
an awfully long runtime.

The first sample of the experiment takes the longest because there is
only one thread being used. The runtime for each subsequent sample will get 
progressively faster from there.

For the mic card experiments, remember to copy the proper mic libiomp5.so
library to your working directory. It will be under the lib/mic path

[twkidd@knightscorner6]$ ls /opt/intel/composerxe/bin/compilervars.^
[twkidd@knightscorner6]$ find /opt/intel -name libiomp5.so
/opt/intel/composer_xe_2015.0.090/compiler/lib/ia32/libiomp5.so
/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/libiomp5.so
/opt/intel/composer_xe_2015.0.090/compiler/lib/mic/libiomp5.so
[twkidd@knightscorner6]$ 

On the mic card, always remember to set LD_LIBRARY_PATH, 
e.g. "export LD_LIBRARY_PATH=."

If you are using the intel compiler, then executing 
composerxe/bin/compilervars.sh will properly setup your *host* LD_LIBRARY_PATH.

The format and type of data collected for KNC and Xeon are very *different*.
Make sure you use the correct parsing and extraction scripts.


pwr-read-xeon.c will only work for recent generations of Intel x89 
architectures, SNB/Sandy Bridge or newer. Older versions use different
and non-standard PMU event names. If you want to modify pwr-read-xeon.c 
for an older architecture, look at "Intel 64 and IA-32 architectures 
software developer's manual combined volumes 3A, 3B, and 3C: System 
programming guide". Volume 3C will be most relevant. The PMU events 
for the different architectures are in the tables. E.g. Table 13-11
lists the NHM events are listed.
