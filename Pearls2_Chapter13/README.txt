To build and run type the following command and follow the instructions:
sh pearls-example.sh 


Launching the "pearls-example.sh" script:
  - Checks out the specific revision of the Uintah Computational Framework (UCF) utilized for
    the experiments discussed within this chapter
  - Creates the appropriate directory structure
  - Sources Intel Parallel Studio XE 2015 (Cluster Edition)
  - Downloads, extracts, and cross-compiles libxml2 and zlib for the coprocessor
  - Cross compiles the UCF for the coprocessor
  - Transfers required Intel and UCF files to mic0

NOTE 1: To anonymously checkout the UCF, use the username "anonymous" with a blank password.
NOTE 2: This scripts assumes source scripts and library paths for Intel Parallel Studio XE 2015 
        (Cluster Edition). Please update these per your specific build environment.
NOTE 3: It is required that multi-threaded versions of Intel MPI are used with this example.
NOTE 4: This script assumes that you have a local account with a matching
        $HOME directory on mic0.

To launch your first Uintah Computational Framework simulations:
    ssh mic0
    export LD_LIBRARY_PATH=$HOME/intel/lib64:$LD_LIBRARY_PATH
    export PATH=$HOME/intel/bin:$PATH
    cd tmp.example
    ./run-pearls-example.sh

Launching the "run-pearls-example.sh" script:
  - Executes a RMCRT simulation with the fastest run configuration identified during these
    experiments for 2 threads per physical core
  - Executes a RMCRT simulation with the fastest run configuration identified during these
    experiments for 3 threads per physical core
  - Executes a RMCRT simulation with the fastest run configuration identified during these
    experiments for 4 threads per physical core

NOTE: A list of environment variables (e.g. those used to specify thread affinity) can be found
      within the "/uintah/branches/pearls/src/environmentalFlags.txt" file.

____________________
For More Information

To learn more about the Uintah Computational Framework, please refer to the link below:
  http://www.uintah.utah.edu/

To explore Uintah and C-SAFE-related publications, please refer to the link below:
  http://www.uintah.utah.edu/pubs/pubs.html

To learn more about radiative heat transfer and the discrete ordinates method-based radiation model, please refer to the publication linked below:
  Spinti, J.P., Thornock, J.N., Eddings, E.G., Smith, P.J., and Sarofim, A.F. ``Heat transfer
  to objects in pool fires.'' In Transport Phenomena in Fires, WIT Press, Southampton,
  UK (2008).

To learn more about the reverse Monte-Carlo ray tracing-based radiation model, please refer to the publication linked below:
  Humphrey, A., Meng, Q., Berzins, M., and Harman, T. ``Radiation modeling using the Uintah
  heterogeneous CPU/GPU runtime system.'' In Proceedings of the 1st Conference of the Extreme
  Science and Engineering Discovery Environment (XSEDE12). ACM, 2012.

To download the latest version of Uintah, please refer to the link below:
  http://www.sci.utah.edu/uintah-survey.html

To learn more about installing and using the Uintah Computational Framework, please refer to the link below:
  http://uintah-build.sci.utah.edu/trac/wiki
