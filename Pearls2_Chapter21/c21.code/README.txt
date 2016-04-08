Directions to build and run codes for "Prefetch Tuning Optimizations"

####### MIC versions #######

### SHOC MD ###
1. cd MIC
2. extract shoc-mic.tar.gz
2. cd shock-mic
3. Build with "make md"
3. Run with ./RUN_MD.sh

### Smith-Waterman
1. cd MIC extract SWMIC.tar.gz
2. cd SWMIC
3. Build with "make"
3. Run with "./RUN.sh"

### Stream ###
1. cd MIC extract STREAM.tar.gz
2. cd STREAM
3. Build with "make"
4. Run with "./RUN.sh"

###############


####### Xeon versions #######

### SHOC MD ###
1. cd Xeon and extract shoc-xeon.tar.gz
2. cd shoc-xeon
3. Run "make md" to build MD.
4. Execute ./RUN.sh to run the benchmark.

### Smith-Waterman ###
The Xeon implementation of Smith Waterman algorithm is included in the ssearch36 tool from FASTA package(http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml).

1. cd Xeon and extract fasta-36.3.7a.tar.gz and sw_data.tar.gz archives.
2. cd fasta-36.3.7a
2. Build with ./BUILD.sh.
3. Run with ./RUN.sh.


### Stream ###
1. cd Xeon  and extract stream-xeon.tar.gz
2. cd stream-xeon
3. Build with "make"
4. Run with ./RUN.sh

################
