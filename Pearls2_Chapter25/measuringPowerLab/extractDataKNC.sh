#!/bin/sh.exe
# This script extracts data from the log created by the 
#collectData.sh for KNC data. 
#
# Note: this only works for KNC data. The log file for Xeon data is different 
#    because the format of the power data collect differs. See the script 
#    extractDataXeon.sh for the Xeon processor.

echo "This script parses and extracts only log data obtained from"
echo "    the Intel Xeon Phi coprocessor"

if [ $# -ne 2 ]
  then
    echo
    echo "Usage: $0 <source_log_file> <generated_extracted_csv_data_file>"
    exit -1
fi

# put header entries for data
echo '"(blank)","CPU time","Threads Total","Threads Avail","Mops/s","Mops/s/thread","Win1 Pwr (uW)","Win2 Pwr (uW)"' > tmp.csv
# parse and extract power and execution data from Xeon log
sed -n '/^OMP/ p; /Completed/{N;N;N;N;N;N; p;}; /Mop/ p; /^PWR/ {N;N;p}' $1 \
    | sed -n '/Time/{s/[^0-9.]*/\n/; p}; /threads/{s/[^0-9.]*//g; p}; /Mop/{s/[^0-9.]*//g; p}; /^PWR/{N;N;s/^PWR[*]*\n//;p}' \
    | sed -n '/^$/{N;N;N;N;N;N;N; s/\n/,/g; p}' \
    >> tmp.csv;
mv tmp.csv $2    #send extracted data to specified file
echo
echo "$1 data extracted and placed into $2 in csv form"
