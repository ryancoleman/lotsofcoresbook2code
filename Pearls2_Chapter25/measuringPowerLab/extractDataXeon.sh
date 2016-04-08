#!/bin/sh.exe
# This script extracts data from the log created by the collectData.sh for 
#    Xeon data. 
#
# Note: this only works for Xeon data. The log file for Phi data is different 
#    because the format of the power data collect differs. See the script 
#    extractDataMic.sh for the MIC/KNC coprocessor.

echo "This script parses and extracts only log data obtained from"
echo "    the Intel Xeon processor"

if [ $# -ne 2 ]
  then
    echo
    echo "Usage: $0 <source_log_file> <generated_extracted_csv_data_file>"
    exit -1
fi

# put header entries for data
echo '"Epoch time (s)","energy (J)", "thread number","Mops/s","Mops/s/thread"'> tmp.csv
# parse and extract power and execution data from Xeon log
sed -n "/^>>>/ p; /^ Mop/ p; /^OMP/ p;" $1 \
   | sed 's/Epoch.*$//g;/^>>>time/{N;N;N;N;s/\n/,/g};s/[>\/:_()a-zA-Z= ]\+//g' \
   | sed '/^>>>/{N; s/sec/,/; s/[>time:nrgy\nJ ]*//g; s/$/,,,/}' \
   >> tmp.csv
mv tmp.csv $2    #send extracted data to specified file
#cat tmp.csv     #send extracted data to standard out
echo "$1 data extracted and placed into $2 in csv form"
