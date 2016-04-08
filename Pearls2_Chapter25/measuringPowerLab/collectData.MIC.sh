
#This script is very simple intentionally. Making a bullet 
#    proof script would have obscured the key statements

#change to reflect the benchmark you are using (e.g. ep.D.x)
#    and the naming of the log files (e.g. epC_${sys_name}_...)

affinity=scatter
#affinity=compact

#for KNC
sys_name=knc
num_expts=4 #for KNC
num_threads=244 #for KNC
num_threads_core=4
benchmark=./ep.B.knc
log_prescript=epB

#for Xeon
#sys_name=xeon
#num_expts=4
#num_threads=32
#num_threads_core=2
#benchmark=./ep.C.x
#log_prescript=epD

#sample the nth thread per core
#  Xeon (hyperthreading) n = {1|2}
#  KNC (4 threads/core) n = {1|2|3|4}
sample_thread=1

echo >&2 "Experiment Parameters"
echo >&2 "   affinity($affinity)"
echo >&2 "   number of experiments($num_expts)"
echo >&2 "   total number of threads per experiment($num_threads)"
echo >&2 "   number of threads per core($num_threads_core)"
echo >&2 "   sampling thread number $sample_thread of $num_threads_core"
echo

#script to collect the data
export KMP_AFFINITY="$affinity,granularity=fine"
  #sample for experiments
  for j in $(seq 1 $num_expts)
    do echo
    >&2 echo
    >&2 echo EXPT $j
    #collect for experiment
    for i in $(seq $sample_thread $num_threads_core $num_threads)
      do export OMP_NUM_THREADS=$i
      echo OMP_NUM_THREADS\($OMP_NUM_THREADS\) KMP_AFFINITY\($KMP_AFFINITY\)
      >&2 echo -n $OMP_NUM_THREADS
      #execution of NPB EP benchmark class C (may use class D for Xeon)
      #./ep.C.x
      ${benchmark}
      echo PWR********
      #collect power data; may need to be sudo
      echo cat /sys/class/micras/power #for Xeon Phi
      #sudo ./pwr-read             #for Xeon
    #done > epC_${sys_name}_${sample_thread}_${num_threads}_${num_threads_core}_${affinity}_$j.log
    done > ${log_prescript}_${sys_name}_${sample_thread}_${num_threads}_${num_threads_core}_${affinity}_$j.log
  done

