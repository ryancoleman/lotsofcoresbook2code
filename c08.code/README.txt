 Source files
  original source file
  	am_call_original.cpp

  Optimized source file using #pragma simd
       am_call.cpp

  Optimized source file using C/C++ vector extension
  	am_call_vector.cpp
  
 Executive files
  Original Version - Scalar code runs in one thread
       am_call_original

  Scalar Serial - Scalar, serial  Optimized code
       am_call_ss

  Vector Serial - Scalar Serial Version with #pragma SIMD vectorization
       am_call_vs

  Scalar Parallel - Scalar Serial code with OpenMP parallelization
       am_call_sp

  Vector Parallel - Vector Serial code with OpenMP parallelization
       am_call_vp

  Vector Extension Vectorization
  	am_call_vector Vector version runs on AVX2 Haswell only

  KNC Vector Parallel -- vector parallel version for KNC target
       am_call_knc

To run am_call_knc on system call myhost as device 0, here is what you have to do or call 'runknc.sh'
scp am_call_knc myhost-mic0:
scp /opt/intel/composer_xe_2015.1.133/tbb/lib/mic/libtbbmalloc.so.2 myhost-mic0:
scp /opt/intel/composer_xe_2015.1.133/tbb/lib/mic/libtbbmalloc.so myhost-mic0:
scp /opt/intel/composer_xe_2015.1.133/compiler/lib/mic/libiomp5.so myhost-mic0:
ssh myhost-mic0 "export LD_LIBRARY_PATH=. ; export OMP_NUM_THREADS=244 ; export KMP_AFFINITY='compact,granularity=fine' ; ./am_call_knc"
