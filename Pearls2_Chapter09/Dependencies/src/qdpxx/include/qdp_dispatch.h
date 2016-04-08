#ifndef QDP_DISPATCH_H
#define QDP_DISPATCH_H


#include "qdp_config.h"

#if defined(QDP_USE_OMP_THREADS)
#warning QDP using OpenMP threading
#include <omp.h>
namespace QDP {    
/* OpenMP threading version of the dispatch.*/

  
  inline
    int qdpNumThreads()
  {
    return omp_get_max_threads();
  }


template<class Arg>
void dispatch_to_threads(int numSiteTable, Arg a, void (*func)(int,int,int, Arg*)){
   
  int threads_num;
  int myId;
  int low = 0;
  int high = numSiteTable;
   
#pragma omp parallel shared(numSiteTable, threads_num, a) private(myId, low, high) default(shared)
    {
     
      threads_num = omp_get_num_threads();
      myId = omp_get_thread_num();
      low = (numSiteTable*myId)/threads_num;

      // NB: high can never be too high since for the last thread
      // myId + 1 = (threads_num - 1) + 1 = threads_num.
      // So numSiteTable*(myId+1) will always be a strict multiple of threads_num
      // and so truncation issues will not bite.
      // I am addig in the parentheses tho to force the precedence
      high =(numSiteTable*(myId+1))/threads_num;
 
      func(low, high, myId, &a);
    }
}
}

#else 

#if defined(QDP_USE_QMT_THREADS)
#warning QDP using QMT threading

 /* QMT threading version of the dispatch. Call the qmt_call routine
     with userfunc, numSiteTable, and argument */

#include <qmt.h>

namespace QDP { 

 inline
   int qdpNumThreads()
 {
   return qmt_num_threads();
 }


template<class Arg>
void dispatch_to_threads(int numSiteTable, Arg a, void (*func)(int,int,int,Arg*)){
 
   qmt_call((qmt_userfunc_t)func, numSiteTable, &a);
 
}
}
#else
namespace QDP {

 inline
 int qdpNumThreads()
 {
   return 1;
 }

template<class Arg>
void dispatch_to_threads(int numSiteTable, Arg a, void (*func)(int,int,int,Arg*)){
   
  int low = 0;
  int high = numSiteTable;
   
  func(low, high, 0, &a);
    
 }

} 
#endif
#endif

#endif
