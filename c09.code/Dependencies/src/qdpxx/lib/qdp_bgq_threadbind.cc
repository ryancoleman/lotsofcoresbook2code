#include "qdp.h"
#include "qdp_threadbind.h"


#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>
#include <spi/include/kernel/location.h>

#ifdef QDP_USE_OMP_THREADS
#include <omp.h>
#endif

namespace QDP 
{

#ifdef QDP_USE_OMP_THREADS
  void setThreadAffinity(int nCores, int threadsPerCore)
  {
    #pragma omp parallel 
    {

      // Get the OpenMP thread number
      int tid = omp_get_thread_num();

      // Split into core and SIMT ID (assuming SIMT runs fastest) 
      int core  = tid/threadsPerCore;
      int smtid = tid - core*threadsPerCore;

      // Convert to hardware processor ID, basically using same scheme as 
      // Table 3-2 of the IBM redbook: 
      // http://www.redbooks.ibm.com/redbooks/pdfs/sg247948.pdf            

      // NB: 4 is hardwired here for BG/Q.  (NB: This would let us run with 
      // 'gaps' i.e. even 2 threads per core but still get the binding right.)       

      int hw_procid = smtid + 4*core;   
        
      cpu_set_t set;

      CPU_ZERO(&set);
      CPU_SET(hw_procid, &set);

      // Bind the OMP threads to hardware threads
      if((sched_setaffinity(0, sizeof(set), &set)) == -1) {
	QDPIO::cout << "WARN: Cound not do sched_setaffinity for thread " << tid << "\n";
      }
    } // pragma omp parallel
  } // function



  // BlueGene/Q reporting via SPI Locations interface
  void reportAffinity()
  {

    uint32_t cids[64], htids[64];

#pragma omp parallel
    {
      htids[omp_get_thread_num()] = Kernel_ProcessorThreadID();
      cids[omp_get_thread_num()] = Kernel_ProcessorCoreID();
    }
  
    QDPIO::cout << "Thread Bindings\n";
    for (int i = 0; i < omp_get_max_threads(); ++i)
      QDPIO::cout <<"OMP thread " << i << ": core = " << cids[i] << " hardware thread = " << htids[i] << std::endl;
  }

#else
  void setThreadAffinity(int nCores, int threadsPerCore)
  {
 	QDPIO::cout << "Thread binding only implemented for OpenMP Threads\n";
  }

  void reportAffinity()
  {
        QDPIO::cout << "reportAffinity: only implemented for OpenMP Threads \n";
  }
#endif 
  }; // END namespace
