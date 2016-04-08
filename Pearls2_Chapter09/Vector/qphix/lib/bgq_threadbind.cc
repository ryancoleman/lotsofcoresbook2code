#include "qphix/threadbind.h"
#include "qphix/print_utils.h"


#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>
#include <spi/include/kernel/location.h>


#include <omp.h>

namespace QPhiX {

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

      pid_t pid = (pid_t) syscall(SYS_gettid);
      // Bind the OMP threads to hardware threads
      if((sched_setaffinity(pid, sizeof(set), &set)) == -1) {
	masterPrintf("WARN: Cound not do sched_setaffinity\n");
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
  
    masterPrintf("ThreadBindings\n");
    for (int i = 0; i < omp_get_max_threads(); ++i)
      masterPrintf("OMP thread %d: core = %d, hardware thread = %d\n",
			  i, cids[i], htids[i]);
  }

}; // END namespace

