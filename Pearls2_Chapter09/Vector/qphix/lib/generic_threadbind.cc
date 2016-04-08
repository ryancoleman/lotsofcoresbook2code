#include "qphix/threadbind.h"
#include "qphix/print_utils.h"

namespace QPhiX {


  // Default: User should use runtime
  void setThreadAffinity(int nCores, int threadsPerCore) 
  {
    masterPrintf("Direct thread affinity not implemented for this architecture. Please use GOMP_CPU_AFFINITY, KMP_AFFINITY or other appropriate method\n");
  }

  // Default: User should use runtiem
  void reportAffinity()
  {
    masterPrintf("Affinity reporting not implemented for this architecture\n");
  }

}; // END namespace

