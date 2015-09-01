#include "qdp.h"
#include "qdp_threadbind.h"

namespace QDP {

  // Default: User should use runtime
  void setThreadAffinity(int nCores, int threadsPerCore) 
  {
    QDPIO::cout << "Direct thread affinity not implemented for this architecture. Please use GOMP_CPU_AFFINITY, KMP_AFFINITY or other appropriate method\n";
  }

  // Default: User should use runtiem
  void reportAffinity()
  {
    QDPIO::cout << "Affinity reporting not implemented for this architecture\n";
  }


}; // END namespace

