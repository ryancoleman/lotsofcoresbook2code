#ifndef QPHIX_THREADBIND_H
#define QPHIX_THREADBIND_H

namespace QPhiX
{
  void setThreadAffinity(int nCores, int threadsPerCore);
  void reportAffinity(void);
}




#endif
