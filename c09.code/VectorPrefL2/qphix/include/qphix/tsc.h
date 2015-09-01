#ifndef QPHIX_TSC_H_
#define QPHIX_TSC_H_

// Stuff for rdtsc
#include <sys/types.h>
#include <sys/sysctl.h>

typedef long long TSC_tick;
#define CLOCK_NOW(a)     (a) = __rdtsc()
#endif
