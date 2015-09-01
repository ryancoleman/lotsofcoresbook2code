#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

#ifndef CLOCK_SGI_CYCLE
#define CLOCK_SGI_CYCLE CLOCK_PROCESS_CPUTIME_ID
#endif

#define seci(t) ((double)t.tv_sec+(double)t.tv_nsec*1.e-9)

double csecond(void )
{
         struct timespec tt;
         double tot;
         int ret = clock_gettime(CLOCK_SGI_CYCLE, &tt);
         if(ret<0) perror("clock_gettime");
         tot = seci(tt);
         return tot;
}
