#ifndef __DCAP_OVERLOAD_H__
#define __DCAP_OVERLOAD_H__

#ifdef HAVE_DCAP
  #include <stdio.h>
  #include <unistd.h>
  #include <dcap.h>
  #define DCAP(token)  dc_ ## token
  #if defined(_FILE_OFFSET_BITS) && _FILE_OFFSET_BITS == 64
     #define DCAPL(token) dc_ ## token ## 64
  #else
     #define DCAPL(token) dc_ ## token
  #endif
#else
  #define DCAP(token)  token
  #define DCAPL(token) token
#endif

#endif
