/* Concurrent read version */

/* $Id: htime.c 867 2011-10-30 14:44:04Z wrp $ */
/* $Revision: 867 $  */

#include <stdio.h>
#include <time.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/time.h>
#ifdef TIMES
#include <sys/times.h>
#else
#undef TIMES
#endif
#endif

#ifndef HZ
#define HZ 100
#endif

#define GETTIMEOFDAY

long s_time ()			/* returns time in milliseconds */
{

#ifdef GETTIMEOFDAY
struct timeval tv;
gettimeofday(&tv, NULL);
return tv.tv_sec*1000 + tv.tv_usec/1000;
#else
#ifndef TIMES
  time_t time(), tt;
  return time(&tt)*1000;
#else
  struct tms tt;
  times(&tt);
#ifdef CLK_TCK
  return tt.tms_utime*1000/CLK_TCK;
#else
  return tt.tms_utime*1000/HZ;
#endif
#endif
#endif
}

void ptime (FILE *fp, long time)		/* prints the time */
{
  fprintf (fp, "%6.3f",(double)(time)/1000.0);
}

