/* FILE LICENSE TAG: SAMPLE */

//=======================timer stuff ===============
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#else
#include <sys/time.h>
#endif

// Timing Prototypes: Include where needed
#include "dtime.h"
//void  dtimeInit();
//double dtimeGet();
//double dtimeSince(double t1,char *str);
//double dtimeElapse(double t1);

// Timing Definitions: Included once.
static double Global_Tickrate = 0.0;
#define DTIMEINITMSG _T("You need to call dtimeInit(); once on Windows")
#define DTIMEINITERR _T("Coding Error")

double dtimeGet()
{
    double t;
#ifdef _WIN32
    LARGE_INTEGER tvalue;
    QueryPerformanceCounter(&tvalue);
    t = (double)(tvalue.QuadPart) / Global_Tickrate;
#else
    struct timeval tv1;
    gettimeofday(&tv1, NULL);
    t = (tv1.tv_sec) + 1.0e-6 * tv1.tv_usec;
#endif
    return (t);
}

double dtimeSince(double t1, char *str)
{
    double t2;
    double telapsed;
#ifdef _WIN32
    LARGE_INTEGER tvalue;
    QueryPerformanceCounter(&tvalue);
    t2 = (double)(tvalue.QuadPart);
    if (Global_Tickrate > 0.0) {
        telapsed = (t2 - t1) / Global_Tickrate;
    } else {
        telapsed = -1.0;
        MessageBox(NULL, DTIMEINITMSG, DTIMEINITERR, MB_OK);
    }
#else
    struct timeval tv2;
    gettimeofday(&tv2, NULL);
    t2 = (tv2.tv_sec) + 1.0e-6 * tv2.tv_usec;
    telapsed = t2 - t1;
#endif
    printf("%.5g secs <-- Elapsed Time for: '%s'\r\n", telapsed, str);
    fflush(stdout);
    return (telapsed);
}

double dtimeElapse(double t1)
{
    double t2;
    double telapsed;
#ifdef _WIN32
    LARGE_INTEGER tvalue;
    QueryPerformanceCounter(&tvalue);
    t2 = (double)(tvalue.QuadPart);
    if (Global_Tickrate > 0.0) {
        telapsed = (t2 - t1) / Global_Tickrate;
    } else {
        telapsed = -1.0;
        MessageBox(NULL, DTIMEINITMSG, DTIMEINITERR, MB_OK);
    }
#else
    struct timeval tv2;
    gettimeofday(&tv2, NULL);
    t2 = (tv2.tv_sec) + 1.0e-6 * tv2.tv_usec;
    telapsed = t2 - t1;
#endif
    return (telapsed);
}


void dtimeInit()
{
#ifdef _WIN32
    LARGE_INTEGER cpufreq;
    QueryPerformanceFrequency(&cpufreq);
    Global_Tickrate = (double)(cpufreq.QuadPart);
#endif
    return;
}
//===================================== end of timer stuff
