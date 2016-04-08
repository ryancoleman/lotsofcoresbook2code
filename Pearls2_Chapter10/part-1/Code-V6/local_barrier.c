#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <immintrin.h>
#include "local_barrier.h"

__declspec(target(mic)) void cpu_pause()
{
    #ifdef __MIC__
    _mm_delay_64(50);
    #else
    _mm_pause();
    #endif
}

// local barrier bar is a thread-local variable.
__declspec(target(mic)) void userCoreBarrier(barrier_t *bar)
{
    int mysense = bar->usersense;
    int coretid = bar->mycoretid;
    int mycoresense = mysense ? bar->coreval : 0;
    corebarrier_t *me = bar->me;

    // signal my arrival
    ((char *)&me->userbarrier_arrive)[coretid] = mysense;
    // wait for others to arrive
    while (me->userbarrier_arrive != mycoresense)
        cpu_pause();

    // signal my departure
    ((char *)&me->userbarrier_depart)[coretid] = mysense;
    // wait for others to depart
    while (me->userbarrier_depart != mycoresense)
        cpu_pause();

    bar->usersense = 1 - mysense;
}
