#ifndef H_LOCALBARRIER
#define H_LOCALBARRIER

typedef struct corebarrier_t
{
    volatile int userbarrier_arrive;
    volatile int userbarrier_depart;
    int padding[14];    // pad things out to fill the cacheline.
} corebarrier_t;

typedef struct barrier_t
{
    int usersense;
    int mycoretid;
    int coreval;
    corebarrier_t* me;
} barrier_t;

__declspec(target(mic)) void userCoreBarrier(barrier_t *bar);
#endif
