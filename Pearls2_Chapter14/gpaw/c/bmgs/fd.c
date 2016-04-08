/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include "bmgs.h"
#include <pthread.h>
#include "../extensions.h"

struct Z(fds){
  int thread_id;
  int nthds;
  const bmgsstencil* s;
  const T* a;
  T* b;
};

void *Z(bmgs_fd_worker)(void *threadarg)
{
  struct Z(fds) *args = (struct Z(fds) *) threadarg;
  const T* a = args->a;
  T* b = args->b;
  const bmgsstencil* s = args->s;

  int chunksize = s->n[0] / args->nthds + 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= s->n[0])
    return NULL;
  int nend = nstart + chunksize;
  if (nend > s->n[0])
    nend = s->n[0];


  for (int i0 = nstart; i0 < nend; i0++)
  {
    const T* aa = a + i0 * (s->j[1] + s->n[1] * (s->j[2] + s->n[2]));
    T* bb = b + i0 * s->n[1] * s->n[2];

    for (int i1 = 0; i1 < s->n[1]; i1++)
      {
        for (int i2 = 0; i2 < s->n[2]; i2++)
          {
            T x = 0.0;
            for (int c = 0; c < s->ncoefs; c++)
              x += aa[s->offsets[c]] * s->coefs[c];
            *bb++ = x;
            aa++;
          }
        aa += s->j[2];
      }
  }
  return NULL;
}



void Z(bmgs_fd)(const bmgsstencil* s, const T* a, T* b)
{
  a += (s->j[0] + s->j[1] + s->j[2]) / 2;

  int nthds = 1;
#ifdef GPAW_OMP_MONLY
  if (getenv("OMP_NUM_THREADS") != NULL)
    nthds = atoi(getenv("OMP_NUM_THREADS"));
#endif
  struct Z(fds) *wargs = GPAW_MALLOC(struct Z(fds), nthds);
  pthread_t *thds = GPAW_MALLOC(pthread_t, nthds);

  for(int i=0; i < nthds; i++)
    {
      (wargs+i)->thread_id = i;
      (wargs+i)->nthds = nthds;
      (wargs+i)->s = s;
      (wargs+i)->a = a;
      (wargs+i)->b = b;
    }
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_create(thds + i, NULL, Z(bmgs_fd_worker), (void*) (wargs+i));
#endif
  Z(bmgs_fd_worker)(wargs);
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_join(*(thds+i), NULL);
#endif
  free(wargs);
  free(thds);

}
