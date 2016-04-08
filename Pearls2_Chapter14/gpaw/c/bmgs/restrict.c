/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include "bmgs.h"
#include <pthread.h>
#include "../extensions.h"

#ifdef K
struct RST1DA{
  int thread_id;
  int nthds;
  const T* a;
  int n;
  int m;
  T* b;
};

void *RST1DW(void *threadarg)
{
  struct RST1DA *args = (struct RST1DA *) threadarg;
  int m = args->m;
  int chunksize = m / args->nthds + 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= m)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > m)
    nend = m;

  for (int j = 0; j < m; j++)
    {
      const T* aa = args->a + j * (args->n * 2 + K * 2 - 3);
      T* bb = args->b + j;

      for (int i = 0; i < args->n; i++)
        {
          if      (K == 2)
            bb[0] = 0.5 * (aa[0] +
              0.5 * (aa[1] + aa[-1]));

          else if (K == 4)
            bb[0] = 0.5 * (aa[0] +
               0.5625 * (aa[1] + aa[-1]) +
              -0.0625 * (aa[3] + aa[-3]));

          else if (K == 6)
            bb[0] = 0.5 * (aa[0] +
               0.58593750 * (aa[1] + aa[-1]) +
              -0.09765625 * (aa[3] + aa[-3]) +
               0.01171875 * (aa[5] + aa[-5]));

          else
            bb[0] = 0.5 * (aa[0] +
               0.59814453125 * (aa[1] + aa[-1]) +
              -0.11962890625 * (aa[3] + aa[-3]) +
               0.02392578125 * (aa[5] + aa[-5]) +
              -0.00244140625 * (aa[7] + aa[-7]));
          aa += 2;
          bb += m;
        }
    }
  return NULL;
}

void RST1D(const T* a, int n, int m, T* b)
{
  a += K - 1;

  int nthds = 1;
#ifdef GPAW_OMP_MONLY
  if (getenv("OMP_NUM_THREADS") != NULL)
    nthds = atoi(getenv("OMP_NUM_THREADS"));
#endif
  struct RST1DA *wargs = GPAW_MALLOC(struct RST1DA, nthds);
  pthread_t *thds = GPAW_MALLOC(pthread_t, nthds);

  for(int i=0; i < nthds; i++)
    {
      (wargs+i)->thread_id = i;
      (wargs+i)->nthds = nthds;
      (wargs+i)->a = a;
      (wargs+i)->n = n;
      (wargs+i)->m = m;
      (wargs+i)->b = b;
    }
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_create(thds + i, NULL, RST1DW, (void*) (wargs+i));
#endif
  RST1DW(wargs);
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_join(*(thds+i), NULL);
#endif
  free(wargs);
  free(thds);
}

#else
#  define K 2
#  define RST1D Z(bmgs_restrict1D2)
#  define RST1DA Z(bmgs_restrict1D2_args)
#  define RST1DW Z(bmgs_restrict1D2_worker)
#  include "restrict.c"
#  undef RST1D
#  undef RST1DA
#  undef RST1DW
#  undef K
#  define K 4
#  define RST1D Z(bmgs_restrict1D4)
#  define RST1DA Z(bmgs_restrict1D4_args)
#  define RST1DW Z(bmgs_restrict1D4_worker)
#  include "restrict.c"
#  undef RST1D
#  undef RST1DA
#  undef RST1DW
#  undef K
#  define K 6
#  define RST1D Z(bmgs_restrict1D6)
#  define RST1DA Z(bmgs_restrict1D6_args)
#  define RST1DW Z(bmgs_restrict1D6_worker)
#  include "restrict.c"
#  undef RST1D
#  undef RST1DA
#  undef RST1DW
#  undef K
#  define K 8
#  define RST1D Z(bmgs_restrict1D8)
#  define RST1DA Z(bmgs_restrict1D8_args)
#  define RST1DW Z(bmgs_restrict1D8_worker)
#  include "restrict.c"
#  undef RST1D
#  undef RST1DA
#  undef RST1DW
#  undef K

void Z(bmgs_restrict)(int k, T* a, const int n[3], T* b, T* w)
{
  void (*plg)(const T*, int, int, T*);

  if (k == 2)
    plg = Z(bmgs_restrict1D2);
  else if (k == 4)
    plg = Z(bmgs_restrict1D4);
  else if (k == 6)
    plg = Z(bmgs_restrict1D6);
  else
    plg = Z(bmgs_restrict1D8);

  int e = k * 2 - 3;
  plg(a, (n[2] - e) / 2, n[0] * n[1], w);
  plg(w, (n[1] - e) / 2, n[0] * (n[2] - e) / 2, a);
  plg(a, (n[0] - e) / 2, (n[1] - e) * (n[2] - e) / 4, b);
}

#endif
