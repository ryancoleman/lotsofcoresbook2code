/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include "bmgs.h"
#include <pthread.h>
#include "../extensions.h"

#ifdef K
struct IP1DA{
  int thread_id;
  int nthds;
  const T* a;
  int n;
  int m;
  T* b;
  int *skip;
};

void *IP1DW(void *threadarg)
{
  struct IP1DA *args = (struct IP1DA *) threadarg;
  int m = args->m;
  int chunksize = m / args->nthds + 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= m)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > m)
    nend = m;

  for (int j = nstart; j < nend; j++)
    {
      const T* aa = args->a + j * (K - 1 - args->skip[1] + args->n);
      T* bb = args->b + j;
      for (int i = 0; i < args->n; i++)
        {
          if (i == 0 && args->skip[0])
            bb -= m;
          else
            bb[0] = aa[0];

          if (i == args->n - 1 && args->skip[1])
            bb -= m;
          else
            {
              if (K == 2)
                bb[m] = 0.5 * (aa[0] + aa[1]);
              else if (K == 4)
                bb[m] = ( 0.5625 * (aa[ 0] + aa[1]) +
                         -0.0625 * (aa[-1] + aa[2]));
              else if (K == 6)
                bb[m] = ( 0.58593750 * (aa[ 0] + aa[1]) +
                         -0.09765625 * (aa[-1] + aa[2]) +
                          0.01171875 * (aa[-2] + aa[3]));
              else
                bb[m] = ( 0.59814453125 * (aa[ 0] + aa[1]) +
                         -0.11962890625 * (aa[-1] + aa[2]) +
                          0.02392578125 * (aa[-2] + aa[3]) +
                         -0.00244140625 * (aa[-3] + aa[4]));
            }
          aa++;
          bb += 2 * m;
        }
    }
  return NULL;
}

void IP1D(const T* a, int n, int m, T* b, int skip[2])
{
  a += K / 2 - 1;

  int nthds = 1;
#ifdef GPAW_OMP_MONLY
  if (getenv("OMP_NUM_THREADS") != NULL)
    nthds = atoi(getenv("OMP_NUM_THREADS"));
#endif
  struct IP1DA *wargs = GPAW_MALLOC(struct IP1DA, nthds);
  pthread_t *thds = GPAW_MALLOC(pthread_t, nthds);

  for(int i=0; i < nthds; i++)
    {
      (wargs+i)->thread_id = i;
      (wargs+i)->nthds = nthds;
      (wargs+i)->a = a;
      (wargs+i)->n = n;
      (wargs+i)->m = m;
      (wargs+i)->b = b;
      (wargs+i)->skip = skip;
    }
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_create(thds + i, NULL, IP1DW, (void*) (wargs+i));
#endif
  IP1DW(wargs);
#ifdef GPAW_OMP_MONLY
  for(int i=1; i < nthds; i++)
    pthread_join(*(thds+i), NULL);
#endif
  free(wargs);
  free(thds);
}

#else
#  define K 2
#  define IP1D Z(bmgs_interpolate1D2)
#  define IP1DA Z(bmgs_interpolate1D2_args)
#  define IP1DW Z(bmgs_interpolate1D2_worker)
#  include "interpolate.c"
#  undef IP1D
#  undef IP1DA
#  undef IP1DW
#  undef K
#  define K 4
#  define IP1D Z(bmgs_interpolate1D4)
#  define IP1DA Z(bmgs_interpolate1D4_args)
#  define IP1DW Z(bmgs_interpolate1D4_worker)
#  include "interpolate.c"
#  undef IP1D
#  undef IP1DA
#  undef IP1DW
#  undef K
#  define K 6
#  define IP1D Z(bmgs_interpolate1D6)
#  define IP1DA Z(bmgs_interpolate1D6_args)
#  define IP1DW Z(bmgs_interpolate1D6_worker)
#  include "interpolate.c"
#  undef IP1D
#  undef IP1DA
#  undef IP1DW
#  undef K
#  define K 8
#  define IP1D Z(bmgs_interpolate1D8)
#  define IP1DA Z(bmgs_interpolate1D8_args)
#  define IP1DW Z(bmgs_interpolate1D8_worker)
#  include "interpolate.c"
#  undef IP1D
#  undef IP1DA
#  undef IP1DW
#  undef K

void Z(bmgs_interpolate)(int k, int skip[3][2],
       const T* a, const int size[3], T* b, T* w)
{
  void (*ip)(const T*, int, int, T*, int[2]);
  if (k == 2)
    ip = Z(bmgs_interpolate1D2);
  else if (k == 4)
    ip = Z(bmgs_interpolate1D4);
  else if (k == 6)
    ip = Z(bmgs_interpolate1D6);
  else
    ip = Z(bmgs_interpolate1D8);

  int e = k - 1;

  ip(a, size[2] - e + skip[2][1],
     size[0] *
     size[1],
     b, skip[2]);
  ip(b, size[1] - e + skip[1][1],
     size[0] *
     ((size[2] - e) * 2 - skip[2][0] + skip[2][1]),
     w, skip[1]);
  ip(w, size[0] - e + skip[0][1],
     ((size[1] - e) * 2 - skip[1][0] + skip[1][1]) *
     ((size[2] - e) * 2 - skip[2][0] + skip[2][1]),
     b, skip[0]);
}
#endif
