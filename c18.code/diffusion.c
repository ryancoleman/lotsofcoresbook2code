#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>

#define __USE_MISC
#include <sys/mman.h>
#include "crew.h"      /* Header file for using Crew */

#define REAL float

#define CACHE_LINE_SIZE 64
#define N_REALS_PER_CACHE_LINE (CACHE_LINE_SIZE / sizeof(REAL))

// default NX
#define NX (256)
//#define NX (512)

// default count
#define COUNT 1200

#define NXP_DELTA 16

//#define NP_MINOR 4

#ifndef M_PI
#define M_PI (3.1415926535897932384626)
#endif

#if defined(__MIC__)
#define WAIT_A_BIT _mm_delay_32(10)
#else
#define WAIT_A_BIT _mm_pause();
#endif

int np_major;
int np_minor;
int majorIndex[256];
int minorIndex[256];
int nCores = 0;
int nHTs = 0;

/* divide the Y or Z for load balance */
void divideBlocks(int ntot, int npart, int* adist)
{
  int bat=ntot/npart;
  int residue = npart-ntot%npart;
  adist[0] = 0;
  for(int i=0; i<npart; i++)
  {
    if(i<residue)
      adist[i+1] = adist[i] + bat;
    else
      adist[i+1] = adist[i] + bat+1;
  }
}

void init(REAL *buff, const int nx, const int ny, const int nz,
          const REAL kx, const REAL ky, const REAL kz,
          const REAL dx, const REAL dy, const REAL dz,
          const REAL kappa, const REAL time) {
  REAL ax, ay, az;
  //int jz, jy, jx;
  ax = exp(-kappa*time*(kx*kx));
  ay = exp(-kappa*time*(ky*ky));
  az = exp(-kappa*time*(kz*kz));

//#pragma omp parallel for collapse(2)
#pragma omp parallel for num_threads(60)
  for (int jy = 0; jy < ny; jy++) {
    for (int jz = 0; jz < nz; jz++) {
      for (int jx = 0; jx < nx; jx++) {
        int j = jz*(nx + NXP_DELTA)*ny + jy*(nx + NXP_DELTA) + jx;
        REAL x = dx*((REAL)(jx + 0.5));
        REAL y = dy*((REAL)(jy + 0.5));
        REAL z = dz*((REAL)(jz + 0.5));
        REAL f0 = (REAL)0.125
          *(1.0 - ax*cos(kx*x))
          *(1.0 - ay*cos(ky*y))
          *(1.0 - az*cos(kz*z));
        buff[j] = f0;
      }
    }
  }

}

void inityz(REAL *buff, const int nx, const int ny, const int nz,
          const REAL kx, const REAL ky, const REAL kz,
          const REAL dx, const REAL dy, const REAL dz,
          const REAL kappa, const REAL time) {
#pragma omp parallel 
  {
    REAL ax = exp(-kappa*time*(kx*kx));
    REAL ay = exp(-kappa*time*(ky*ky));
    REAL az = exp(-kappa*time*(kz*kz));
    const int ip=omp_get_thread_num()/np_minor;
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];
    const int ipy=omp_get_thread_num()%np_minor;
    const int y0=minorIndex[ipy];
    const int y1=minorIndex[ipy+1];

      for (int y = y0; y < y1; y++) {
        for (int z = z0; z < z1; z++) {
          int c =  0 + y * (nx + NXP_DELTA) + z * (nx + NXP_DELTA) * ny;
          for(int x=0; x<nx; x++,c++)
          {
            REAL vx = dx*((REAL)(x + 0.5));
            REAL vy = dy*((REAL)(y + 0.5));
            REAL vz = dz*((REAL)(z + 0.5));
            REAL f0 = (REAL)0.125
              *(1.0 - ax*cos(kx*vx))
              *(1.0 - ay*cos(ky*vy))
              *(1.0 - az*cos(kz*vz));
            buff[c] = f0;
          }
        }
      }
  }
}


REAL accuracy(const REAL *b1, REAL *b2, const int len) {
  REAL err = 0.0;
  int i;
  for (i = 0; i < len; i++) {
    err += (b1[i] - b2[i]) * (b1[i] - b2[i]);
  }
  return (REAL)sqrt(err/len);
}

/*
Write the entire output array to out
*/

void write_out(REAL *f1, int nx, int ny, int nz, FILE* out) {
  if (out) {
    if(NXP_DELTA) {
      // When NXP_DELTA non-zero, collect data into larger buffer
      // This will reduce the number of writes by a factor of ny
      REAL buff[nx*ny] __attribute__ ((aligned(CACHE_LINE_SIZE)));
      size_t nitems = nx * ny;
      int src = 0;
      for (int z = 0; z < nz; z++) {
        int dst = 0;
        for (int y = 0; y < ny; y++) {
#pragma vector nontemporal
#pragma simd
          for (int x = 0; x < nx; x++) {
            buff[dst + x] = f1[src + x];
          } // 
          src += nx + NXP_DELTA;
          dst += nx;
        } // for (int y = 0; y < ny; y++)
        fwrite(buff, sizeof(REAL), nitems, out);
      } // for (int z = 0; z < nz; z++)
    } else {
      // When NXP_DELTA zero, write the buffer in place
      size_t nitems = nx * ny * nz;
      fwrite(f1, sizeof(REAL), nitems, out);
    }
  } // if (out)
} // void write_out(REAL *f1, int nx, int ny, int nz, FILE* out)

/*
Write an X/Y plane of output array to out
*/

volatile int write_xy_out_next_sequence = 0;

void write_xy_out(REAL *f1, int nx, int ny, int nz, int sequence, FILE* out) {
  if (out) {
    int z = sequence % nz;
    if(NXP_DELTA) {
      // When NXP_DELTA non-zero, collect data into larger buffer
      // This will reduce the number of writes by a factor of ny
      REAL buff[nx*ny] __attribute__ ((aligned(CACHE_LINE_SIZE)));
      size_t nitems = nx * ny;
      int src = (nx + NXP_DELTA) * ny * z;
      int dst = 0;
      for (int y = 0; y < ny; y++) {
#pragma vector nontemporal
#pragma simd
        for (int x = 0; x < nx; x++) {
          buff[dst + x] = f1[src + x];
        } // 
        src += nx + NXP_DELTA;
        dst += nx;
      } // for (int y = 0; y < ny; y++)
      while(sequence != write_xy_out_next_sequence)
        WAIT_A_BIT;
      fwrite(buff, sizeof(REAL), nitems, out);
      write_xy_out_next_sequence++;
    } else {
      // When NXP_DELTA zero, write the buffer in place
      size_t nitems = nx * ny;
      while(sequence != write_xy_out_next_sequence)
        WAIT_A_BIT;
      fwrite(f1+(nx * ny * z), sizeof(REAL), nitems, out);
      write_xy_out_next_sequence++;
    }
  } // if (out)
} // void write_xy_out(REAL *f1, int nx, int ny, int nz, int sequence, FILE* out)

typedef void (*diffusion_loop_t)(REAL *f1, REAL *f2, int nx, int ny, int nz,
                                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                                 REAL cb, REAL cc, REAL dt,
                                 int count);

static void
diffusion_baseline(REAL *f1, REAL *f2, int nx, int ny, int nz,
                   REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                   REAL cb, REAL cc, REAL dt,
                   int count) {
  int i;
  for (i = 0; i < count; ++i) {
    int z;
    for (z = 0; z < nz; z++) {
      int y;
      for (y = 0; y < ny; y++) {
        int x;
        for (x = 0; x < nx; x++) {
          int c, w, e, n, s, b, t;
          c =  x + y * (nx + NXP_DELTA) + z * (nx + NXP_DELTA) * ny;
          w = (x == 0)    ? c : c - 1;
          e = (x == nx-1) ? c : c + 1;
          n = (y == 0)    ? c : c - (nx + NXP_DELTA);
          s = (y == ny-1) ? c : c + (nx + NXP_DELTA);
          b = (z == 0)    ? c : c - (nx + NXP_DELTA) * ny;
          t = (z == nz-1) ? c : c + (nx + NXP_DELTA) * ny;
          f2[c] = cc * f1[c] + cw * f1[w] + ce * f1[e]
              + cs * f1[s] + cn * f1[n] + cb * f1[b] + ct * f1[t];
        }
      }
    }
    REAL *t = f1;
    f1 = f2;
    f2 = t;
  }
  return;
}

static void
diffusion_openmp(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {
#pragma omp parallel
  {
    REAL *f1_t = f1;
    REAL *f2_t = f2;

    for (int i = 0; i < count; ++i) {
#pragma omp for collapse(3)
      for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
          for (int x = 0; x < nx; x++) {
            int c, w, e, n, s, b, t;
            c =  x + y * (nx + NXP_DELTA) + z * (nx + NXP_DELTA) * ny;
            w = (x == 0)    ? c : c - 1;
            e = (x == nx-1) ? c : c + 1;
            n = (y == 0)    ? c : c - (nx + NXP_DELTA);
            s = (y == ny-1) ? c : c + (nx + NXP_DELTA);
            b = (z == 0)    ? c : c - (nx + NXP_DELTA) * ny;
            t = (z == nz-1) ? c : c + (nx + NXP_DELTA) * ny;
            f2_t[c] = cc * f1_t[c] + cw * f1_t[w] + ce * f1_t[e]
                + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
          }
        }
      }
      REAL *t = f1_t;
      f1_t = f2_t;
      f2_t = t;
    }
  }

  return;
}


static void
diffusion_mic(REAL *restrict f1, REAL *restrict f2, int nx, int ny, int nz,
              REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
              REAL cb, REAL cc, REAL dt, int count) {

  unsigned long (*pmc1)[2], (*pmc2)[2];
  unsigned long pmcs[2];
  unsigned long tsc;
  int nthreads;
  tsc = _rdtsc();
#pragma omp parallel
  {
    REAL *f1_t = f1;
    REAL *f2_t = f2;
    int mythread;

#if defined(PMU)
#pragma omp master
    {
    nthreads = omp_get_num_threads();
#if defined(PMU)
    pmc1 = malloc(nthreads * sizeof(pmc1[0]));
    pmc2 = malloc(nthreads * sizeof(pmc1[0]));
#endif

#if 0
    printf("%d threads running\n", nthreads);
#endif
    }
#pragma omp barrier

    mythread = omp_get_thread_num();
    pmc1[mythread][0] = _rdpmc(0);
    pmc1[mythread][1] = _rdpmc(1);
#endif
    for (int i = 0; i < count; ++i) {
#define YBF 14
#pragma omp for collapse(2)
      for (int yy = 0; yy < ny; yy += YBF) {
      for (int z = 0; z < nz; z++) {
        int ymax = yy + YBF;
        if (ymax >= ny) ymax = ny;
        for (int y = yy; y < ymax; y++) {
          int x;
          int c, w, e, n, s, b, t;
          x = 0;
          c =  x + y * (nx + NXP_DELTA) + z * (nx + NXP_DELTA) * ny;
          w = (x == 0)    ? c : c - 1;
          e = (x == nx-1) ? c : c + 1;
          n = (y == 0)    ? c : c - (nx + NXP_DELTA);
          s = (y == ny-1) ? c : c + (nx + NXP_DELTA);
          b = (z == 0)    ? c : c - (nx + NXP_DELTA) * ny;
          t = (z == nz-1) ? c : c + (nx + NXP_DELTA) * ny;
          f2_t[c] = cc * f1_t[c] + cw * f1_t[w] + ce * f1_t[e]
              + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
#pragma vector nontemporal
#pragma simd
//#pragma ivdep
          for (x = 1; x < nx-1; x++) {
            ++c;
            ++n;
            ++s;
            ++b;
            ++t;
            f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c+1]
                + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
          }
          ++c;
          ++n;
          ++s;
          ++b;
          ++t;
          f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c]
              + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
        } // tile ny
      } // tile nz
      } // block ny
      REAL *t = f1_t;
      f1_t = f2_t;
      f2_t = t;
    } // count
#if defined(PMU)
    pmc2[mythread][0] = _rdpmc(0);
    pmc2[mythread][1] = _rdpmc(1);
#endif
  } // parallel
  tsc = _rdtsc() - tsc;
#if 0
  printf("%lu cycles %g secs\n", tsc, ((double)tsc/1200e6));
#endif

#if defined(PMU)
  pmcs[0] = 0; pmcs[1] = 0;
  for (int j = 0; j < nthreads; ++j) {
    for (int i = 0; i < 2; ++i)
        pmcs[i] += pmc2[j][i] - pmc1[j][i];
  }
  printf("%lu, %lu\n", pmcs[0], pmcs[1]);
#endif

  return;
}

/* The loop over x in diffusion_mic is encapsulated in an inline function.
 */
static inline void diffusion_x_loop(const REAL* f1_t, REAL* f2_t,
    int nx, int ny, int nz, int y, int z,
    REAL cc, REAL ce, REAL cw, REAL cn, REAL cs, REAL ct, REAL cb)
{
  int c, w, e, n, s, b, t;
  int x = 0;
  c =  x + y * (nx + NXP_DELTA) + z * (nx + NXP_DELTA) * ny;
  w = (x == 0)    ? c : c - 1;
  e = (x == nx-1) ? c : c + 1;
  n = (y == 0)    ? c : c - (nx + NXP_DELTA);
  s = (y == ny-1) ? c : c + (nx + NXP_DELTA);
  b = (z == 0)    ? c : c - (nx + NXP_DELTA) * ny;
  t = (z == nz-1) ? c : c + (nx + NXP_DELTA) * ny;
  f2_t[c] = cc * f1_t[c] + cw * f1_t[w] + ce * f1_t[e]
    + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
#pragma vector nontemporal
#pragma simd
  for (x = 1; x < nx-1; x++) {
    ++c;
    ++n;
    ++s;
    ++b;
    ++t;
    f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c+1]
      + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
  }
  ++c;
  ++n;
  ++s;
  ++b;
  ++t;
  f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c]
    + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
}


/** naive nested loops
 * z = over cores
 * y = over threads
 */
static void
diffusion_nested(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {
  for (int i = 0; i < count; ++i) 
  {
    REAL *f1_t = (i&1)?f2:f1;
    REAL *f2_t = (i&1)?f1:f2;
#pragma omp parallel 
    {
#pragma omp for
      for (int z =0; z < nz; z++) {
#pragma omp parallel for 
        for(int y=0; y<ny; y++){
          diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
        }
      }
    }
  }
  return;
} 

/** nested loops, partition is the same as diffusion_nested * z = over cores
 * y = over threads
 * Similar to diffusion_omp2d
 */
static void
diffusion_nested_lb(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {

#pragma omp parallel 
  for (int i = 0; i < count; ++i) 
  {
    REAL *f1_t = (i&1)?f2:f1;
    REAL *f2_t = (i&1)?f1:f2;
    const int ip=omp_get_thread_num();
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];

#pragma omp parallel for firstprivate(z0,z1)
    for(int y=0; y<ny; y++){
      for (int z = z0; z < z1; z++) {
        diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
      }
    }
#pragma omp barrier
  }

  return;
}

static void
diffusion_nested_task(REAL *f1, REAL *f2, int nx, int ny, int nz,
    REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
    REAL cb, REAL cc, REAL dt, int count) {

#pragma omp parallel 
  for (int i = 0; i < count; ++i) 
  {
    REAL *f1_t = (i&1)?f2:f1;
    REAL *f2_t = (i&1)?f1:f2;
    const int ip=omp_get_thread_num();
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];
#pragma omp parallel
    {
#pragma omp for
      for(int yb=0; yb<np_minor; yb++){
#pragma omp task default(shared) firstprivate(z0,z1,yb)
        {
          for (int z = z0; z < z1; z++) {
            for(int y=minorIndex[yb];y<minorIndex[yb+1]; ++y) {
              diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
            }
          }
        }//end-of-task
      }//end-of-omp-for
#pragma omp taskwait
    }//parallel to do task
#pragma omp barrier
  }
  return;
}

static void
diffusion_task(REAL *f1, REAL *f2, int nx, int ny, int nz,
    REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
    REAL cb, REAL cc, REAL dt, int count) {

#pragma omp parallel 
  for (int i = 0; i < count; ++i) {
    REAL *f1_t = (i&1)?f2:f1;
    REAL *f2_t = (i&1)?f1:f2;
    const int nthreads=np_major*np_minor;
#pragma omp for
    for(int t=0; t<nthreads; ++t)
    {
      const int ip=t/np_minor;
      const int ipy=t%np_minor;
      const int z0=majorIndex[ip];
      const int z1=majorIndex[ip+1];
      const int y0=minorIndex[ipy];
      const int y1=minorIndex[ipy+1];
#pragma omp task default(shared) firstprivate(z0,z1,y0,y1)
      {
        for(int y=y0;y<y1; ++y) {
          for (int z = z0; z < z1; z++) {
            diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
          }
        }
      }//end-of-task
    } //endof-omp-for
#pragma omp taskwait
  }

  return;
}


#if defined(__MIC__)
static void
diffusion_crew(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {
#pragma omp parallel
  {
    const int ip=omp_get_thread_num();
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];

    for (int i = 0; i < count; ++i) {
      const REAL *restrict f1_t = (i&1)?f2:f1;
      REAL *restrict f2_t = (i&1)?f1:f2;
#pragma intel_crew parallel for
      for(int y=0; y<ny; y++){
        for (int z = z0; z < z1; z++) {
          diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
        }
      }

#pragma omp barrier

    }
  }
  return;
}
#endif

static void
diffusion_omp2d_yz(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {
#pragma omp parallel 
  {
    const int ip=omp_get_thread_num()/np_minor;
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];
    const int ipy=omp_get_thread_num()%np_minor;
    const int y0=minorIndex[ipy];
    const int y1=minorIndex[ipy+1];

    for (int i = 0; i < count; ++i) {
      const REAL *restrict f1_t = (i&1)?f2:f1;
      REAL *restrict f2_t = (i&1)?f1:f2;
      //reproducing diffusion_mic with collapse(2) 
      for (int y = y0; y < y1; y++) {
        for (int z = z0; z < z1; z++) {
          diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
        }
      }

#pragma omp barrier

    }
  }

  return;
}

static void
diffusion_omp2d_zy(REAL *f1, REAL *f2, int nx, int ny, int nz,
                 REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                 REAL cb, REAL cc, REAL dt, int count) {
#pragma omp parallel 
  {
    const int ip=omp_get_thread_num()/np_minor;
    const int z0=majorIndex[ip];
    const int z1=majorIndex[ip+1];
    const int ipy=omp_get_thread_num()%np_minor;
    const int y0=minorIndex[ipy];
    const int y1=minorIndex[ipy+1];

    for (int i = 0; i < count; ++i) {
      const REAL *restrict f1_t = (i&1)?f2:f1;
      REAL *restrict f2_t = (i&1)?f1:f2;
      for (int z = y0; z < y1; z++) {
        for (int y = z0; y < z1; y++) {
          diffusion_x_loop(f1_t,f2_t,nx,ny,nz,y,z,cc,ce,cw,cn,cs,ct,cb);
        }
      }

#pragma omp barrier

    }
  }

  return;
}

static double cur_second(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}



static void dump_result(REAL *f, int nx, int ny, int nz, char *out_path) {
  FILE *out = fopen(out_path, "w");
  assert(out);
  size_t nitems = nx * ny * nz;
  fwrite(f, sizeof(REAL), nitems, out);
  fclose(out);
}

void help()
{
  printf("diffusion_pad2 <method> [nx=n] [count=n] [nf=n] [np=n] [out=xxx]\n\n");
  printf("Where:\n\n");
  printf("  <method> is required and is one of:\n");
  printf("    base = single threaded base line\n");
  printf("    openmp = simplified threaded using OpenMP\n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("    base, openmp, mic, nested, nested_lb, twoyz, task, crew, \n\n");
  printf("Optional args:\n\n");
  printf("nx=n    The x, y, and z dimensions, default(%d)\n", NX);
  printf("count=n The number of iterations (per frame), default(%d)\n", COUNT);
  printf("nf=n    The number of frames, default(1)\n");
  printf("np=n    Optional for twoyz to specify np_minor, default(%d)\n",nHTs);
  printf("out=xxx Specifies optional output file, default is <method>.out\n");
  printf("        Specifying out= with no file name indicates no output\n");
}

int main(int argc, char *argv[]) 
{
  if(argc == 1)
  {
    help();
    return 0;
  }
int nthreads_1 = omp_get_max_threads();
  int nthreads_2 = 1;
  int nthreads = nthreads_1;
#pragma omp parallel
  {
#pragma omp master
    nthreads_2=omp_get_max_threads();
  }

#if defined(__MIC__)
  nHTs = 4;
#else
  nHTs = 2;
#endif

  if(nthreads_1 != nthreads_2) //nested openmp
  {
    np_major = nthreads_1;
    np_minor = nthreads_2;
    nthreads=np_major*np_minor;
  }
  else
  {
    np_minor=nHTs;
    np_major = nthreads / np_minor;
  }

  struct timeval time_begin, time_end;

  int    nx    = NX;
  int    ny    = NX;
  int    nz    = NX;
  int    nxp   = nx + NXP_DELTA;

  int use_crew=0;
  int omp2d=-1;
  int    count = COUNT;

 

  diffusion_loop_t diffusion_loop = NULL;
  int n_frames = 1;
  char* p_out = argv[1];	// default is <method>.out

  for(int i=2; i<argc; ++i)
  {
    if (strncmp(argv[i], "nx=", 3) == 0) {
      nx=ny=nz=atoi(argv[i]+3);
      nxp = nx + NXP_DELTA;
      if((nx <= 0) || (nx%16)) {
        printf("nx must be multiple of 16\n");
        return -1;
      }
      continue;
    }
    if (strncmp(argv[i], "count=", 6) == 0) {
      count=atoi(argv[i]+6);
      if(count <= 0) {
        printf("count must be > 0\n");
        return -1;
      }
      continue;
    }
    if (strncmp(argv[i], "nf=", 3) == 0) {
      n_frames=atoi(argv[i]+3);
      if(n_frames <= 0) {
        printf("nf must be > 0\n");
        return -1;
      }
      continue;
    }
    if (strncmp(argv[i], "np=", 3) == 0) {
      np_minor=atoi(argv[i]+3);
      np_major=nthreads/np_minor;
      if(np_minor < 1) {
        printf("np must be > 0\n");
        return -1;
      }
      continue;
    }
    if (strncmp(argv[i], "out=", 3) == 0) {
      p_out=argv[i]+4;
      if(p_out[0] == 0)
        p_out = NULL;
      continue;
    }
  }

  if((np_major>256) || (np_minor>256))
  {
    printf("Major %d minor %d\n", np_major, np_minor);
    return -1;
  }

  printf("Num of threads: Total %d Major %d minor %d\n", nthreads,np_major, np_minor);

  divideBlocks(nx,np_major,majorIndex);
  divideBlocks(nx,np_minor,minorIndex);

  if (strcmp(argv[1], "base") == 0) {
    diffusion_loop = diffusion_baseline;
  }
  else
  if (strcmp(argv[1], "openmp") == 0) {
    diffusion_loop = diffusion_openmp;
  }
  else
  if (strcmp(argv[1], "mic") == 0) {
    printf("MIC\n");
    diffusion_loop = diffusion_mic;
  }
  else
  if (strcmp(argv[1], "nested_lb") == 0) {
    printf("Load-balanced nested OMP\n");
    diffusion_loop = diffusion_nested_lb;
  }
  else
  if (strcmp(argv[1], "twoyz") == 0) {
    printf("2D decomposition YZ\n");
    diffusion_loop = diffusion_omp2d_yz;
  }
  else
  if (strcmp(argv[1], "task") == 0) {
    printf("Use task in the outer loop\n");
    diffusion_loop = diffusion_task;
  }
  else
  if (strcmp(argv[1], "task_inner") == 0) {
    printf("Use task in the inner loop\n");
    diffusion_loop = diffusion_nested_task;
  }
#if defined(__MIC__)
  else
  if (strcmp(argv[1], "crew") == 0) {
    printf("Use crew in the inner loop\n");
    use_crew=1;
    kmp_crew_create();
    diffusion_loop = diffusion_crew;
  }
#endif
  else
  {
    help();
    return -1;
  }

  //np_major=nthreads/np_minor;


#if defined(LARGE_PAGES)
  REAL *f1_padded = (REAL *)mmap(0, sizeof(REAL)*nx*nx*nxp + N_REALS_PER_CACHE_LINE*2 + 2*1024*1024,
                            PROT_READ|PROT_WRITE, MAP_HUGETLB|MAP_ANON|MAP_PRIVATE, -1, 0);
  REAL *f2_padded = (REAL *)mmap(0, sizeof(REAL)*nx*nx*nxp + N_REALS_PER_CACHE_LINE*2 + 2*1024*1024,
                            PROT_READ|PROT_WRITE, MAP_HUGETLB|MAP_ANON|MAP_PRIVATE, -1, 0);
#else
  // align the allocations to cache line
  // increase allocation size by 2 cache lines
  REAL *f1_padded = (REAL *)_mm_malloc(
    sizeof(REAL)*(nxp*ny*nz + N_REALS_PER_CACHE_LINE*2),
    CACHE_LINE_SIZE);

  // assure allocation succeeded
  assert(f1_padded != NULL);

  // align the allocations to cache line
  // increase allocation size by 2 cache lines
  REAL *f2_padded = (REAL *)_mm_malloc(
    sizeof(REAL)*(nxp*ny*nz + N_REALS_PER_CACHE_LINE*2),
    CACHE_LINE_SIZE);

  // assure allocation succeeded
  assert(f2_padded != NULL);
#endif
  
  // advance one cache line into buffer
  REAL *f1 = f1_padded + N_REALS_PER_CACHE_LINE;
  
  f1[-1] = 0.0;       // assure cell prior to array not Signaling NaN
  f1[nx*ny*nz] = 0.0; // assure cell following array not Signaling NaN
  
  // advance one cache line into buffer
  REAL *f2 = f2_padded + N_REALS_PER_CACHE_LINE;
  
  f2[-1] = 0.0;       // assure cell prior to array not Signaling NaN
  f2[nx*ny*nz] = 0.0; // assure cell following array not Signaling NaN

  REAL *answer = (REAL *)_mm_malloc(sizeof(REAL) * nxp*ny*nz, CACHE_LINE_SIZE);
  assert(answer != NULL);

  REAL *f_final = NULL;

  REAL   time  = 0.0;

  REAL l, dx, dy, dz, kx, ky, kz, kappa, dt;
  REAL ce, cw, cn, cs, ct, cb, cc;

  l = 1.0;
  kappa = 0.1;
  dx = dy = dz = l / nx;
  kx = ky = kz = 2.0 * M_PI;
  dt = 0.1*dx*dx / kappa;
  int count_total = 0;


  ce = cw = kappa*dt/(dx*dx);
  cn = cs = kappa*dt/(dy*dy);
  ct = cb = kappa*dt/(dz*dz);
  cc = 1.0 - (ce + cw + cn + cs + ct + cb);

  FILE *out = NULL;
  if(p_out) {
    out = fopen(p_out, "w");
    assert(out);
  }

  //printf("Using 2D partition %d %d\n",np_major, np_minor);
  printf("Running diffusion kernel %d times with nx=%d\n", count, nx);
  printf("for %d number of frames\n", n_frames);
  printf("output %s\n", (out==NULL) ? "(no output)" : p_out);

  init(f1, nx, ny, nz, kx, ky, kz, dx, dy, dz, kappa, time);

  gettimeofday(&time_begin, NULL);
  for (int i_frame = 0; i_frame < n_frames; ++i_frame) {
    diffusion_loop(f1, f2, nx, ny, nz, ce, cw, cn, cs, ct, cb, cc, dt, count);
    count_total += count;
    f_final = (count_total % 2)? f2 : f1;
    if (out) write_out(f_final, nx, ny, nz, out);
  } // for (int i_frame = 0; i_frame < n_frames; ++i_frame)
  gettimeofday(&time_end, NULL);

  time = count_total * dt;
  if (out) fclose(out);
  init(answer, nx, ny, nz, kx, ky, kz, dx, dy, dz, kappa, time);
  REAL err = accuracy(f_final, answer, nx*ny*nz);
  double elapsed_time = (time_end.tv_sec - time_begin.tv_sec)
      + (time_end.tv_usec - time_begin.tv_usec)*1.0e-6;
  REAL mflops = (nx*ny*nz)*13.0*count_total/elapsed_time * 1.0e-06;
  double thput = (nx * ny * nz) * sizeof(REAL) * 3.0 * count_total
      / elapsed_time * 1.0e-09;

  fprintf(stderr, "Elapsed time : %.3f (s)\n", elapsed_time);
  fprintf(stderr, "FLOPS        : %.3f (MFlops)\n", mflops);
  fprintf(stderr, "Throughput   : %.3f (GB/s)\n", thput);  
  fprintf(stderr, "Accuracy     : %e\n", err);
  
#if defined(__MIC__)
  if(use_crew) kmp_crew_destroy();
#endif
  //free(f1);
  //free(f2);
  return 0;
}
