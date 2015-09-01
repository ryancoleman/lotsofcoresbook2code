#ifndef QPHIX_DSLASH_UTILS_H
#define QPHIX_DSLASH_UTILS_H



#include <cerrno>
#include "qphix/qphix_config.h"
#include "qphix/tsc.h"


#ifndef QPHIX_USE_MM_MALLOC
#include <cstdlib>
#else 

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
#include <immintrin.h>
#endif

#endif

using namespace std;

// Options
#undef SCATTERED_THREAD_IDS

#define S0     (0)
#define S1     (1)
#define S2     (2)
#define S3     (3)

#define C0     (0)
#define C1     (1)
#define C2     (2)

#define RE     (0)
#define IM     (1)

// Use _mm_malloc on non MIC for now.


namespace QPhiX { 

  struct BlockPhase {
    int by;
    int bz;
    int bt;
    int nt;
    int group_tid;
    int cid_t;
    int cid_yz;

  };


 inline
  void *aligned_malloc(int size, int alignment)
    {
      void *v;
      int ok;
      ok = posix_memalign((void **)&v, alignment, size);
      if( ok != 0) {
        return NULL;
      }
      return v;
    }




};

#ifdef __MIC__

// MIC BUFFERS are HUGE PAGES
#include <fcntl.h>
#include <sys/mman.h>
#include "qphix/print_utils.h"

#define HUGE_PAGE_SIZE (2*1024*1024)
inline
void* BUFFER_MALLOC(size_t size, int alignment) 
{
  // Round size up to nearest huge_page size
  int n_pages = size/HUGE_PAGE_SIZE;
  if ( size % HUGE_PAGE_SIZE != 0 ) { 
    n_pages++;
  }
  
  void *ret_val = mmap(NULL,size,PROT_READ|PROT_WRITE, 
		       MAP_ANONYMOUS | MAP_SHARED | MAP_HUGETLB | MAP_POPULATE, -1, 0);

  // OK if MMAP Fails it returns -1 rather than NULL, so trap on that and return NULL if needed
  if( (long int)ret_val == -1 ) return NULL;
  
  return ret_val;
}

inline
void BUFFER_FREE(void *addr,size_t length)
{
  int n_pages = length/HUGE_PAGE_SIZE;
  if ( length % HUGE_PAGE_SIZE != 0 ) { 
    n_pages++;
  }
  
  int status = munmap(addr,n_pages*HUGE_PAGE_SIZE);
  switch (status) {
  case 0:
    break;
  case EINVAL :
    QPhiX::localPrintf("munmap returned EINVAL\n");
    abort();
    break;
  default:
    break;
  };
}

inline 
void* ALIGNED_MALLOC(size_t size, unsigned int alignment) {
  void* ret_val =  _mm_malloc(size, alignment);
  return ret_val;
}

inline void
ALIGNED_FREE(void *addr)
{
  _mm_free(addr);
}

#else 
// NON MIC Architectures
// We can implement large pages if needed. Otherwise use user alignment
inline
void* BUFFER_MALLOC(size_t size, int alignment) 
{

#ifdef QPHIX_USE_MM_MALLOC
  void* ret_val =  _mm_malloc(size, alignment);
#else
  void* ret_val =  QPhiX::aligned_malloc(size, alignment);
#endif
  return ret_val;
}

inline
void BUFFER_FREE(void *addr,size_t length)
{
#ifdef QPHIX_USE_MM_MALLOC
  _mm_free(addr);
#else
  free(addr);
#endif
}

inline 
void* ALIGNED_MALLOC(size_t size, unsigned int alignment) 
{
#ifdef QPHIX_USE_MM_MALLOC
  void* ret_val =  _mm_malloc(size, alignment);
#else
  void* ret_val =  QPhiX::aligned_malloc(size, alignment);
#endif
  return ret_val;
}

inline void
ALIGNED_FREE(void *addr)
{
#ifdef QPHIX_USE_MM_MALLOC
  _mm_free(addr);
#else
  free(addr);
#endif
}

#endif // __MIC__

#if defined(QPHIX_MIC_SOURCE)
#warning including barrier
#include "qphix/Barrier_mic.h"
#else
#include "qphix/Barrier_stubs.h"
#endif

#define MIN(a,b)   ( (a) < (b) ? (a) : (b) )

#define BARRIER_TSLICES 16
#define N_PROBES 8


#endif
