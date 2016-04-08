// -*- C++ -*-

/*! \file
 * \brief Default memory allocator for QDP
 *
 */

#ifndef QDP_DEFAULT_ALLOCATOR
#define QDP_DEFAULT_ALLOCATOR

#include "qdp_allocator.h"
#include "qdp_stdio.h"
#include "qdp_singleton.h"
#include <string>
#include <map>

namespace QDP
{
  namespace Allocator
  {

    // Specialise allocator to the default case
    class QDPDefaultAllocator {
    private:
      // Disallow Copies
      QDPDefaultAllocator(const QDPDefaultAllocator& c) {}

      // Disallow assignments (copies by another name)
      void operator=(const QDPDefaultAllocator& c) {}

      // Disallow creation / destruction by anyone except 
      // the singleton CreateUsingNew policy which is a "friend"
      // I don't like friends but this follows Alexandrescu's advice
      // on p154 of Modern C++ Design (A. Alexandrescu)
      QDPDefaultAllocator() {init();}
      ~QDPDefaultAllocator() {}

      friend class QDP::CreateUsingNew<QDP::Allocator::QDPDefaultAllocator>;
    public:

      // Pusher
      void pushFunc(const char* func, int line);
  
      // Popper
      void popFunc();
  
      //! Allocator function. Allocates n_bytes, into a memory pool
      //! This is a default implementation, with only 1 memory pool
      //! So we simply ignore the memory pool hint.
      void*
      allocate(size_t n_bytes,const MemoryPoolHint& mem_pool_hint);

      //! Free an aligned pointer, which was allocated by us.
      void 
      free(void *mem);

      //! Dump the map
      void
      dump();

    protected:
      void init();
    };

    // Turn into a Singleton. Create with CreateUsingNew
    // Has NoDestroy lifetime, as it may be needed for 
    // the destruction policy is No Destroy, so the 
    // Singleton is not cleaned up on exit. This is so 
    // that static objects can refer to it with confidence
    // in their own destruction, not having to worry that
    // atexit() may have destroyed the allocator before
    // the static objects need to feed memory. 
    typedef SingletonHolder<QDP::Allocator::QDPDefaultAllocator,
			    QDP::CreateUsingNew,
			    QDP::NoDestroy,
			    QDP::SingleThreaded> theQDPAllocator;

  } // namespace Allocator
} // namespace QDP

#endif
