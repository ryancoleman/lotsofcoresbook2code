#ifndef SSE_PREFETCH_H
#define SSE_PREFETCH_H

#include "qdp_config.h"

/* Read Prefetch a 64 byte (8 doubles) line */
#define PREFETCHNTA(a) _mm_prefetch((a), _MM_HINT_NTA)

#define PREFETCH(a) _mm_prefetch((a), _MM_HINT_T0)

/* Write prefetch a 64 byte (8 doubles line) */
#if defined(QDP_USE_3DNOW)
#warning "Using 3DNOW Write Prefetch"
#include <mm3dnow.h>
#define PREFETCHW(a) _m_prefetchw((a))

#else

#define PREFETCHW(a) _mm_prefetch((a), _MM_HINT_T0)

#endif

#endif
