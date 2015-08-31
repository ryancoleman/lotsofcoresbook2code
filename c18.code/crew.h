/*
    crew.h -- header file for Intel Crew runtime interfaces.

    Copyright 2012-2013 Intel Corporation.  All Rights Reserved.

    The source code contained or described herein and all documents related
    to the source code ("Material") are owned by Intel Corporation or its
    suppliers or licensors.  Title to the Material remains with Intel
    Corporation or its suppliers and licensors.  The Material is protected
    by worldwide copyright laws and treaty provisions.  No part of the
    Material may be used, copied, reproduced, modified, published, uploaded,
    posted, transmitted, distributed, or disclosed in any way without
    Intel's prior express written permission.

    No license under any patent, copyright, trade secret or other
    intellectual property right is granted to or conferred upon you by
    disclosure or delivery of the Materials, either expressly, by
    implication, inducement, estoppel or otherwise.  Any license under such
    intellectual property rights must be express and approved by Intel in
    writing.
*/
#ifndef __CREW_H
# define __CREW_H

# ifdef _INTEL_CREW
# error Preprocessor symbol _INTEL_CREW defined.  Please define __INTEL_CREW instead.
# endif

# if (__INTEL_CREW)&&defined(__MIC__)

/* Crew is enabled, declare functions and define macros. */
# ifdef __cplusplus
extern "C" {
# endif

#if defined(__INTEL_OFFLOAD)
#pragma offload_attribute(push,target(mic))
#endif

#if __INTEL_COMPILER_BUILD_DATE >= 20130301
   /* Revision uses kmp_ instead of kmpc_ as prefix. */
   extern void kmp_crew_create();
   extern void kmp_crew_destroy();
   extern int kmp_crew_get_max_size();
#else
   /* Backwards compatibility for icc 2013u1 and 2023u2. 
      To be removed when icc 2013u3 is ubiquitous. */ 
   extern void kmpc_crew_create();
   extern void kmpc_crew_destroy();
   extern int kmpc_crew_get_max_size();
#  define kmp_crew_create kmpc_crew_create
#  define kmp_crew_destroy kmpc_crew_destroy
#  define kmp_crew_get_max_size kmpc_crew_get_max_size
#endif

#if defined(__INTEL_OFFLOAD)
#pragma offload_attribute(pop)
#endif

# ifdef __cplusplus
}
# endif
# define CREW_FOR_LOOP _Pragma("intel_crew parallel for")

# else

/* Crew disabled. Define macros to remove Crew functionality. */
#  define kmp_crew_create()  ((void)0)
#  define kmp_crew_destroy() ((void)0)
#  define kmp_crew_get_max_size() (1)
#  define CREW_FOR_LOOP

#endif  /* Crew enabled */

#endif /* __CREW_H */
