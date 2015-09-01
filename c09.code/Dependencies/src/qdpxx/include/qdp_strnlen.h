#ifndef QDP_STRNLEN_H
#define QDP_STRNLEN_H

#include "qdp_config.h"

#include <stdlib.h>
#include <sys/types.h>
#include <string.h>

// Declare this. If it is already declared in string.h it should
// still be ok.
// On qcdoc it appears to have linkage but no prototype in which 
// case it should also be declared
// Likewise if a prototype were not to exist and no linkage were
// to exist declaring it is still OK as qdp_strnlen.cc will provide
// the linkage
extern "C" { 
  size_t strnlen(const char *s, size_t maxlen);
};

#endif
