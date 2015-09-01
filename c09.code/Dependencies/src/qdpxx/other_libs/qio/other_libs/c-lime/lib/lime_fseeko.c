#include <lime_config.h>
#include <stdio.h>

#ifndef HAVE_FSEEKO

#ifdef HAVE_DCAP
lime_misconfiguration_fseeko_not_declared;
#endif

#include "lime_fseeko.h"

int fseeko(FILE *stream, off_t offset, int whence) {
  return fseek(stream, (long)offset, whence);
}
 
off_t ftello(FILE *stream) {
  return (off_t)ftell(stream);
}

#endif

