#ifndef LIME_FSEEKO_H
#define LIME_FSEEKO_H

#include "lime_config.h"

#ifndef HAVE_FSEEKO

/* If fseeko, ftell is not defined then define our own versions of it
   which just call fseek and ftell but keep the same interface */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

int fseeko(FILE *stream, off_t offset, int whence);
off_t ftello(FILE *stream);

#else
#include <stdio.h>
#endif

#endif
