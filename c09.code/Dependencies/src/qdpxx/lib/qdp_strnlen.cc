#include "qdp_config.h"

#ifndef HAVE_STRNLEN
#include "qdp_strnlen.h"
#include <stdlib.h>

extern "C" { 

  /* Drop in replacement for strnlen */
  size_t strnlen(const char *s, size_t maxlen)
  {
    size_t pos;

    /* What to do if s is NULL? Not defined in man page */
    if (s == NULL) { 
      return 0;
    }
    else { 

      /* Scan through looking for a terminator */
      for(pos=0; pos < maxlen; pos++) {
	if( s[pos] == '\0' ) {

	  /* The length of the string including the terminator
	     is pos+1, but I am not supposed to count the terminator
	     according to the man page */

	  return pos;
	}
      }
      /* I did not find a terminator so return maxlen */
      return maxlen;
    }
  }
};

#endif
