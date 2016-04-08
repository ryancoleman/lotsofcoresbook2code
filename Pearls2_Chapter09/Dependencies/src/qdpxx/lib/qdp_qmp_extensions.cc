#include "qmp.h"
#include "qdp_qmp_extensions.h"
#include <stdlib.h>
#include <stdio.h>

extern "C" { 

#ifndef HAVE_QMP_VERBOSE
void   QMP_verbose(QMP_bool_t verbose) 
{
   fprintf(stderr, "QMP_verbose does nothing for me or you\n");
}
#endif

};
 
