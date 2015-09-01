//
/*! 
 * @file
 * @brief Utility routines like info and error reporting
 */

#include <cstdarg>
#include <unistd.h>

#include "qdp.h"


// QCDOC Hack. QCDOC doesn't have gethostname(char *hostname, size_t size)
// Must provide a suitable stub
#ifndef HAVE_GETHOSTNAME

#include "qdp_strnlen.h"

int gethostname(char *hostname, size_t size)
{
   char *my_host_name="node";

   // My Host Size returns either the length of my_host_name
   // or the size of the supplied array - 1 char
   // (this means that if size -1 is passed it wastes 1 char
   size_t my_host_size=strnlen(my_host_name, size-1);
  
   strncpy(hostname, my_host_name, my_host_size);

   // Null terminate. This should be ok, because my_host_size is at most
   // size-1;
   hostname[my_host_size]='\0'; 
   return 0;
} 
#endif

namespace QDP {

/**
 * Simple information display routine
 */
int
QDP_info (const char* format, ...)
{
  va_list argp;
  char    info[128], hostname[256];
  char    buffer[1024];

  /* get machine host name */
  gethostname (hostname, sizeof (hostname) - 1);
  int rank = Layout::nodeNumber();
  int size = Layout::numNodes();
  sprintf (info, "QDP m%d,n%d@%s info: ", 
	   rank, size, hostname);

  va_start (argp, format);
  int status = vsprintf (buffer, format, argp);
  va_end (argp);

  fprintf (stdout, "%s %s\n", info, buffer);
  return status;
}


/**
 * Simple error display routine
 */
int
QDP_error (const char* format, ...)
{
  va_list argp;
  char    info[128], hostname[256];
  char    buffer[1024];

  /* get machine host name */
  gethostname (hostname, sizeof (hostname) - 1);
  int rank = Layout::nodeNumber();
  int size = Layout::numNodes();
  sprintf (info, "QDP m%d,n%d@%s error: ", rank, size, hostname);
	   
  va_start (argp, format);
  int status = vsprintf (buffer, format, argp);
  va_end (argp);

  fprintf (stderr, "%s %s\n", info, buffer);
  return status;
}

/**
 * Simple error display and abort routine
 */
void
QDP_error_exit (const char* format, ...)
{
  va_list argp;
  char    info[128], hostname[256];
  char    buffer[1024];

  /* get machine host name */
  gethostname (hostname, sizeof (hostname) - 1);
  int rank = Layout::nodeNumber();
  int size = Layout::numNodes();
  sprintf (info, "QDP m%d,n%d@%s fatal error: ", rank, size, hostname);
	   
  va_start (argp, format);
  vsprintf (buffer, format, argp);
  va_end (argp);

  fprintf (stderr, "%s %s\n", info, buffer);

  QDP_abort(1);  // Abort knows how to shutdown the machine
}

/**
 * Enable or Disable verbose mode.
 */
void
QDP_verbose (bool verbose)
{
  // Currently a NOP
}


//-----------------------------------------------------------------------------
//! Unique-ify a list
/*! Given a list of ints, collapse it to a list of the unique ints */
multi1d<int>
uniquify_list(const multi1d<int>& ll)
{
  multi1d<int> d(ll.size());

  // Enter the first element as unique to prime the search
  int ipos = 0;
  int num = 0;
  int prev_node;
  
  d[num++] = prev_node = ll[ipos++];

  // Find the unique source nodes
  while (ipos < ll.size())
  {
    int this_node = ll[ipos++];

    if (this_node != prev_node)
    {
      // Has this node occured before?
      bool found = false;
      for(int i=0; i < num; ++i)
	if (d[i] == this_node)
	{
	  found = true;
	  break;
	}

      // If this is the first time this value has occurred, enter it
      if (! found)
	d[num++] = this_node;
    }

    prev_node = this_node;
  }

  // Copy into a compact size array
  multi1d<int> dd(num);
  for(int i=0; i < num; ++i)
    dd[i] = d[i];

  return dd;
}


} // namespace QDP;
