/**
 * Create Test Strings in random size
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#define MIN_LEN 5
#define FIRST_CHAR (char)'0'
#define LAST_CHAR  (char)'~'

int
main (int argc, char** argv)
{
  struct timeval tv;
  int num, maxksize, maxdsize, i, tmp, len, k, fixed_size;
  int dsize, ksize;
  char *kstr, *dstr;

  fixed_size = 0;
  if (argc < 4) {
    fprintf (stderr, "Usage: %s numstring largestkey largestdata [fixed size 0/1]\n", argv[0]);
    return 1;
  }
  
  num = atoi  (argv[1]);
  maxksize = atoi (argv[2]);
  maxdsize = atoi (argv[3]);

  if (argc >= 4) 
    fixed_size = atoi (argv[4]);
  if (fixed_size != 0 && fixed_size != 1) {
    fprintf (stderr, "Fixed size must be either 0 or 1.\n");
    return -1;
  }
  
  gettimeofday (&tv, 0);
  srand (tv.tv_sec);

  if (!fixed_size) {
    ksize = maxksize + 1 + 32;
    dsize = maxdsize + 1 + 32;

    kstr = (char *)malloc(maxksize + 1 + 32);
    dstr = (char *)malloc(maxdsize + 1 + 32);
  }
  else {
    ksize = maxksize;
    dsize = maxdsize;

    kstr = (char *)malloc(maxksize);
    dstr = (char *)malloc(maxdsize);
  }

  for (i = 0; i < num; i++) {
    /* create key */
    memset (kstr, 0, ksize);
    if (!fixed_size)
      len = rand() % maxksize;
    else
      len = ksize - 1;
    if (len < MIN_LEN)
      len = MIN_LEN;
    

    for (k = 0; k < len; k++) {
      tmp = rand() % LAST_CHAR;
      if (tmp < FIRST_CHAR)
	tmp += FIRST_CHAR;
      kstr[k] = (char)tmp;
    }
    if (!fixed_size)
      printf ("%s-%d\n",kstr, i);
    else
      printf ("%s\n", kstr);


    /* create data */
    memset (dstr, 0, dsize);
    if (!fixed_size)
      len = rand() % maxdsize;
    else
      len = dsize - 1;
    if (len < MIN_LEN)
      len = MIN_LEN;
    

    for (k = 0; k < len; k++) {
      tmp = rand() % LAST_CHAR;
      if (tmp < FIRST_CHAR)
	tmp += FIRST_CHAR;
      dstr[k] = (char)tmp;
    }
    if (!fixed_size)
      printf ("%s-%d\n",dstr, i);
    else
      printf ("%s\n", dstr);
  }
  
  free (kstr);
  free (dstr);
    
  return 0;
}
