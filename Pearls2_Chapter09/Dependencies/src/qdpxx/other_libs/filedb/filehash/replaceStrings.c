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

#define MAX_LEN 32768

int
main (int argc, char** argv)
{
  struct timeval tv;
  char kstr[MAX_LEN], dstr[MAX_LEN], numstr[32];
  char *p1, *p2, tmp;
  FILE* fd;
  unsigned int len, k, i;

  if (argc < 2) {
    fprintf (stderr, "Usage: %s file\n", argv[0]);
    return 1;
  }

  fd = fopen (argv[1], "w+");
  if (!fd) {
    fprintf (stderr, "Cannot open replacement string file %s\n", argv[1]);
    return -1;
  }

  gettimeofday (&tv, 0);
  srand (tv.tv_sec);

  i = 0;
  while (fgets(kstr, MAX_LEN, stdin) &&
	 fgets(dstr, MAX_LEN, stdin)) {
    /* create number string */
    sprintf (numstr, "-%d", i);

    p1 = kstr;
    while (p1 && *p1 != '\n')
      p1++;
    *p1 = '\0';

    p2 = dstr;
    while (p2 && *p2 != '\n')
      p2++;
    *p2 = '\0';
    len = strlen(dstr);

    /* create random string to replace data */
    memset (dstr, 0, MAX_LEN);
    len = rand() % (len - strlen(numstr));
    if (len < MIN_LEN)
      len = MIN_LEN;
    
    for (k = 0; k < len; k++) {
      tmp = rand() % LAST_CHAR;
      if (tmp < FIRST_CHAR)
	tmp += FIRST_CHAR;
      dstr[k] = (char)tmp;
    }
    
    fprintf (fd, "%s\n",kstr);
    fprintf (fd, "%s%s\n",dstr,numstr);
    i++;
  }
  
  fclose (fd);
    
  return 0;
}
