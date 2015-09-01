/**
 * Simple code to test read database
 */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <ffdb_db.h>

#if 0
#define INITIAL	500000
#define MAXWORDS 500000	       /* # of elements in search table */
#endif

#define INITIAL	50000
#define MAXWORDS 50000	       /* # of elements in search table */


static void
read_user_info (FFDB_DB* dbp, unsigned char string[], unsigned int* len)
{
  ffdb_get_user_info (dbp, string, len);
  fprintf (stderr, "User information %d \n", *len);
  fprintf (stderr, "%s\n", string);
}

static void
read_config_info (FFDB_DB* dbp, ffdb_all_config_info_t* configs)
{
  unsigned int i;

  ffdb_get_all_configs (dbp, configs);

  fprintf (stderr, "Number of configurations = %d\n", configs->numconfigs);
  for (i = 0; i < configs->numconfigs; i++) {
    fprintf (stderr, "config %d index %d inserted %d type %d file %s\n",
	     configs->allconfigs[i].config,
	     configs->allconfigs[i].index,
	     configs->allconfigs[i].inserted,
	     configs->allconfigs[i].type,
	     configs->allconfigs[i].fname);
  }

  free (configs->allconfigs);
}

typedef struct {		       /* info to be stored */
  int num, siz;
} info;

#define MAX_LEN 32768

int main(int argc, char** argv)
{
  FFDB_DBT key, res;
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  ffdb_cursor_t* cur;
  int	stat, i, numkey, numpairs;
  char *dbase;
  unsigned char userdata[4096];
  unsigned int len = 4096;
  ffdb_all_config_info_t acf;
  char *mykey[MAXWORDS], *mydata[MAXWORDS], tmp1[MAX_LEN], tmp2[MAX_LEN];
  FILE *fd;

  if (argc < 4) {
    fprintf (stderr, "Usage: %s cachesize dbase stringfiles\n", argv[0]);
    exit (1);
  }

  fd = fopen (argv[3], "r");
  if (!fd) {
    fprintf (stderr, "Cannot open string file %s\n", argv[3]);
    return -1;
  }

  i = 0;
  while (fscanf (fd, "%s", tmp1) >= 1 &&
	 fscanf (fd, "%s", tmp2) >= 1 && i < MAXWORDS) {
    mykey[i] = (char *)malloc(strlen(tmp1) + 1);
    strcpy (mykey[i], tmp1);
    mydata[i] = (char *)malloc(strlen(tmp2) + 1);
    strcpy (mydata[i], tmp2);
    i++;
  }
  numpairs = i;
  fclose (fd);


  i = 0;
  argv++;
  ctl.nbuckets = INITIAL;
  ctl.hash = NULL;
  ctl.cmp = NULL;
  ctl.bsize = 64;
  ctl.cachesize = atoi(*argv++);
  ctl.rearrangepages = 0;
  dbase = *argv++;
  fprintf (stderr, "dbase = %s\n", dbase);
  if (!(dbp = ffdb_dbopen(dbase, O_RDONLY, 0600, &ctl))) {
    /* create table */
    fprintf(stderr, "cannot open: hash table\n" );
    exit(1);
  }

  /* Read user and configuration information */
  read_user_info (dbp, userdata, &len);

  read_config_info (dbp, &acf);

  /* create a new cursor */
  stat = dbp->cursor (dbp, &cur, FFDB_KEY_CURSOR);
  if (stat != 0) {
    fprintf (stderr, "Cannot open a cursor\n");
    (dbp->close)(dbp);
    return -1;
  }

#if 0
  key.data = 0;
  key.size = 0;
  res.data = 0;
  res.size = 0;
#endif
  key.data = tmp1;
  key.size = MAX_LEN;
  res.data = tmp2;
  res.size = MAX_LEN;

  numkey = 0;
  while ((stat = cur->get (cur, &key, &res, FFDB_NEXT)) == FFDB_SUCCESS) {
#if 0
    fprintf (stderr, "Key len %d = %s\n", strlen((char *)(key.data)),(char *)(key.data));
    fprintf (stderr, "Data %s\n", (char *)(res.data));
#endif
    numkey++;

    /* Look for the key in the string table */
    for (i = 0; i < numpairs; i++) {
      if (strcmp (mykey[i], (char *)(key.data)) == 0) {
	if (strcmp (mydata[i], (char *)(res.data)) != 0) {
	    fprintf (stderr, "Data mismatch.\n");
	    return -1;
	}
      }
      break;
    }
    if (i >= numpairs) {
      fprintf (stderr, "Cannot match any keys\n");
      return -1;
    }

#if 0	  
    free (key.data);
    free (res.data);

    key.data = 0;
    key.size = 0;
    res.data = 0;
    res.size = 0;
#endif

    key.data = tmp1;
    key.size = MAX_LEN;
    res.data = tmp2;
    res.size = MAX_LEN;
  }

  cur->close (cur);

  (dbp->close)(dbp);

  fprintf (stderr, "Number of pairs found = %d\n", numkey);

  for (i = 0; i < numpairs; i++) {
    free (mykey[i]);
    free (mydata[i]);
  }

  return 0;
}
