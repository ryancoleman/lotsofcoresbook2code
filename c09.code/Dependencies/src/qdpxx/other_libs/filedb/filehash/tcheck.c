/**
 * Simple code to check database
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

#define INITIAL	500000
#define MAXWORDS 500000	       /* # of elements in search table */


static void
read_user_info (FFDB_DB* dbp, unsigned char string[], unsigned int* len)
{
  int status = ffdb_get_user_info (dbp, string, len);
  if (status != 0)
    exit (1);
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

char	wp1[MAX_LEN];
char	wp2[MAX_LEN];

int main(int argc, char** argv)
{
  FFDB_DBT item, key, res;
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  int	stat, i, numkey;
  char *p1, *p2, *dbase;
  unsigned char userdata[4096];
  unsigned int len = 4096;
  unsigned char line[MAX_LEN];
  ffdb_all_config_info_t acf;
  ffdb_cursor_t* cur;

  if (argc < 3) {
    fprintf (stderr, "Usage: %s cachesize dbase\n", argv[0]);
    exit (1);
  }

  i = 0;
  
  argv++;
  ctl.nbuckets = INITIAL;
  ctl.hash = NULL;
  ctl.cmp = NULL;
  ctl.bsize = 8192;
  ctl.cachesize = atoi(*argv++);
  ctl.rearrangepages = 0;
  dbase = *argv++;
  fprintf (stderr, "dbase = %s\n", dbase);
#if 1
  if (!(dbp = ffdb_dbopen(dbase, O_RDONLY, 0400, &ctl))) {
    /* create table */
    fprintf(stderr, "cannot open: hash table\n" );
    exit(1);
  }
#endif

#if 0
  if (!(dbp = ffdb_dbopen(dbase, O_RDWR, 0600, &ctl))) {
    /* create table */
    fprintf(stderr, "cannot open: hash table\n" );
    exit(1);
  }
#endif

#if 0
  /* Read user and configuration information */
  read_user_info (dbp, userdata, &len);

  read_config_info (dbp, &acf);
#endif

  /* Walk through all keys */
  /* create a new cursor */
  stat = dbp->cursor (dbp, &cur, FFDB_KEY_CURSOR);
  if (stat != 0) {
    fprintf (stderr, "Cannot open a cursor\n");
    (dbp->close)(dbp);
    return -1;
  }

  key.data = 0;
  key.size = 0;
  res.data = 0;
  res.size = 0;

  numkey = 0;
  while ((stat = cur->get (cur, &key, &res, FFDB_NEXT)) == FFDB_SUCCESS) {
#if 0
    fprintf (stderr, "Key len %d = %s\n", strlen((char *)(key.data)),(char *)(key.data));
    fprintf (stderr, "Data %s\n", (char *)(res.data));
#endif
    numkey++;  

#if 0
    fprintf (stderr, "Key len = %d result len = %d\n", key.size, res.size);
#endif
    if (key.size > 0)
      free (key.data);
    if (res.size > 0)
      free (res.data);

    key.data = 0;
    key.size = 0;
    res.data = 0;
    res.size = 0;
  }

  fprintf (stderr, "Number of keys = %d\n", numkey);

  (dbp->close)(dbp);

  return 0;
}
