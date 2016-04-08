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

#define MAX_LEN 100000

char	wp1[MAX_LEN];
char	wp2[MAX_LEN];

int main(int argc, char** argv)
{
  FFDB_DBT item, key, res;
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  int	stat, i;
  char *p1, *p2, *dbase;
  unsigned char userdata[4096];
  unsigned int len = 4096;
  unsigned char line[MAX_LEN];
  ffdb_all_config_info_t acf;

  if (argc < 4) {
    fprintf (stderr, "Usage: %s cachesize rearrange(0|1) dbase\n", argv[0]);
    exit (1);
  }

  i = 0;
  
  argv++;
  ctl.nbuckets = INITIAL;
  ctl.hash = NULL;
  ctl.cmp = NULL;
  ctl.bsize = 64;
  ctl.cachesize = atoi(*argv++);
  ctl.rearrangepages = atoi(*argv++);
  dbase = *argv++;
  fprintf (stderr, "dbase = %s\n", dbase);
  if (!(dbp = ffdb_dbopen(dbase, O_RDWR, 0600, &ctl))) {
    /* create table */
    fprintf(stderr, "cannot open: hash table\n" );
    exit(1);
  }

  /* Read user and configuration information */
  read_user_info (dbp, userdata, &len);

  read_config_info (dbp, &acf);

  key.data = wp1;
  item.data = wp2;
  while ( fgets(wp1, MAX_LEN, stdin) &&
	  fgets(wp2, MAX_LEN, stdin) &&
	  i++ < MAXWORDS) {
    /*
     * put info in structure, and structure in the item
     */
    p1 = wp1;
    while (p1 && *p1 != '\n')
      p1++;
    *p1 = '\0';

    p2 = wp2;
    while (p2 && *p2 != '\n')
      p2++;
    *p2 = '\0';

    key.size = strlen(wp1) + 1;
    item.size = strlen(wp2) + 1;

#if 0
    /* clear out result */
    res.data = 0;
    res.size = 0;
#endif
    res.data = line;
    res.size = MAX_LEN;

    stat = (dbp->get)(dbp, &key, &res, 0);
    if (stat < 0) {
      fprintf ( stderr, "Error retrieving %s\n", (char *)key.data );
      exit(1);
    } else if ( stat > 0 ) {
      fprintf ( stderr, "%s not found\n", (char *)key.data );
      exit(1);
    }
    else {
      if (strcmp ((char *)item.data, (char *)res.data) != 0) 
	fprintf (stderr, "Retriving data mismatch %s != %s (expeced)\n",
		 (char *)res.data, (char *)item.data);
#if 0
      free (res.data);
#endif
    }
  }
  (dbp->close)(dbp);

  return 0;
}
