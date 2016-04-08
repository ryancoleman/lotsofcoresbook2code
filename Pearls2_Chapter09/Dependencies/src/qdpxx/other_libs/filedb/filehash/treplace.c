/**
 * Simple Test to insert pairs into database
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <ffdb_db.h>

#if 0
#define INITIAL	    500000
#define MAXWORDS    500000	              /* # of elements in search table */
#endif

#define INITIAL	    500000
#define MAXWORDS    500000	              /* # of elements in search table */

#define MAX_LEN 32768

#define USER_STRING "Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCDHello threre I am the new user string for QCD  Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD"

static void
insert_user_info (FFDB_DB* dbp, char* string)
{
  ffdb_set_user_info (dbp, (unsigned char *)string, strlen (string) + 1);
}

static void
init_config_info (FFDB_DB* dbp, unsigned int numconfig)
{
  ffdb_all_config_info_t acf;
  unsigned int i;

  acf.numconfigs = numconfig;
  acf.allconfigs = (ffdb_config_info_t *)calloc(1, numconfig * sizeof(ffdb_config_info_t));
  
  for (i = 0; i < numconfig; i++) {
    acf.allconfigs[i].config = i;
    acf.allconfigs[i].index = i;
    acf.allconfigs[i].inserted = 0;
    acf.allconfigs[i].type = 0;
    acf.allconfigs[i].mtime = 0;
    strcpy(acf.allconfigs[i].fname, "N/A");
  }

  ffdb_set_all_configs (dbp, &acf);

  free (acf.allconfigs);
}

char	wp1[MAX_LEN];
char	wp2[MAX_LEN];
int main(int argc, char** argv)
{
  FFDB_DBT item, key;
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  char *p1, *p2, *dbase;

  if (argc < 4) {
    fprintf (stderr, "Usage: %s cachesize rearrange(0|1) dbasename\n", argv[0]);
    exit(1);
  }

  int i = 0;
  
  argv++;
  ctl.hash = NULL;
  ctl.cmp = 0;
  ctl.cachesize = atoi(*argv++);
  ctl.bsize = 0;
  ctl.nbuckets = 0;
  ctl.rearrangepages = atoi(*argv++);
  ctl.numconfigs = 100;
  dbase = *argv++;

  if (!(dbp = ffdb_dbopen( dbase,
			   O_CREAT|O_RDWR, 0600, &ctl))){
    /* create table */
    fprintf(stderr, "cannot create: hash table (size %d)\n",
	    INITIAL);
    exit(1);
  }

  /* Set user and configuration information */
  insert_user_info (dbp, USER_STRING);
  init_config_info (dbp, 100);

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
    
    /*
     * enter key/data pair into the table
     */
    if ((dbp->put)(dbp, &key, &item, 0) != 0) {
      fprintf(stderr, "cannot enter: key %s\n",
	      (char *)item.data);
      exit(1);
    }
  }
  
  (dbp->close)(dbp);

  fprintf (stderr, "We put %d item in the database %s\n", i, dbase);

  exit(0);
}
