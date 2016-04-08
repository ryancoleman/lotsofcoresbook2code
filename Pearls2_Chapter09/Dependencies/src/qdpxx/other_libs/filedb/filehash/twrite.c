/**
 * Write random strings into database
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <ffdb_db.h>


#define MIN_LEN 5
#define FIRST_CHAR (char)'0'
#define LAST_CHAR  (char)'~'

/**
 * create a new pair of keys and data with maximum size of keysize
 * of maxksize and maximum of datasize of maxdsize
 */
static void key_and_data (char key[], unsigned int maxksize,
			  char val[], unsigned int maxdsize)
{
  unsigned int len, k;
  int tmp;

  memset (key, 0, maxksize);
  len = rand() % maxksize;
  if (len < MIN_LEN)
    len = MIN_LEN + 1;

  for (k = 0; k < len - 1; k++) {
    tmp = rand () % LAST_CHAR;
    if (tmp < FIRST_CHAR)
      tmp += FIRST_CHAR;
    
    key[k] = (char)tmp;
  }
  key[k] = '\0';

  memset (val, 0, maxdsize);
  len = rand () % maxdsize;

  if (len < MIN_LEN)
    len = MIN_LEN + 1;

  for (k = 0; k < len - 1; k++) {
    tmp = rand() % LAST_CHAR;
    if (tmp < FIRST_CHAR)
      tmp += FIRST_CHAR;
    val[k] = (char)tmp;
  }
  val[k] = '\0';
}

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

#define MAX_LEN 100000

int main(int argc, char** argv)
{
  struct timeval tv;
  int num, maxksize, maxdsize, k;
  FFDB_DBT item, key;
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  char *p1, *p2, *dbase;
  char kstr[MAX_LEN], vstr[MAX_LEN];

  if (argc < 8) {
    fprintf (stderr, "Usage: %s bucketsize nbuckets rearrange(0|1) dbasename numkeys maxkeysize maxdatasize\n", argv[0]);
    exit(1);
  }

  int i = 0;
  
  argv++;
  ctl.hash = NULL;
  ctl.cmp = 0;
  ctl.cachesize = 0;
  ctl.bsize = atoi(*argv++);
  ctl.nbuckets = atoi(*argv++);
  ctl.rearrangepages = atoi(*argv++);
  ctl.numconfigs = 100;
  ctl.userinfolen = 100000;
  dbase = *argv++;

  num = atoi(*argv++);
  maxksize = atoi(*argv++);
  maxdsize = atoi(*argv++);

  if (maxksize > MAX_LEN) {
    fprintf (stderr, "Key string size exceeds maximum allowed value %d\n",
	     MAX_LEN);
    return -1;
  }

  if (maxdsize > MAX_LEN) {
    fprintf (stderr, "Data string size exceeds maximum allowed value %d\n",
	     MAX_LEN);
    return -1;
  }


  gettimeofday (&tv, 0);
  srand (tv.tv_sec);

#if 0
  if (!(dbp = ffdb_dbopen( dbase,
			   O_CREAT|O_TRUNC|O_RDWR, 0600, &ctl))){
    /* create table */
    fprintf(stderr, "cannot create: hash table (size %d)\n",
	    ctl.nbuckets);
    exit(1);
  }
#endif
  if (!(dbp = ffdb_dbopen( dbase,
			   O_CREAT|O_RDWR, 0600, &ctl))){
    /* create table */
    fprintf(stderr, "cannot create: hash table (size %d)\n",
	    ctl.nbuckets);
    exit(1);
  }

  /* Set user and configuration information */
  insert_user_info (dbp, USER_STRING);
  init_config_info (dbp, 100);

  key.data = kstr;
  item.data = vstr;
  for (i = 0; i < num; i++) {
    key_and_data (kstr, maxksize, vstr, maxdsize);
#if 0
    printf ("%s\n", kstr);
    printf ("%s\n", vstr);
#endif

    key.size = strlen(kstr) + 1;
    item.size = strlen(vstr) + 1;

    /**
     * enter key/data pair into the table
     */
    if ((dbp->put)(dbp, &key, &item, FFDB_NOOVERWRITE) != 0) {
      fprintf(stderr, "cannot enter: key %s\n",
	      (char *)item.data);
      exit(1);
    }
  }
  
  (dbp->close)(dbp);

  fprintf (stderr, "We put %d item in the database %s\n", i, dbase);

  exit(0);
}
