/**
 * Simple code to test threaded read database
 */
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <pthread.h>
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

typedef struct _thread_info_
{
  FFDB_DB* dbp;
  int      numthread;
  char*    fname;
  int      tnum;
}thread_info_t;

#define MAX_LEN 32768

#define MAX_THREADS 100

static void *reader_thread (void* arg)
{
  thread_info_t *tinfo = (thread_info_t *)arg;
  FILE* fd;
  char	wp1[MAX_LEN], wp2[MAX_LEN], recv[MAX_LEN];
  FFDB_DBT key, item, res;
  int i, stat;


  fd = fopen (tinfo->fname, "r");
  if (!fd) {
    fprintf (stderr, "Cannot open string file %s\n", tinfo->fname);
    free (tinfo);
    return 0;
  }

  i = 0;
  while (fscanf (fd, "%s", wp1) >= 1 &&
	 fscanf (fd, "%s", wp2) >= 1 && i < MAXWORDS) {
    key.data = wp1;
    item.data = wp2;
    key.size = strlen(wp1) + 1;
    item.size = strlen(wp2) + 1;

#if 0
    /* clear out result */
    res.data = 0;
    res.size = 0;
#endif
    res.data = recv;
    res.size = MAX_LEN;

    stat = (tinfo->dbp->get)(tinfo->dbp, &key, &res, 0);
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
    i++;
  }
  fprintf (stderr, "Done thread %d\n", tinfo->tnum);

  free (tinfo);
  return 0;
}

int main(int argc, char** argv)
{
  FFDB_DB	*dbp;
  FFDB_HASHINFO ctl;
  int  i, numthread;
  char *dbase, *strfile;
  unsigned char userdata[4096];
  unsigned int len = 4096;
  ffdb_all_config_info_t acf;
  thread_info_t* tinfo;
  pthread_t tid[MAX_THREADS];
  void* status;

  if (argc < 5) {
    fprintf (stderr, "Usage: %s cachesize numthreads dbase stringfile\n", argv[0]);
    exit (1);
  }

  i = 0;
  
  argv++;
  ctl.nbuckets = INITIAL;
  ctl.hash = NULL;
  ctl.cmp = NULL;
  ctl.bsize = 64;
  ctl.cachesize = atoi(*argv++);
  ctl.rearrangepages = 0;
  numthread = atoi(*argv++);
  dbase = *argv++;
  strfile = *argv++;
  fprintf (stderr, "dbase = %s number thread = %d\n", dbase, numthread);
  if (!(dbp = ffdb_dbopen(dbase, O_RDONLY, 0600, &ctl))) {
    /* create table */
    fprintf(stderr, "cannot open: hash table\n" );
    exit(1);
  }

  /* Read user and configuration information */
  read_user_info (dbp, userdata, &len);

  read_config_info (dbp, &acf);

  /* Fireup several threads */
  for (i = 0; i < numthread; i++) {
    tinfo = (thread_info_t *)malloc(sizeof(thread_info_t));
    tinfo->dbp = dbp;
    tinfo->fname = strfile;
    tinfo->numthread = numthread;
    tinfo->tnum = i;
    pthread_create (&tid[i], 0, reader_thread, (void *)tinfo);
  }


  /* Wait for all threads to finish */
  for (i = 0; i < numthread; i++)
    pthread_join (tid[i], &status);

  dbp->close (dbp);

  return 0;
}

