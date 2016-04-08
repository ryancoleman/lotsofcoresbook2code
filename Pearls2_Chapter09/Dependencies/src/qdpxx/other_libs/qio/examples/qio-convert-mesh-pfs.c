/* Host file conversion for grid machines.  FOR SINGLE PROCESSOR ONLY! */
/* Converts from SINGLEFILE to PARTFILE format. */
/* Here we treat the simplest case that each node reads its own file
   from its own directory path */

/* Usage

   qio-convert-mesh-pfs <part_sing> <filename> [<ildgLFN>]< layoutfile
 
    where

      <part_sing> = 0 to convert single to partition format
                  = 1 to convert partition to singlefile ILDG format
                  = 2 to convert partition to singlefile SciDAC native format
      <filename> is the base name of the file

    and

      layoutfile has the following format

      line 1: machdim            Number of machine dimensions allocated
      line 2: mx my mz ...       Dimensions of the allocated machine

         The remaining lines specify the host path to the file system
         for each logical node.  The first value is the logical node
         number and the second is the path.

      line 3: 0 /pfs/r22c0/R22/C0/B0/M0/D0/A0  
      line 4: 1 /pfs/r22c0/R22/C0/B0/M0/D0/A1
      line 5:       etc

*/

#include <qio.h>
#include <stdio.h>
#include "qio-convert-mesh.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define PATHLENGTH 256
#define	BASE_DIRMODE	0775

/* One-to-one */
static int self_io_node(int node){return node;}

static int zero_master_io_node(){return 0;}

static char *errmsg(void)
{
  return (strerror(errno) != NULL)?
    strerror(errno): "unknown error";
}

static QIO_Filesystem *create_multi_pfs(int numnodes){
  QIO_Filesystem *fs;
  int i, k;
  struct stat dir_stat;
  mode_t dir_mode = BASE_DIRMODE;

  /* Build the QIO file system structure */
  fs = (QIO_Filesystem *)malloc(sizeof(QIO_Filesystem));
  if(!fs){
    printf("Can't malloc fs\n");
    return NULL;
  }
  fs->number_io_nodes = numnodes;
  fs->type = QIO_MULTI_PATH;
  fs->my_io_node = self_io_node;
  fs->master_io_node = zero_master_io_node;
  fs->io_node = NULL;
  fs->node_path = NULL;

  /* Make room for the table of path names */
  fs->node_path = (char **)calloc(numnodes, sizeof(char *));
  if(!fs->node_path){
    printf("Path table malloc failed\n");
    return NULL;
  }
  
  for(i = 0; i < numnodes; i++){
    fs->node_path[i] = (char *)calloc(PATHLENGTH, sizeof(char));
    if(!fs->node_path[i]){
      printf("Path table malloc failed\n");
      return NULL;
    }
    fs->node_path[i][0] = '\0';
  }

  /* Read the table from stdin */
  /* Table entries have node number followed by directory path */
  for(i = 0; i < numnodes; i++){
    /* Read and check logical node number */
    if(scanf("%d",&k) != 1 || k < 0 || k >= numnodes){
      printf("Error reading path list\n");
      return NULL;
    }
    /* Read and create the corresponding directory path if need be */
    scanf("%s",fs->node_path[k]);
    if( mkdir(fs->node_path[k], dir_mode) < 0){
      if( errno != EEXIST ){
	printf("Can't make %s: %s\n", fs->node_path[k], errmsg());
	return NULL;
      }
    }
  }

  /* Check that we have an entry for each node */
  for(i = 0; i < numnodes; i++){
    if(strlen(fs->node_path[i]) == 0){
      printf("Missing path for node %d\n",i);
      return NULL;
    }
  }

  return fs;
}

static void destroy_multi_pfs(QIO_Filesystem *fs){
  int i;

  if(fs != NULL){
    if(fs->io_node != NULL)
      free(fs->io_node);
    if(fs->node_path != NULL){
      for(i = 0; i < fs->number_io_nodes; i++)
	if(fs->node_path[i] != NULL)
	  free(fs->node_path[i]);
      free(fs->node_path);
    }
    free(fs);
  }
}

int main(int argc, char *argv[]){
  QIO_Filesystem *fs;
  int status;
  QIO_Mesh_Topology *mesh;
  
  /* Check arguments and process layout parameters */

  if(argc < 3){
    fprintf(stderr,"Usage %s <0 (sing to part) | 1 (part to sing ILDG) | 2 (part to sing SciDAC) > <filename> [<ildgLFN>] < layoutfile\n",argv[0]);
    return 1;
  }

  /* Read topology */
  mesh = qio_read_topology(1);

  /* Create layout and file system structure */
  fs = create_multi_pfs(mesh->numnodes);
  if(!fs)return 1;

  /* Do the file conversion */
  status = qio_mesh_convert(fs, mesh, argc, argv);

  destroy_multi_pfs(fs);

  qio_destroy_topology(mesh);

  return status;
}
