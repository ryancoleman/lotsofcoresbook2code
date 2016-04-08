/* Host file conversion for grid machines.  FOR SINGLE PROCESSOR ONLY! */
/* Converts from SINGLEFILE to PARTFILE format. */
/* Here we treat the simplest case that each node reads its own file
   from the same directory path */

/* Usage

   qio-convert-mesh-singlefs <part_sing> <filename> [<ildgLFN>]< layoutfile
 
    where

      <part_sing> = 0 to convert single to partition format
                  = 1 to convert partition to singlefile ILDG format
                  = 2 to convert partition to singlefile SciDAC native format
      <filename> is the base name of the file

    and

      layoutfile has the following format

      line 1: machdim            Number of machine dimensions allocated
      line 2: mx my mz ...       Dimensions of the allocated machine

*/

#include <qio.h>
#include <stdio.h>
#include "qio-convert-mesh.h"

static int self_io_node(int node){return node;}

static int zero_master_io_node(){return 0;}

static QIO_Filesystem *create_simple_fs(int numnodes){
  QIO_Filesystem *fs;

  fs = (QIO_Filesystem *)malloc(sizeof(QIO_Filesystem));
  if(!fs){
    printf("Can't malloc fs\n");
    return NULL;
  }
  fs->number_io_nodes = numnodes;
  fs->type = QIO_SINGLE_PATH;
  fs->my_io_node = self_io_node;
  fs->master_io_node = zero_master_io_node;
  fs->io_node = NULL;
  fs->node_path = NULL;

  return fs;
}

static void destroy_simple_fs(QIO_Filesystem *fs){
  free(fs);
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
  fs = create_simple_fs(mesh->numnodes);
  if(!fs)return 1;

  status = qio_mesh_convert(fs, mesh, argc, argv);

  destroy_simple_fs(fs);

  qio_destroy_topology(mesh);

  return status;
}
