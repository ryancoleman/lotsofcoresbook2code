/* Host file conversion for grid machines.  FOR SINGLE PROCESSOR ONLY! */
/* Converts from SINGLEFILE to PARTFILE format and back. */
/* Here we treat the more general case that nodes are grouped into I/O
   partitions and there is one file per partition from its own
   directory path. */

/* Usage

   qio-convert-mesh-ppfs <part_sing> <filename> [<ildgLFN>]< layoutfile
 
    where

      <part_sing> = 0 to convert single to partition format
                  = 1 to convert partition to singlefile ILDG format
                  = 2 to convert partition to singlefile SciDAC native format
      <filename> is the base name of the file

    and

      layoutfile has the following format

      line 1: machdim            Number of machine dimensions allocated
      line 2: mx my mz ...       Dimensions of the node-partitioned machine
      line 3: px py pz ...       Dimensions of the I/O partitioned machine

         The remaining lines specify the host path to the file system
         for each logical node.  The first value is the logical node
         number and the second is the path.  The path can be empty.

      line 4:  0 path1
      line 5:  k pathk
      line 6: 2k
      line 7:    etc

   Explanation of lines 2 and 3

      A lattice of dimension nx ny nz ... is divided into hypecubes
      with mx hypercubes in the x direction, my in the y direction, etc.
      Thus the sublattice dimension is nx/mx, ny/my, ... and the
      number of nodes is nx * ny * nz * ....

      For purposes of I/O the lattice is divided into possibly larger
      hypercubes with px in the x direction, py in the y direction, etc.
      The sublattice dimension for the I/O families is nx/px, ny/py, etc.
      Thus number of files is px * py * pz * ...

      We require that all nodes in an I/O family have their sites
      on a single file.  So px must be a factor of nx, py of ny, etc.  

      The factor k above is the number of nodes per file.
      It is px/nx * py/ny * pz/nz * ...

   Example:

      48^3 x 144 lattice with 

         mx, my, mz, mt = 2, 2, 2, 16

      and 

         px, py, pz, pt = 1, 1, 2, 16

      is laid out on 128 nodes with a sublattice dimension of 
      24^3 x 9.  There are 32 files with sublattice dimension
      48^2 x 24 x 9.

*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <qio.h>
#include "qio-convert-mesh.h"

#define PATHLENGTH 256
#define LINELENGTH PATHLENGTH+16
#define	BASE_DIRMODE	0775

static QIO_Mesh_Topology *mesh;
static int *nodes_per_ionode;
static int *io_node_coords;

/* Initialize my_io_node data */

static int init_my_io_node(){
  int i;
  char myname[] = "init_my_io_node";

  /* Make space for node coordinates */
  io_node_coords = lex_allocate_coords(mesh->machdim, myname);
  if(io_node_coords == NULL)return 1;

  /* Make space for consolidation factor along each direction */
  nodes_per_ionode = lex_allocate_coords(mesh->machdim, myname);
  if(nodes_per_ionode == NULL)return 1;

  /* Compute the number of nodes per I/O node along each direction */
  for(i = 0; i < mesh->machdim; i++)
    nodes_per_ionode[i] = mesh->machsize[i]/mesh->iomachsize[i];

  return 0;
}

/* Map any node to its I/O node */
static int my_io_node(int node){
  int i; 

  /* Get the machine coordinates for the specified node */
  lex_coords(io_node_coords, mesh->machdim, mesh->machsize, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(i = 0; i < mesh->machdim; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, mesh->machdim, mesh->machsize);
}

static int zero_master_io_node(){return 0;}

static char *errmsg(void)
{
  return (strerror(errno) != NULL)?
    strerror(errno): "unknown error";
}

static QIO_Filesystem *create_multi_ppfs(void){
  QIO_Filesystem *fs;
  int i, j, k, d, numnodes;
  int *io_part_coords;
  mode_t dir_mode = BASE_DIRMODE;
  char myname[] = "create_multi_ppfs";

  numnodes = mesh->numnodes;
  char line[LINELENGTH];
  char path[PATHLENGTH+1];
  char *p;

  /* Build the QIO file system structure */
  fs = (QIO_Filesystem *)malloc(sizeof(QIO_Filesystem));
  if(!fs){
    printf("Can't malloc fs\n");
    return NULL;
  }
  fs->number_io_nodes = mesh->number_io_nodes;
  fs->type = QIO_MULTI_PATH;
  fs->my_io_node = my_io_node;
  fs->master_io_node = zero_master_io_node;
  fs->io_node = NULL;
  fs->node_path = NULL;

  /* Make room for the table of I/O nodes */
  fs->io_node = (int *)malloc(mesh->number_io_nodes*sizeof(int));
  if(!fs->io_node){
    printf("Path table malloc failed\n");
    return NULL;
  }

  /* Make room for the table of path names */
  fs->node_path = (char **)calloc(mesh->number_io_nodes, sizeof(char *));
  if(!fs->node_path){
    printf("Path table malloc failed\n");
    return NULL;
  }
  
  for(i = 0; i < mesh->number_io_nodes; i++){
    fs->node_path[i] = (char *)calloc(PATHLENGTH+1, sizeof(char));
    if(!fs->node_path[i]){
      printf("Path table malloc failed\n");
      return NULL;
    }
    fs->node_path[i][0] = '\0';
  }

  /* Read the table from stdin */
  /* Table entries have node number followed by directory path */
  /* Directory path could be empty, so we read line-by-line */
  /* The table  could be empty, in which case we create it automatically */
  i = 0;
  while(i < fs->number_io_nodes){

    /* Read and check logical node number */
    p = fgets(line, LINELENGTH, stdin);
    if(p == NULL)break;  /* Usually EOF */

    /* Parse the line.  Skip it if we can't read a value. */
    j = sscanf(line, "%d %s",&k,path);
    if( j <= 0 )continue;

    /* Put the I/O node number in the table */
    fs->io_node[i] = k;

    /* If the path is missing, set the string to zero length */
    if(j == 1) fs->node_path[i][0] = '\0';
    /* Otherwise, copy the path to the table */
    else{
      strncpy(fs->node_path[i], path, PATHLENGTH);
      fs->node_path[i][PATHLENGTH] = '\0';
    }

    /* Create the directory */
    if(j > 1)
      if( mkdir(fs->node_path[i], dir_mode) < 0){
	if( errno != EEXIST ){
	  printf("Can't make %s: %s\n", fs->node_path[k], errmsg());
	  return NULL;
	}
      }
    i++;
  }

  /* If the table is empty, create it automatically */
  if(i == 0){
    /* Iterate over I/O partition coordinates */
    io_part_coords = lex_allocate_coords(mesh->machdim, myname);
    lex_init(&d, io_part_coords, mesh->machdim);
    i = 0;
    do{
      /* Convert I/O partition coordinate to node coordinate */
      for(j = 0; j < mesh->machdim; j++)
	io_node_coords[j] = io_part_coords[j]*nodes_per_ionode[j];
      /* Convert node coordinate to rank and store. Set path to null. */
      fs->io_node[i] = lex_rank(io_node_coords, mesh->machdim, 
				mesh->machsize);
      fs->node_path[i][0] = '\0';
      i++;
    } while(lex_next(&d, io_part_coords, mesh->machdim, mesh->iomachsize));

    free(io_part_coords);
  }

  /* Consistency check */
  if(i != fs->number_io_nodes){
    printf("Created %d node paths but expected %d\n",i,fs->number_io_nodes);
    return NULL;
  }

  return fs;
}

static void destroy_multi_ppfs(QIO_Filesystem *fs){
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
  
  /* Check arguments and process layout parameters */

  if(argc < 3){
    fprintf(stderr,"Usage %s <0 (sing to part) | 1 (part to sing ILDG) | 2 (part to sing SciDAC) > <filename> [<ildgLFN>] < layoutfile\n",argv[0]);
    return 1;
  }

  /* Read topology */
  mesh = qio_read_topology(0);

  /* Initialize the my_io_node function */
  status = init_my_io_node();
  if(status != 0)return status;

  /* Create layout and file system structure */
  fs = create_multi_ppfs();
  if(!fs)return 1;

  /* Do the file conversion */
  status = qio_mesh_convert(fs, mesh, argc, argv);

  destroy_multi_ppfs(fs);

  qio_destroy_topology(mesh);

  return status;
}
