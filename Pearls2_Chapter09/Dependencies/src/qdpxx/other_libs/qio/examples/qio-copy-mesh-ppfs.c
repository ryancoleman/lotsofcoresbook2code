/* MPP file copy program */
/* Copies partfiles to the appropriate nodes from source directories. */
/* Here we treat the more general case that nodes are grouped into I/O
   partitions and there is one file per partition from its own
   directory path. */

/* Usage

    qio-copy-mesh-ppfs <srcfilename> <destdir> <delay> < layoutfile
 
    where

      <srcfilename> is the base name of the file
      <destdir> is the destination directory
      <delay> is the delay in seconds between file copy initiation

    and

      layoutfile has the following format

      line 1: machdim            Number of machine dimensions allocated
      line 2: mx my mz ...       Dimensions of the node-partitioned machine
      line 3: px py pz ...       Dimensions of the I/O partitioned machine

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

   Example:

      A lattice with 

         mx, my, mz, mt = 2, 2, 2, 16

      and 

         px, py, pz, pt = 1, 1, 2, 16

      is laid out on 128 nodes.  There are 32 files.  The first
      file has data for nodes 0 1 2 3 and will be read by node 0.
      The secone file has data for nodes 4 5 6 7 and will be read
      by node 4, etc.

*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <qio.h>
#include <qmp.h>
#include "qio-convert-mesh.h"

#define PATHLENGTH 256
#define LINELENGTH PATHLENGTH+16
#define	BASE_DIRMODE	0775

static QIO_Mesh_Topology *mesh;
static int *nodes_per_ionode;
static int *io_node_coords;

/* Initialize my_io_node data */

static int init_my_io_node(void){
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
static int
my_io_node(int node)
{
  /* Get the machine coordinates for this node's IO node */
  lex_coords(io_node_coords, mesh->machdim, mesh->machsize, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(int i = 0; i < mesh->machdim; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);

  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, mesh->machdim, mesh->machsize);
}

static int zero_master_io_node(){return 0;}

static QIO_Filesystem *
create_abbrev_multi_ppfs(void)
{
  /* Build a lean QIO file system structure */
  QIO_Filesystem *fs = (QIO_Filesystem *) malloc(sizeof(QIO_Filesystem));
  if(fs==NULL) {
    printf("Can't malloc fs\n"); fflush(stdout);
    return NULL;
  }
  fs->number_io_nodes = mesh->number_io_nodes;
  fs->type = QIO_MULTI_PATH;
  fs->my_io_node = my_io_node;
  fs->master_io_node = zero_master_io_node;
  fs->io_node = NULL;
  fs->node_path = NULL;
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

static void
init_qmp(int *argc, char **argv[])
{
  QMP_status_t status;
  QMP_thread_level_t req, prv;

  /* Start QMP */
  req = QMP_THREAD_SINGLE;
  status = QMP_init_msg_passing (argc, argv, req, &prv);

  if (status != QMP_SUCCESS) {
    QMP_error ("QMP_init failed: %s\n", QMP_error_string(status));
    QMP_abort(1);
  }
}

static void quit_qmp(void){
  QMP_finalize_msg_passing ();
}

#define MAXBUF  1048576

static int qio_file_copy(QIO_Filesystem *fs, int *argc, char **argv[])
{

  char *srcfilestem, *srcfilepath, *srcfilename;
  char *dstdir, *dstfilepath;
  double delay;
  int this_node;
  int status = 0;
  size_t rbytes, wbytes;
  FILE *srcfp, *dstfp;
  char *buf;
  char myname[] = "qio_file_copy";

  /* Command line args */
  sscanf((*argv)[2],"%lf",&delay);
  srcfilestem = (*argv)[3];
  dstdir = (*argv)[4];

  /* Initialize QMP */
  init_qmp(argc, argv);

  this_node = QMP_get_node_number();

  /* Delay according to the node number */
  QIO_wait(this_node*delay);

  /* If I am my own I/O node, then I copy the file */
  if(this_node == fs->my_io_node(this_node)){

    /* Create the full source file name with the volume number suffix */
    srcfilepath = QIO_filename_edit(srcfilestem, QIO_PARTFILE, this_node);

    /* Create the destination file name */
    srcfilename = strrchr(srcfilepath,'/') + 1;
    if(srcfilename == NULL)srcfilename = srcfilepath;
    dstfilepath = (char *)malloc(strlen(dstdir) + 1 + 
				  strlen(srcfilename) + 1);
    sprintf(dstfilepath,"%s/%s",dstdir,srcfilename);

    /* Copy the file */

    printf("%s(%d) copying %s to %s\n",myname,this_node,srcfilepath,dstfilepath);fflush(stdout);

    buf = (char *)malloc(MAXBUF);
    if(buf == NULL){
      printf("qio_file_copy: Can't malloc buf\n");fflush(stdout);
      status = 1;
    }

    srcfp = DCAPL(fopen)(srcfilepath,"r");
    if(srcfp == NULL){
      printf("qio_file_copy: Can't open %s\n",srcfilepath);fflush(stdout);
      status = 1;
    }

    dstfp = DCAPL(fopen)(dstfilepath,"w");
    if(dstfp == NULL){
      printf("qio_file_copy: Can't open %s\n",dstfilepath);fflush(stdout);
      status = 1;
    }

    /* Copy loop */
    while(status == 0){

      rbytes = DCAP(fread)(buf, sizeof(unsigned char), (size_t)MAXBUF, srcfp);
      if(rbytes == 0)break;
      
      wbytes = DCAP(fwrite)(buf, sizeof(unsigned char), rbytes, dstfp);

      if(rbytes != wbytes){
	printf("qio_file_copy: Error writing %ld bytes\n",(long)rbytes);
	fflush(stdout);
	status = 1; 
	break;
      }
      if(rbytes < MAXBUF)break;
    }

    if(srcfilepath != NULL) free(srcfilepath);
    if(dstfilepath != NULL) free(dstfilepath);
    if(srcfp != NULL)       DCAP(fclose)(srcfp);
    if(dstfp != NULL)       DCAP(fclose)(dstfp);
  }

  /* Synchronize */
  QMP_sum_int(&status);

  /* Finalize QMP */
  quit_qmp();

  return status;
}

int main(int argc, char *argv[]){
  QIO_Filesystem *fs;
  int status;
  FILE *fp;
  
  /* Check arguments and process layout parameters */

  if(argc < 5){
    fprintf(stderr,"Usage %s <layoutfile> <delay> <srcfilename> <destdir> \n",
	    argv[0]);
    return 1;
  }

  /* Take stdin from file given by the 1st arg on command line */
  fp = freopen(argv[1],"r",stdin);
  if(fp == NULL){
    fprintf(stderr,"Can't open stdin file %s for reading.\n",argv[1]);
    return 1;
  }

  /* Read topology */
  mesh = qio_read_topology(0);

  /* Initialize the my_io_node function */
  status = init_my_io_node();
  if(status != 0)return status;

  /* Create layout and file system structure */
  fs = create_abbrev_multi_ppfs();
  if(!fs)return 1;

  /* Do the file copy */
  status = qio_file_copy(fs, &argc, &argv);

  destroy_multi_ppfs(fs);

  qio_destroy_topology(mesh);

  return status;
}
