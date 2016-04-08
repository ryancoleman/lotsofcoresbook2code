/* Convert between SINGLEFILE and PARTFILE format.
   on single processor. */

#include <qio.h>
#include <stdio.h>
#define MAIN
#include "qio-convert-mesh.h"

/*--------------------------------------------------------------------*/

QIO_Layout *create_mpp_layout(int numnodes, int *latsize_in, int latdim){
  int i;
  size_t volume;
  QIO_Layout *layout;
  int *latsize;

  layout = (QIO_Layout *)malloc(sizeof(QIO_Layout));
  if(!layout){
    printf("Can't malloc layout\n");
    return NULL;
  }

  /* Make local copy of lattice size */
  if(latsize_in != NULL){
    latsize = (int *)malloc(latdim*sizeof(int));
    if(!latsize){
      printf("Can't malloc latsize\n");
      return NULL;
    }
    for(i = 0; i < latdim; i++)
      latsize[i] = latsize_in[i];
  }
  else
    latsize = NULL;

  volume = 1;
  for(i = 0; i < latdim; i++)
    volume *= latsize[i];

  layout->node_number = node_number;
  layout->node_index = node_index;
  layout->get_coords = get_coords;
  layout->num_sites = num_sites;
  layout->latsize = latsize;
  layout->latdim = latdim;
  layout->volume = volume;
  layout->sites_on_node = 0;  /* Ignored */
  layout->this_node = 0;      /* Reset */
  layout->number_of_nodes = numnodes;
  return layout;
}

void destroy_mpp_layout(QIO_Layout *layout){
  free(layout);
}

/* If onetoone = 1 then each node is its own I/O node and we don't
   read the I/O machine dimensions */

int qio_mesh_convert(QIO_Filesystem *fs, QIO_Mesh_Topology *mesh,
		     int argc, char *argv[])
{

  /* For the case that each nodes reads its own file from the same
     directory path */

  QIO_Layout *mpp_layout;

  int n;
  int status;
  int numnodes = mesh->numnodes;
  int part_to_single;
  int latdim;
  int *latsize;
  QIO_Reader *qio_in;
  char *filename;
  char *newfilename;
  char *stringLFN;
  QIO_String *ildgLFN;

  QIO_verbose(QIO_VERB_REG);

  /* Command line options */

  n = 1;

  /* Which direction to convert? */
  part_to_single = atoi(argv[n++]);

  /* File name */
  filename = argv[n++];

  /* ildgLFN, if specified. Used only when recombining a file. */
  if(argc > n)
    stringLFN = argv[n++];
  else
    stringLFN = NULL;

  ildgLFN = QIO_string_create();
  QIO_string_set(ildgLFN, stringLFN);

  /* Start from a dummy layout */
  mpp_layout = create_mpp_layout(numnodes, NULL, 0);
  if(!mpp_layout)return 1;

  if(part_to_single == 0){
    /* If we are converting single to partfile the input master file is
       the current single file */
    qio_in = QIO_open_read_master(filename, mpp_layout, 0, fs->my_io_node,
				  fs->master_io_node);
  }
  else{
    /* Otherwise the input master file is a partfile */
    /* Set input path for file according to MULTI/SINGLE PATH flag */
    newfilename = QIO_set_filepath(fs,filename,fs->master_io_node());
    
    /* Get lattice dimensions from file */
    qio_in = QIO_open_read_master(newfilename, mpp_layout, 0, fs->my_io_node,
				  fs->master_io_node);
    free(newfilename);
    if(!qio_in)return 1;
  }

  latdim = QIO_get_reader_latdim(qio_in);
  latsize = QIO_get_reader_latsize(qio_in);

  /* Now create the real layout functions */
  mpp_layout = create_mpp_layout(numnodes, latsize, latdim);
  if(setup_layout(latsize, latdim,mesh->machsize, mesh->machdim)){
    printf("Error in setup_layout\n");
    return 1;
  }

  if(!mpp_layout)return 1;

  /* Close the file */
  QIO_close_read(qio_in);

  /* Do the conversion */
  if(part_to_single == 0)
    {
      printf("Converting %s from SINGLEFILE to PARTFILE\n",filename);
      status = QIO_single_to_part(filename, fs, mpp_layout);
    }
  else if(part_to_single == 1)
    {
      /* ILDG compatible format */
      printf("Converting %s from PARTFILE to SINGLEFILE ILDG\n",filename);
      status = QIO_part_to_single(filename, QIO_ILDGLAT, ildgLFN, 
				  fs, mpp_layout);
    }
  else
    {
      /* SciDAC native format */
      printf("Converting %s from PARTFILE to SINGLEFILE SciDAC\n",filename);
      status = QIO_part_to_single(filename, QIO_ILDGNO, NULL, fs, mpp_layout);
    }

  /* Clean up */
  destroy_mpp_layout(mpp_layout);

  return status;
}

