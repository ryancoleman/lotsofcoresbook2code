/* Convert between SINGLEFILE and PARTFILE format.
   on single processor. */

#include <stdio.h>
#include <stdlib.h>
#include <qio.h>
#define MAIN
#include "qio-convert-mesh.h"

/*---------------------------------------------------------------------*/
/* A set of utilities for counting coordinates in lexicographic order  */
/*---------------------------------------------------------------------*/

/* Initialize the coordinate counter */
void lex_init( int *dimp, int coords[], int dim )
{
  int d;
  for(d = 0; d < dim; d++)coords[d] = 0;
  *dimp = 0;
}

/* Increase the coordinate counter by one */
int lex_next(int *dimp, int coords[], int dim, int size[])
{
  if(++coords[*dimp] < size[*dimp]){
    *dimp = 0;
    return 1;
  }
  else{
    coords[*dimp] = 0;
    if(++(*dimp) < dim)return lex_next(dimp, coords, dim, size);
    else return 0;
  }
}

/*------------------------------------------------------------------*/
/* Convert linear lexicographic rank to lexicographic coordinate */

void lex_coords(int coords[], const int dim, const int size[], 
		    const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   lex_coords) */

size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Make temporary space for coords */

int *lex_allocate_coords(int dim, char *myname){
  int *coords;

  coords = (int *)malloc(dim*sizeof(int));
  if(!coords)printf("%s can't malloc coords\n",myname);
  return coords;
}

QIO_Mesh_Topology *
qio_read_topology(int onetoone)
{
  QIO_Mesh_Topology *mesh =
    (QIO_Mesh_Topology *) malloc(sizeof(QIO_Mesh_Topology));

  /* Read the number of lattice dimensions */
  int nscan = scanf("%d", &mesh->machdim);
  if(nscan!=1) {
    printf("Can't parse machdim from stdin\n");
    return NULL;
  }
  mesh->machsize = (int *)malloc(mesh->machdim*sizeof(int));
  if(mesh->machsize == NULL) {
    printf("Can't malloc machsize with dim %d\n",mesh->machdim);
    return NULL;
  }

  /* Read the machine size */
  mesh->numnodes = 1;
  for(int i=0; i < mesh->machdim; i++){
    nscan = scanf("%d",&mesh->machsize[i]);
    if(nscan!=1) {
      printf("Can't parse machsize[%i] from stdin\n", i);
      return NULL;
    }
    mesh->numnodes *= mesh->machsize[i];
  }

  /* Read the I/O machine size if requested */
  if(onetoone){
    mesh->number_io_nodes = mesh->numnodes;
    mesh->iomachsize = mesh->machsize;
  }
  else{
    mesh->iomachsize = (int *)malloc(mesh->machdim*sizeof(int));
    if(mesh->iomachsize == NULL){
      printf("Can't malloc machsize with dim %d\n",mesh->machdim);
      return NULL;
    }
    mesh->number_io_nodes = 1;
    for(int i = 0; i < mesh->machdim; i++){
      if(scanf("%d",&mesh->iomachsize[i]) != 1){
	printf("Missing I/O machine dimension\n");
	return NULL;
      }
      mesh->number_io_nodes *= mesh->iomachsize[i];
      if(mesh->machsize[i] % mesh->iomachsize[i] != 0){
	printf("Machine size %d not commensurate with I/O size %d\n",
	       mesh->machsize[i], mesh->iomachsize[i]);
	return NULL;
      }
    }
  }

  return mesh;
}

void qio_destroy_topology(QIO_Mesh_Topology *mesh){
  if(mesh->machsize != NULL)
    free(mesh->machsize); 
  if(mesh->iomachsize != mesh->machsize && mesh->iomachsize != NULL){
    free(mesh->iomachsize);
  }
  free(mesh);  mesh = NULL;
}


