/******** layout_hyper_mesh.c *********/
/* adapted from MIMD version 6 and SciDAC API */
/* C. DeTar Jan 30, 2005 */

/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES 
   ON A MESH ARCHITECTURE */

/* On a mesh architecture the layout is fixed by the machine topology
   (or machine dimension) */

/* The API:

   setup_layout()  sets up layout
   node_number()   returns the node number on which a site lives
   node_index()    returns the index of the site on the node
   get_coords()    gives lattice coords from node & index
*/

#include <stdlib.h>
#include <stdio.h>

static int sites_on_node;
static int *squaresize;   /* dimensions of hypercubes */
static int *nsquares;     /* number of hypercubes in each direction */
static int ndim;
static int *size1[2], *size2;

int
setup_layout(
  int len[],         /* Lattice size */
  int nd,            /* Number of lattice dimensions */
  int nsquares2[],   /* Machine size: require len[i] % nsquares2[i] == 0 */
  int ndim2)         /* Number of machine dimensions: require ndim2 <= nd */
{
  int i;

  ndim = nd;
  squaresize = (int *) malloc(ndim*sizeof(int));
  nsquares = (int *) malloc(ndim*sizeof(int));

  for(i=0; i<ndim; ++i) {
    squaresize[i] = len[i];
    nsquares[i] = 1;
  }

  for(i=0; i<ndim; i++) {
    if(i<ndim2) nsquares[i] = nsquares2[i];
    else nsquares[i] = 1;
  }
  for(i=0; i<ndim; i++) {
    if(len[i]%nsquares[i] != 0) {
      printf("LATTICE SIZE DOESN'T FIT GRID\n");
      return 1;
    }
    squaresize[i] = len[i]/nsquares[i];
  }

  sites_on_node = 1;
  for(i=0; i<ndim; ++i) {
    sites_on_node *= squaresize[i];
  }

  size1[0] = (int *)malloc(2*(ndim+1)*sizeof(int));
  size1[1] = size1[0] + ndim + 1;
  size2 = (int *)malloc((ndim+1)*sizeof(int));

  size1[0][0] = 1;
  size1[1][0] = 0;
  size2[0] = 1;
  for(i=1; i<=ndim; i++) {
    size1[0][i] = size2[i-1]*(squaresize[i-1]/2)
                + size1[0][i-1]*(squaresize[i-1]%2);
    size1[1][i] = size2[i-1]*(squaresize[i-1]/2)
                + size1[1][i-1]*(squaresize[i-1]%2);
    size2[i] = size1[0][i] + size1[1][i];
    /* printf("%i\t%i\t%i\n", size1[0][i], size1[1][i], size2[i]); */
  }
  return 0;
}

int
node_number(const int x[])
{
  int i; 
  int *m;
  int rank;

  m = (int *)malloc(sizeof(int)*ndim);
  if( m==(int *)NULL ) {
    printf("Error allocating m in node_number\n");
    exit(1);
  }  

  /* Get mesh node or "logical" coordinates */
  for(i=0; i<ndim; i++) {
    m[i] = x[i]/squaresize[i];
  }

  /* This algorithm must match that of QMP_get_node_number_from() */
  rank = 0;
  for (i = ndim - 1; i >=0; i--)
    rank = rank * nsquares[i] + m[i];

  free(m);
  return rank;
}

int
node_index(const int x[])
{
  int i, r=0, p=0;

  for(i=ndim-1; i>=0; --i) {
    r = r*squaresize[i] + (x[i]%squaresize[i]);
    p += x[i];
  }

  if( p%2==0 ) { /* even site */
    r /= 2;
  } else {
    r = (r+sites_on_node)/2;
  }
  return r;
}

void
get_coords(int x[], int node, int index)
{
  int i, s, si;
  int *m = (int *)malloc(ndim*sizeof(int));
  int pos;

  si = index;

  /* The algorithm must match that of
     QMP_get_logical_coordinates_from() */
  pos = node;
  for (i = 0; i < ndim; i++) {
    m[i] = pos % nsquares[i];
    pos = pos / nsquares[i];
  }

  s = 0;
  /* p = node; */
  for(i=0; i<ndim; ++i) {
    /* x[i] = (p%nsquares[i])*squaresize[i]; */
    x[i] = m[i] * squaresize[i];
    s += x[i];
    /* p /= nsquares[i]; */
  }
  s &= 1;

  if(index>=size1[s][ndim]) {
    index -= size1[s][ndim];
    s ^= 1;
  }

  for(i=ndim-1; i>0; i--) {
    x[i] += 2*(index/size2[i]);
    index %= size2[i];
    if(index>=size1[s][i]) {
      index -= size1[s][i];
      s ^= 1;
      x[i]++;
    }
  }
  x[0] += 2*index + s;

  free(m);

  if(node_index(x)!=si) {
    fprintf(stderr,"layout_hyper_mesh: error in layout!\n");
    for(i=0; i<ndim; i++) {
      fprintf(stderr,"%i\t%i\t%i\n", size1[0][i], size1[1][i], size2[i]);
    }
    fprintf(stderr,"%i\t%i", node, si);
    for(i=0; i<ndim; i++) fprintf(stderr,"\t%i", x[i]);
    fprintf(stderr,"\n");
    exit(1);
  }
}


/* The number of sites on the specified node */
int num_sites(int node){
  return sites_on_node;
}
