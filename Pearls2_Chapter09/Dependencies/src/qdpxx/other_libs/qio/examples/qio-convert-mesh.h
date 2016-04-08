#ifndef QIO_CONVERT_MESH_H
#define QIO_CONVERT_MESH_H

#include <qio.h>

/* layout_hyper_mesh */
int setup_layout(int len[], int nd, int nsquares2[], int ndim2);
int node_number(const int x[]);
int node_index(const int x[]);
void get_coords(int x[], int node, int index);
int num_sites(int node);

typedef struct
{
  int machdim;
  int *machsize;
  int numnodes;
  int *iomachsize;
  int number_io_nodes;
} QIO_Mesh_Topology;

/* Lexicographic utilities */
void lex_init(int *dimp, int coords[], int dim);
int lex_next(int *dimp, int coords[], int dim, int size[]);
void lex_coords(int coords[], const int dim, const int size[], 
		const size_t rank);
size_t lex_rank(const int coords[], int dim, int size[]);
int *lex_allocate_coords(int dim, char *myname);

int qio_mesh_convert(QIO_Filesystem *fs, QIO_Mesh_Topology *mesh,
		     int argc, char *argv[]);

QIO_Mesh_Topology *qio_read_topology(int onetoone);
void qio_destroy_topology(QIO_Mesh_Topology *mesh);

#endif /* QIO_CONVERT_MESH_H */

