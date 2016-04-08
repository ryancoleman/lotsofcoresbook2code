#include "extensions.h"
#include <stdlib.h>

double distance(double *a, double *b);

// returns the squared distance between a 3d double vector
// and a 3d int vector
double distance3d2_di(double *a, int *b)
{
  double sum = 0;
  double diff;
  for (int c = 0; c < 3; c++) {
    diff = a[c] - (double)b[c];
    sum += diff*diff;
  }
  return sum;
} 

PyObject *exterior_electron_density_region(PyObject *self, PyObject *args)
{
  PyArrayObject* ai;
  PyArrayObject* aatom_c;
  PyArrayObject* beg_c;
  PyArrayObject* end_c;
  PyArrayObject* hh_c;
  PyArrayObject* vdWrad;
  if (!PyArg_ParseTuple(args, "OOOOOO", &ai, &aatom_c, 
			&beg_c, &end_c, &hh_c, &vdWrad))
    return NULL;

  long *aindex = LONGP(ai);
  int natoms = PyArray_DIM(aatom_c, 0);
  double *atom_c = DOUBLEP(aatom_c);
  long *beg = LONGP(beg_c);
  long *end = LONGP(end_c);
  double *h_c = DOUBLEP(hh_c);
  double *vdWradius = DOUBLEP(vdWrad);

  int n[3], ij;
  double pos[3];
  for (int c = 0; c < 3; c++) { n[c] = end[c] - beg[c]; }
  // loop over all points
  for (int i = 0; i < n[0]; i++) {
    pos[0] = (beg[0] + i) * h_c[0];
    for (int j = 0; j < n[1]; j++) {
      pos[1] = (beg[1] + j) * h_c[1];
      ij = (i*n[1] + j)*n[2];
      for (int k = 0; k < n[2]; k++) {
	pos[2] = (beg[2] + k) * h_c[2];
	aindex[ij + k] = (long) 1; /* assume outside the structure */
	// loop over all atoms
	for (int a=0; a < natoms; a++) {
	  double d = distance(atom_c + a*3, pos);
	  if (d < vdWradius[a]) {
	    aindex[ij + k] = (long) 0; /* this is inside */
	    a = natoms;
	  }
	}
      }
    }
  }

  Py_RETURN_NONE;
}

