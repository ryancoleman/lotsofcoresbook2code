#include "extensions.h"
#include <stdlib.h>

PyObject *plane_wave_grid(PyObject *self, PyObject *args)
{
  PyArrayObject* beg_c;
  PyArrayObject* end_c;
  PyArrayObject* h_c;
  PyArrayObject* k_c;
  PyArrayObject* r0_c;
  PyArrayObject* pw_g;
  if (!PyArg_ParseTuple(args, "OOOOOO", &beg_c, &end_c, &h_c, 
			&k_c, &r0_c, &pw_g))
    return NULL;

  long *beg = LONGP(beg_c);
  long *end = LONGP(end_c);
  double *h = DOUBLEP(h_c);
  double *vk = DOUBLEP(k_c);
  double *vr0 = DOUBLEP(r0_c);
  double_complex *pw = COMPLEXP(pw_g);

  double kr[3], kr0[3];
  int n[3], ij;
  for (int c = 0; c < 3; c++) { 
    n[c] = end[c] - beg[c];
    kr0[c] = vk[c] * vr0[c];
  }
  for (int i = 0; i < n[0]; i++) {
    kr[0] = vk[0] * h[0] * (beg[0] + i) - kr0[0];
    for (int j = 0; j < n[1]; j++) {
      kr[1] = kr[0] + vk[1] * h[1] * (beg[1] + j) - kr0[1];
      ij = (i*n[1] + j)*n[2];
      for (int k = 0; k < n[2]; k++) {
	kr[2] = kr[1] + vk[2] * h[2] * (beg[2] + k) - kr0[2];
	pw[ij + k] = cos(kr[2]) + I * sin(kr[2]);
      }
    }
  }
  Py_RETURN_NONE;
}


