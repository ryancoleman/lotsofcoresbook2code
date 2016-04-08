/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include "spline.h"

static void spline_dealloc(SplineObject *xp)
{
  bmgs_deletespline(&xp->spline);
  PyObject_DEL(xp);
}

static PyObject * spline_get_cutoff(SplineObject *self, PyObject *args)
{
  return Py_BuildValue("d", self->spline.dr * self->spline.nbins);
}

static PyObject * spline_get_angular_momentum_number(SplineObject *self,
                                                     PyObject *args)
{
  return Py_BuildValue("i", self->spline.l);
}

static PyObject * spline_get_value_and_derivative(SplineObject *obj, 
                                                  PyObject *args,
                                                  PyObject *kwargs)
{
  double r;
  if (!PyArg_ParseTuple(args, "d", &r))
    return NULL;  
  double f;
  double dfdr;
  bmgs_get_value_and_derivative(&obj->spline, r, &f, &dfdr);
  return Py_BuildValue("(dd)", f, dfdr);
}


// Convert boundary point z-ranges to grid indices for the 2*l+1 boxes
static PyObject * spline_get_indices_from_zranges(SplineObject *self,
                                                      PyObject *args)
{
  PyArrayObject* beg_c_obj;
  PyArrayObject* end_c_obj;
  PyArrayObject* G_b_obj;
  int nm = 2 * self->spline.l + 1;

  if (!PyArg_ParseTuple(args, "OOO", &beg_c_obj, &end_c_obj, &G_b_obj))
    return NULL; 

  long* beg_c = LONGP(beg_c_obj);
  long* end_c = LONGP(end_c_obj);

  int ngmax = ((end_c[0] - beg_c[0]) *
               (end_c[1] - beg_c[1]) *
               (end_c[2] - beg_c[2]));

  int* G_B = INTP(G_b_obj);
  int nB = PyArray_DIMS(G_b_obj)[0];

  int ng = 0;
  for (int b = 0; b < nB; b+=2)
    ng += G_B[b+1]-G_B[b];

  npy_intp gm_dims[2] = {ng, nm};
  PyArrayObject* indices_gm_obj = (PyArrayObject*)PyArray_SimpleNew(2, gm_dims, 
                                                                    NPY_INT);

  int* p = INTP(indices_gm_obj);
  for (int b = 0; b < nB; b += 2) {
    int Ga = G_B[b], Gb = G_B[b+1];
    for (int G = Ga; G < Gb; G++)
      for (int m = 0; m < nm; m++)
        *p++ = m * ngmax + G;
    }

  // PyObjects created in the C code will be initialized with a refcount
  // of 1, for which reason we'll have to decref them when done here
  PyObject* values = Py_BuildValue("(Oii)", indices_gm_obj, ng, nm);
  Py_DECREF(indices_gm_obj);
  return values;
}


static PyMethodDef spline_methods[] = {
    {"get_cutoff",
     (PyCFunction)spline_get_cutoff, METH_VARARGS, 0},
    {"get_angular_momentum_number", 
     (PyCFunction)spline_get_angular_momentum_number, METH_VARARGS, 0},
    {"get_value_and_derivative", 
     (PyCFunction)spline_get_value_and_derivative, METH_VARARGS, 0},
    {"get_indices_from_zranges", 
     (PyCFunction)spline_get_indices_from_zranges, METH_VARARGS, 0},
    {NULL, NULL, 0, NULL}
};

static PyObject* spline_get_attr(PyObject *obj, char *name)
{
    return Py_FindMethod(spline_methods, obj, name);
}

static PyObject * spline_call(SplineObject *obj, PyObject *args,
                              PyObject *kwargs)
{
  double r;
  if (!PyArg_ParseTuple(args, "d", &r))
    return NULL;  
  return Py_BuildValue("d", bmgs_splinevalue(&obj->spline, r));
}

PyTypeObject SplineType = {
  PyObject_HEAD_INIT(NULL) 0,
  "Spline",
  sizeof(SplineObject), 0,
  (destructor)spline_dealloc, 0,
  spline_get_attr, 0, 0, 0, 0, 0, 0, 0,
  (ternaryfunc)spline_call
};

PyObject * NewSplineObject(PyObject *self, PyObject *args)
{
  int l;
  double rcut;
  PyArrayObject* farray;
  if (!PyArg_ParseTuple(args, "idO", &l, &rcut, &farray)) 
    return NULL;
  SplineObject *spline = PyObject_NEW(SplineObject, &SplineType);
  if (spline == NULL)
    return NULL;
  int nbins = PyArray_DIMS(farray)[0] - 1;
  double dr = rcut / nbins;
  spline->spline = bmgs_spline(l, dr, nbins, DOUBLEP(farray));
  return (PyObject*)spline;
}
