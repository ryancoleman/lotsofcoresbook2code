/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Copyright (C) 2005-2009  CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#include <pthread.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "extensions.h"
#include "bc.h"
#include "mympi.h"
#include "bmgs/bmgs.h"

#ifdef GPAW_ASYNC
  #define GPAW_ASYNC_D 3
#else
  #define GPAW_ASYNC_D 1
#endif

typedef struct
{
  PyObject_HEAD
  boundary_conditions* bc;
  int p;
  int k;
  bool interpolate;
  MPI_Request recvreq[2];
  MPI_Request sendreq[2];
  int skip[3][2];
  int size_out[3];          /* Size of the output grid */
} TransformerObject;

static void Transformer_dealloc(TransformerObject *self)
{
  free(self->bc);
  PyObject_DEL(self);
}

struct transapply_args{
  int thread_id;
  TransformerObject *self;
  int ng;
  int ng2;
  int nin;
  int nthds;
  const double* in;
  double* out;
  int real;
  const double_complex* ph;
};

void *transapply_worker(void *threadarg)
{
  struct transapply_args *args = (struct transapply_args *) threadarg;
  boundary_conditions* bc = args->self->bc;
  TransformerObject *self = args->self;
  double* sendbuf = GPAW_MALLOC(double, bc->maxsend * GPAW_ASYNC_D);
  double* recvbuf = GPAW_MALLOC(double, bc->maxrecv * GPAW_ASYNC_D);
  double* buf = GPAW_MALLOC(double, args->ng2);
  int buf2size = args->ng2;
  if (self->interpolate)
    buf2size *= 16;
  else
    buf2size /= 2;
  double* buf2 = GPAW_MALLOC(double, buf2size);
  MPI_Request recvreq[2 * GPAW_ASYNC_D];
  MPI_Request sendreq[2 * GPAW_ASYNC_D];

  int chunksize = args->nin / args->nthds;
  if (!chunksize)
    chunksize = 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= args->nin)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > args->nin)
    nend = args->nin;

  int out_ng = bc->ndouble * self->size_out[0] * self->size_out[1]
               * self->size_out[2];

  for (int n = nstart; n < nend; n++)
    {
      const double* in = args->in + n * args->ng;
      double* out = args->out + n * out_ng;
      for (int i = 0; i < 3; i++)
        {
          bc_unpack1(bc, in, buf, i,
                     recvreq, sendreq,
                     recvbuf, sendbuf, args->ph + 2 * i,
                     args->thread_id, 1);
          bc_unpack2(bc, buf, i,
                     recvreq, sendreq, recvbuf, 1);
        }
      if (args->real)
        {
          if (self->interpolate)
            bmgs_interpolate(self->k, self->skip, buf, bc->size2,
                             out, buf2);
          else
            bmgs_restrict(self->k, buf, bc->size2,
                          out, buf2);
        }
      else
        {
          if (self->interpolate)
            bmgs_interpolatez(self->k, self->skip, (double_complex*)buf,
                              bc->size2, (double_complex*)out,
                              (double_complex*) buf2);
          else
            bmgs_restrictz(self->k, (double_complex*) buf,
                           bc->size2, (double_complex*)out,
                           (double_complex*) buf2);
        }
    }
  free(buf2);
  free(buf);
  free(recvbuf);
  free(sendbuf);
  return NULL;
}



static PyObject* Transformer_apply(TransformerObject *self, PyObject *args)
{
  PyArrayObject* input;
  PyArrayObject* output;
  PyArrayObject* phases = 0;
  if (!PyArg_ParseTuple(args, "OO|O", &input, &output, &phases))
    return NULL;

  int nin = 1;
  if (PyArray_NDIM(input) == 4)
    nin = PyArray_DIMS(input)[0];

  boundary_conditions* bc = self->bc;
  const int* size1 = bc->size1;
  const int* size2 = bc->size2;
  int ng = bc->ndouble * size1[0] * size1[1] * size1[2];
  int ng2 = bc->ndouble * size2[0] * size2[1] * size2[2];

  const double* in = DOUBLEP(input);
  double* out = DOUBLEP(output);
  bool real = (PyArray_DESCR(input)->type_num == NPY_DOUBLE);
  const double_complex* ph = (real ? 0 : COMPLEXP(phases));

  int nthds = 1;
#ifdef GPAW_OMP
  if (getenv("OMP_NUM_THREADS") != NULL)
    nthds = atoi(getenv("OMP_NUM_THREADS"));
#endif
  struct transapply_args *wargs = GPAW_MALLOC(struct transapply_args, nthds);
  pthread_t *thds = GPAW_MALLOC(pthread_t, nthds);

  for(int i=0; i < nthds; i++)
    {
      (wargs+i)->thread_id = i;
      (wargs+i)->nthds = nthds;
      (wargs+i)->self = self;
      (wargs+i)->ng = ng;
      (wargs+i)->ng2 = ng2;
      (wargs+i)->nin = nin;
      (wargs+i)->in = in;
      (wargs+i)->out = out;
      (wargs+i)->real = real;
      (wargs+i)->ph = ph;
    }

#ifdef GPAW_OMP
  for(int i=1; i < nthds; i++)
    pthread_create(thds + i, NULL, transapply_worker, (void*) (wargs+i));
#endif
  transapply_worker(wargs);
#ifdef GPAW_OMP
  for(int i=1; i < nthds; i++)
    pthread_join(*(thds+i), NULL);
#endif
  free(wargs);
  free(thds);

  Py_RETURN_NONE;
}

static PyObject * Transformer_get_async_sizes(TransformerObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

#ifdef GPAW_ASYNC
  return Py_BuildValue("(ii)", 1, GPAW_ASYNC_D);
#else
  return Py_BuildValue("(ii)", 0, GPAW_ASYNC_D);
#endif
}

static PyMethodDef Transformer_Methods[] = {
    {"apply", (PyCFunction)Transformer_apply, METH_VARARGS, NULL},
    {"get_async_sizes",
     (PyCFunction)Transformer_get_async_sizes, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static PyObject* Transformer_getattr(PyObject *obj, char *name)
{
    return Py_FindMethod(Transformer_Methods, obj, name);
}

PyTypeObject TransformerType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "Transformer",
  sizeof(TransformerObject),
  0,
  (destructor)Transformer_dealloc,
  0,
  Transformer_getattr
};

PyObject * NewTransformerObject(PyObject *obj, PyObject *args)
{
  PyArrayObject* size_in;
  PyArrayObject* size_out;
  int k;
  PyArrayObject* paddings;
  PyArrayObject* npaddings;
  PyArrayObject* skip;
  PyArrayObject* neighbors;
  int real;
  PyObject* comm_obj;
  int interpolate;
  if (!PyArg_ParseTuple(args, "OOiOOOOiOi",
                        &size_in, &size_out, &k, &paddings, &npaddings, &skip,
                        &neighbors, &real, &comm_obj,
                        &interpolate))
    return NULL;

  TransformerObject* self = PyObject_NEW(TransformerObject, &TransformerType);
  if (self == NULL)
    return NULL;

  self->k = k;
  self->interpolate = interpolate;

  MPI_Comm comm = MPI_COMM_NULL;
  if (comm_obj != Py_None)
    comm = ((MPIObject*)comm_obj)->comm;

  const long (*nb)[2] = (const long (*)[2])LONGP(neighbors);
  const long (*pad)[2] = (const long (*)[2])LONGP(paddings);
  const long (*npad)[2] = (const long (*)[2])LONGP(npaddings);
  const long (*skp)[2] = (const long (*)[2])LONGP(skip);
  self->bc = bc_init(LONGP(size_in), pad, npad, nb, comm, real, 0);

  for (int c = 0; c < 3; c++)
      self->size_out[c] = LONGP(size_out)[c];

  for (int c = 0; c < 3; c++)
    for (int d = 0; d < 2; d++)
      self->skip[c][d] = (int)skp[c][d];

  return (PyObject*)self;
}
