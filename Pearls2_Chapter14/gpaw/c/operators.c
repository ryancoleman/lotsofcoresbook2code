#define USE_LEO    1
#define USE_TIMING 0

/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Copyright (C) 2005       CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

//*** The apply operator and some associate structors are imple-  ***//
//*** mented in two version: a original version and a speciel     ***//
//*** OpenMP version. By default the original version will        ***//
//*** be used, but it's possible to use the OpenMP version        ***//
//*** by compiling gpaw with the macro GPAW_OMP defined and       ***//
//*** and the compile/link option "-fopenmp".                     ***//
//*** Author of the optimized OpenMP code:                        ***//
//*** Mads R. B. Kristensen - madsbk@diku.dk                      ***//

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <alloca.h>
#include <pthread.h>
#include "extensions.h"
#include "bc.h"
#include "mympi.h"

#if USE_TIMING
#include <omp.h>
#endif 

#ifdef GPAW_ASYNC
  #define GPAW_ASYNC3 3
  #define GPAW_ASYNC2 2
#else
  #define GPAW_ASYNC3 1
  #define GPAW_ASYNC2 1
#endif

extern int gpaw_offload_enabled;

typedef struct
{
  PyObject_HEAD
  bmgsstencil stencil;
  boundary_conditions* bc;
  MPI_Request recvreq[2];
  MPI_Request sendreq[2];
} OperatorObject;

static void Operator_dealloc(OperatorObject *self)
{
  free(self->bc);
  PyObject_DEL(self);
}


static PyObject * Operator_relax(OperatorObject *self,
                                 PyObject *args)
{
#if USE_TIMING
  double tstart, tend;
#endif
  int relax_method;
  PyArrayObject* func;
  PyArrayObject* source;
  int nrelax;
  double w = 1.0;
  int stencil_size;
  unsigned char *stencil_buffer;
  bmgsstencil *s;
  /* TODO: offload_enabled has been introduced as a hack for now */
  int offload_enabled = 1;
  if (!PyArg_ParseTuple(args, "iOOi|di", &relax_method, &func, &source,
                        &nrelax, &w, &offload_enabled))
    return NULL;

  const boundary_conditions* bc = self->bc;

  double* fun = DOUBLEP(func);
  const double* src = DOUBLEP(source);
  const double_complex* ph;

  const int* size2 = bc->size2;
  double* buf = GPAW_MALLOC(double, size2[0] * size2[1] * size2[2] *
                            bc->ndouble);
  double* sendbuf = GPAW_MALLOC(double, bc->maxsend);
  double* recvbuf = GPAW_MALLOC(double, bc->maxrecv);

  ph = 0;

  /* convert the stencil structure to a flat array */
  s = &self->stencil;
  stencil_size = sizeof(bmgsstencil) + s->ncoefs * sizeof(double) + s->ncoefs * sizeof(long);
  stencil_buffer = (unsigned char*) alloca(stencil_size);
  *(bmgsstencil*) stencil_buffer = *s;
  memcpy(stencil_buffer + sizeof(bmgsstencil), s->coefs, s->ncoefs * sizeof(double));
  memcpy(stencil_buffer + sizeof(bmgsstencil) + s->ncoefs * sizeof(double), s->offsets, s->ncoefs * sizeof(long));

#if 0
  if (gpaw_offload_enabled && !offload_enabled) 
    printf("skipping offload\n");
  if (gpaw_offload_enabled && offload_enabled)
    printf("OFFLOADING!!!\n");
#endif  
  
#if USE_LEO
  /* pre-allocate the buffers */
#pragma offload_transfer target(mic:0) if(gpaw_offload_enabled && offload_enabled) \
                                       in(stencil_buffer:length(stencil_size) alloc_if(1) free_if(0)) \
                                       nocopy(buf:length(size2[0] * size2[1] * size2[2]) alloc_if(1) free_if(0))
                                       /* We assume that fun (func) and src (source) have been allocated through
                                          offload_device.bind() and thus already live on the traget device.  Thus,
                                          we do not need these clauses:
                                            nocopy(fun:length(s->n[0] * s->n[1] * s->n[2]) alloc_if(1) free_if(0)) 
                                            in(src:length(s->n[0] * s->n[1] * s->n[2]) alloc_if(1) free_if(0))
                                        */
#endif  
  for (int n = 0; n < nrelax; n++ )
    {
      for (int i = 0; i < 3; i++)
        {
          bc_unpack1(bc, fun, buf, i,
               self->recvreq, self->sendreq,
               recvbuf, sendbuf, ph + 2 * i, 0, 1);
          bc_unpack2(bc, buf, i,
               self->recvreq, self->sendreq, recvbuf, 1);
        }
#if USE_LEO
#if USE_TIMING
        tstart = omp_get_wtime();
#endif
        /* perform the actual offload, only pass in new data, keep buffers */
#pragma offload target(mic:0) if(gpaw_offload_enabled && offload_enabled) \
                              in(stencil_buffer:length(0) alloc_if(0) free_if(0)) \
                              inout(fun:length(s->n[0] * s->n[1] * s->n[2]) alloc_if(0) free_if(0)) \
                              in(buf:length(size2[0] * size2[1] * size2[2]) alloc_if(0) free_if(0)) \
                              in(src:length(0) alloc_if(0) free_if(0))
#endif                              
      {
#if USE_TIMING
          double mic_tstart, mic_tend;
          mic_tstart = omp_get_wtime();
#endif
          /* deserialize the stencil */
          bmgsstencil *stencil;
          stencil = (bmgsstencil*) stencil_buffer;
          stencil->coefs = (double*) (stencil_buffer + sizeof(bmgsstencil));
          stencil->offsets = (long*) (stencil_buffer + sizeof(bmgsstencil) + stencil->ncoefs * sizeof(double));
          bmgs_relax(relax_method, stencil, buf, fun, src, w);
#if USE_TIMING
          mic_tend = omp_get_wtime();
          printf("MIC time: %lf\n", mic_tend - mic_tstart);
#endif
      }
#if USE_TIMING
      tend = omp_get_wtime();
      printf("HOST time: %lf\n", tend - tstart);
#endif      
    }
#if USE_LEO
  /* clean up the device buffers */
#pragma offload_transfer target(mic:0) if(gpaw_offload_enabled && offload_enabled) \
                                       nocopy(stencil_buffer:length(stencil_size) alloc_if(0) free_if(1)) \
                                       nocopy(buf:length(size2[0] * size2[1] * size2[2]) alloc_if(0) free_if(1)) 
                                       /* We did not pre-allocate fun and src (see above).  Thus, we will keep them
                                          as they are, avoid that the we interfere with the Python code here.
                                            nocopy(fun:length(s->n[0] * s->n[1] * s->n[2]) alloc_if(0) free_if(1)) 
                                            nocopy(src:length(s->n[0] * s->n[1] * s->n[2]) alloc_if(0) free_if(1))
                                        */
#endif  
  free(recvbuf);
  free(sendbuf);
  free(buf);
  Py_RETURN_NONE;
}

struct apply_args{
  int thread_id;
  OperatorObject *self;
  int ng;
  int ng2;
  int nin;
  int nthds;
  int chunksize;
  int chunkinc;
  const double* in;
  double* out;
  int real;
  const double_complex* ph;
};

//Plain worker
void *apply_worker(void *threadarg)
{
  struct apply_args *args = (struct apply_args *) threadarg;
  boundary_conditions* bc = args->self->bc;
  MPI_Request recvreq[2];
  MPI_Request sendreq[2];

  int chunksize = args->nin / args->nthds;
  if (!chunksize)
    chunksize = 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= args->nin)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > args->nin)
    nend = args->nin;
  if (chunksize > args->chunksize)
    chunksize = args->chunksize;

  double* sendbuf = GPAW_MALLOC(double, bc->maxsend * args->chunksize);
  double* recvbuf = GPAW_MALLOC(double, bc->maxrecv * args->chunksize);
  double* buf = GPAW_MALLOC(double, args->ng2 * args->chunksize);

  for (int n = nstart; n < nend; n += chunksize)
    {
      if (n + chunksize >= nend && chunksize > 1)
        chunksize = nend - n;
      const double* in = args->in + n * args->ng;
      double* out = args->out + n * args->ng;
      for (int i = 0; i < 3; i++)
        {
          bc_unpack1(bc, in, buf, i,
                     recvreq, sendreq,
                     recvbuf, sendbuf, args->ph + 2 * i,
                     args->thread_id, chunksize);
          bc_unpack2(bc, buf, i, recvreq, sendreq, recvbuf, chunksize);
        }
      for (int m = 0; m < chunksize; m++)
        if (args->real)
          bmgs_fd(&args->self->stencil, buf + m * args->ng2, out + m * args->ng);
        else
          bmgs_fdz(&args->self->stencil, (const double_complex*) (buf + m * args->ng2),
                                         (double_complex*) (out + m * args->ng));
    }
  free(buf);
  free(recvbuf);
  free(sendbuf);
  return NULL;
}

//Async worker
void *apply_worker_cfd_async(void *threadarg)
{
  struct apply_args *args = (struct apply_args *) threadarg;
  boundary_conditions* bc = args->self->bc;
  MPI_Request recvreq[2 * GPAW_ASYNC3];
  MPI_Request sendreq[2 * GPAW_ASYNC3];

  int chunksize = args->nin / args->nthds;
  if (!chunksize)
    chunksize = 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= args->nin)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > args->nin)
    nend = args->nin;
  if (chunksize > args->chunksize)
    chunksize = args->chunksize;

  double* sendbuf = GPAW_MALLOC(double, bc->maxsend * GPAW_ASYNC3 *
                                args->chunksize);
  double* recvbuf = GPAW_MALLOC(double, bc->maxrecv * GPAW_ASYNC3 *
                                args->chunksize);
  double* buf = GPAW_MALLOC(double, args->ng2 * args->chunksize);

  for (int n = nstart; n < nend; n += chunksize)
    {
      if (n + chunksize >= nend && chunksize > 1)
        chunksize = nend - n;
      const double* in = args->in + n * args->ng;
      double* out = args->out + n * args->ng;
      for (int i = 0; i < 3; i++)
        {
          bc_unpack1(bc, in, buf, i,
                     recvreq + i * 2, sendreq + i * 2,
                     recvbuf + i * bc->maxrecv * chunksize,
                     sendbuf + i * bc->maxsend * chunksize, args->ph + 2 * i,
                     args->thread_id, chunksize);
        }
      for (int i = 0; i < 3; i++)
        {
          bc_unpack2(bc, buf, i,
                     recvreq + i * 2, sendreq + i * 2,
                     recvbuf + i * bc->maxrecv * chunksize, chunksize);
        }
      for (int m = 0; m < chunksize; m++)
        if (args->real)
          bmgs_fd(&args->self->stencil, buf + m * args->ng2, out + m * args->ng);
        else
          bmgs_fdz(&args->self->stencil, (const double_complex*) (buf + m * args->ng2),
                                         (double_complex*) (out + m * args->ng));
    }
  free(buf);
  free(recvbuf);
  free(sendbuf);
  return NULL;
}

//Double buffering async worker
void *apply_worker_cfd(void *threadarg)
{
  struct apply_args *args = (struct apply_args *) threadarg;
  boundary_conditions* bc = args->self->bc;
  MPI_Request recvreq[2 * GPAW_ASYNC3 * GPAW_ASYNC2];
  MPI_Request sendreq[2 * GPAW_ASYNC3 * GPAW_ASYNC2];

  int chunksize = args->nin / args->nthds;
  if (!chunksize)
    chunksize = 1;
  int nstart = args->thread_id * chunksize;
  if (nstart >= args->nin)
    return NULL;
  int nend = nstart + chunksize;
  if (nend > args->nin)
    nend = args->nin;
  if (chunksize > args->chunksize)
    chunksize = args->chunksize;

  int chunk = args->chunkinc;
  if (chunk > chunksize);
    chunk = chunksize;

  double* sendbuf = GPAW_MALLOC(double, bc->maxsend * args->chunksize
                                * GPAW_ASYNC3 * GPAW_ASYNC2);
  double* recvbuf = GPAW_MALLOC(double, bc->maxrecv * args->chunksize
                                * GPAW_ASYNC3 * GPAW_ASYNC2);
  double* buf = GPAW_MALLOC(double, args->ng2 * args->chunksize * GPAW_ASYNC2);

  int odd = 0;
  const double* in = args->in + nstart * args->ng;
  double* out;
  for (int i = 0; i < 3; i++)
    bc_unpack1(bc, in, buf + odd * args->ng2 * chunksize, i,
               recvreq + odd * 2 + i * 4, sendreq + odd * 2 + i * 4,
               recvbuf + odd * bc->maxrecv * chunksize + i * bc->maxrecv * chunksize * GPAW_ASYNC2,
               sendbuf + odd * bc->maxsend * chunksize + i * bc->maxsend * chunksize * GPAW_ASYNC2, args->ph + 2 * i,
               args->thread_id, chunk);
  odd = odd ^ 1;
  int last_chunk = chunk;
  for (int n = nstart+chunk; n < nend; n += chunk)
    {
      last_chunk += args->chunkinc;
      if (last_chunk > chunksize);
        last_chunk = chunksize;

      if (n + last_chunk >= nend && last_chunk > 1)
        last_chunk = nend - n;
      in = args->in + n * args->ng;
      out = args->out + (n-chunk) * args->ng;
      for (int i = 0; i < 3; i++)
        {
          bc_unpack1(bc, in, buf + odd * args->ng2 * chunksize, i,
                     recvreq + odd * 2 + i * 4, sendreq + odd * 2 + i * 4,
                     recvbuf + odd * bc->maxrecv * chunksize + i * bc->maxrecv * chunksize * GPAW_ASYNC2,
                     sendbuf + odd * bc->maxsend * chunksize + i * bc->maxsend * chunksize * GPAW_ASYNC2, args->ph + 2 * i,
                     args->thread_id, last_chunk);
        }
      odd = odd ^ 1;
      for (int i = 0; i < 3; i++)
        {
          bc_unpack2(bc, buf + odd * args->ng2 * chunksize, i,
                     recvreq + odd * 2 + i * 4, sendreq + odd * 2 + i * 4,
                     recvbuf + odd * bc->maxrecv * chunksize + i * bc->maxrecv * chunksize * GPAW_ASYNC2, chunk);
        }
      for (int m = 0; m < chunk; m++)
        if (args->real)
          bmgs_fd(&args->self->stencil, buf + m * args->ng2 + odd * args->ng2 * chunksize,
                                        out + m * args->ng);
        else
          bmgs_fdz(&args->self->stencil, (const double_complex*) (buf + m * args->ng2 + odd * args->ng2 * chunksize),
                                         (double_complex*) (out + m * args->ng));
      chunk = last_chunk;
    }

  odd = odd ^ 1;
  out = args->out + (nend-last_chunk) * args->ng;
  for (int i = 0; i < 3; i++)
    {
      bc_unpack2(bc, buf + odd * args->ng2 * chunksize, i,
                 recvreq + odd * 2 + i * 4, sendreq + odd * 2 + i * 4,
                 recvbuf + odd * bc->maxrecv * chunksize + i * bc->maxrecv * chunksize * GPAW_ASYNC2, last_chunk);
    }
  for (int m = 0; m < last_chunk; m++)
    if (args->real)
      bmgs_fd(&args->self->stencil, buf + m * args->ng2 + odd * args->ng2 * chunksize,
                                    out + m * args->ng);
    else
      bmgs_fdz(&args->self->stencil, (const double_complex*) (buf + m * args->ng2 + odd * args->ng2 * chunksize),
                                     (double_complex*) (out + m * args->ng));

  free(buf);
  free(recvbuf);
  free(sendbuf);
  return NULL;
}

static PyObject * Operator_apply(OperatorObject *self,
                                 PyObject *args)
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
  const double_complex* ph;

  bool real = (PyArray_DESCR(input)->type_num == NPY_DOUBLE);

  if (real)
    ph = 0;
  else
    ph = COMPLEXP(phases);

  int chunksize = 1;
  if (getenv("GPAW_CHUNK_SIZE") != NULL)
    chunksize = atoi(getenv("GPAW_CHUNK_SIZE"));

  int chunkinc = chunksize;
  if (getenv("GPAW_CHUNK_INC") != NULL)
    chunkinc = atoi(getenv("GPAW_CHUNK_INC"));

  int nthds = 1;
#ifdef GPAW_OMP
  if (getenv("OMP_NUM_THREADS") != NULL)
    nthds = atoi(getenv("OMP_NUM_THREADS"));
#endif
  struct apply_args *wargs = GPAW_MALLOC(struct apply_args, nthds);
  pthread_t *thds = GPAW_MALLOC(pthread_t, nthds);

  for(int i=0; i < nthds; i++)
    {
      (wargs+i)->thread_id = i;
      (wargs+i)->nthds = nthds;
      (wargs+i)->chunksize = chunksize;
      (wargs+i)->chunkinc = chunkinc;
      (wargs+i)->self = self;
      (wargs+i)->ng = ng;
      (wargs+i)->ng2 = ng2;
      (wargs+i)->nin = nin;
      (wargs+i)->in = in;
      (wargs+i)->out = out;
      (wargs+i)->real = real;
      (wargs+i)->ph = ph;
    }
#ifndef GPAW_ASYNC
  if (1)
#else
  if (bc->cfd == 0)
#endif
    {
    #ifdef GPAW_OMP
      for(int i=1; i < nthds; i++)
        pthread_create(thds + i, NULL, apply_worker, (void*) (wargs+i));
    #endif
      apply_worker(wargs);
    }
  else
    {
    #ifdef GPAW_OMP
      for(int i=1; i < nthds; i++)
        pthread_create(thds + i, NULL, apply_worker_cfd, (void*) (wargs+i));
    #endif
      apply_worker_cfd(wargs);
    }
#ifdef GPAW_OMP
  for(int i=1; i < nthds; i++)
    pthread_join(*(thds+i), NULL);
#endif
  free(wargs);
  free(thds);

  Py_RETURN_NONE;
}


static PyObject * Operator_get_diagonal_element(OperatorObject *self,
                                              PyObject *args)
{
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  const bmgsstencil* s = &self->stencil;
  double d = 0.0;
  for (int n = 0; n < s->ncoefs; n++)
    if (s->offsets[n] == 0)
      d = s->coefs[n];

  return Py_BuildValue("d", d);
}

static PyObject * Operator_get_async_sizes(OperatorObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

#ifdef GPAW_ASYNC
  return Py_BuildValue("(iii)", 1, GPAW_ASYNC2, GPAW_ASYNC3);
#else
  return Py_BuildValue("(iii)", 0, GPAW_ASYNC2, GPAW_ASYNC3);
#endif
}

#ifndef __MIC__
static PyMethodDef Operator_Methods[] = {
    {"apply",
     (PyCFunction)Operator_apply, METH_VARARGS, NULL},
    {"relax",
     (PyCFunction)Operator_relax, METH_VARARGS, NULL},
    {"get_diagonal_element",
     (PyCFunction)Operator_get_diagonal_element, METH_VARARGS, NULL},
    {"get_async_sizes",
     (PyCFunction)Operator_get_async_sizes, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}

};

static PyObject* Operator_getattr(PyObject *obj, char *name)
{
    return Py_FindMethod(Operator_Methods, obj, name);
}

PyTypeObject OperatorType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "Operator",
  sizeof(OperatorObject),
  0,
  (destructor)Operator_dealloc,
  0,
  Operator_getattr
};

PyObject * NewOperatorObject(PyObject *obj, PyObject *args)
{
  PyArrayObject* coefs;
  PyArrayObject* offsets;
  PyArrayObject* size;
  int range;
  PyArrayObject* neighbors;
  int real;
  PyObject* comm_obj;
  int cfd;
  if (!PyArg_ParseTuple(args, "OOOiOiOi",
                        &coefs, &offsets, &size, &range, &neighbors,
                        &real, &comm_obj, &cfd))
    return NULL;

  OperatorObject *self = PyObject_NEW(OperatorObject, &OperatorType);
  if (self == NULL)
    return NULL;

  self->stencil = bmgs_stencil(PyArray_DIMS(coefs)[0], DOUBLEP(coefs),
                               LONGP(offsets), range, LONGP(size));

  const long (*nb)[2] = (const long (*)[2])LONGP(neighbors);
  const long padding[3][2] = {{range, range},
                             {range, range},
                             {range, range}};

  MPI_Comm comm = MPI_COMM_NULL;
  if (comm_obj != Py_None)
    comm = ((MPIObject*)comm_obj)->comm;

  self->bc = bc_init(LONGP(size), padding, padding, nb, comm, real, cfd);
  return (PyObject*)self;
}
#endif
