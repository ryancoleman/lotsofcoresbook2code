#ifdef GPAW_WITH_FFTW
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include <fftw3.h>


/* Create plan and return pointer to plan as a string */
PyObject * FFTWPlan(PyObject *self, PyObject *args)
{
    PyArrayObject* in;
    PyArrayObject* out;
    int sign;
    unsigned int flags;
    if (!PyArg_ParseTuple(args, "OOiI",
			  &in, &out, &sign, &flags))
        return NULL;

    fftw_plan* plan = (fftw_plan*)malloc(sizeof(fftw_plan));
    
    if (in->descr->type_num == PyArray_DOUBLE)
      *plan = fftw_plan_dft_r2c(in->nd, in->dimensions,
				(double*)in->data,
				(double (*)[2])out->data,
				flags);
    else if (out->descr->type_num == PyArray_DOUBLE)
        *plan = fftw_plan_dft_c2r(in->nd, out->dimensions,
				  (double (*)[2])in->data,
				  (double*)out->data,
				  flags);
    else
        *plan = fftw_plan_dft(in->nd, out->dimensions,
			      (double (*)[2])in->data,
			      (double (*)[2])out->data,
			      sign, flags);
    
    return Py_BuildValue("s#", plan, sizeof(fftw_plan*));
}


PyObject * FFTWExecute(PyObject *self, PyObject *args)
{
    fftw_plan* plan;
    int n;
    if (!PyArg_ParseTuple(args, "s#", &plan, &n))
        return NULL;
    fftw_execute(*plan);
    Py_RETURN_NONE;
}


PyObject * FFTWDestroy(PyObject *self, PyObject *args)
{
    fftw_plan* plan;
    int n;
    if (!PyArg_ParseTuple(args, "s#", &plan, &n))
        return NULL;
    fftw_destroy_plan(*plan);
    free(plan);
    Py_RETURN_NONE;
}

#endif // GPAW_WITH_FFTW
