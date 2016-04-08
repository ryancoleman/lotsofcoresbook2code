.. _c_extension:

============
C extensions
============

The GPAW Python code makes use of some compiled C code in the
dynamically linked extension ``_gpaw.so``.  In the following it is
demonstrated how a C function is made available to Python through a
Python extension (more details can be found in the official `Python
documentation`_.  The wrapper code from ``c/blas.c`` shows how to wrap
the two BLAS functions ``daxpy`` and ``zaxpy`` in Python::
  
 PyObject* axpy(PyObject *self, PyObject *args)
 {
   PyObject* alpha;
   PyArrayObject* x;
   PyArrayObject* y;
   if (!PyArg_ParseTuple(args, "OOO", &alpha, &x, &y)) 
     return NULL;
   integer n = x->dimensions[0];
   for (int d = 1; d < x->nd; d++)
     n *= x->dimensions[d];
   integer incx = 1;
   integer incy = 1;
   if (PyFloat_Check(alpha))
     {
       PyFloatObject* palpha = (PyFloatObject*)alpha;
       daxpy_(&n, &(palpha->ob_fval), 
             DOUBLEP(x), &incx,
             DOUBLEP(y), &incy);
     }
   else
     {
       PyComplexObject* palpha = (PyComplexObject*)alpha;
       zaxpy_(&n, (doublecomplex*)(&(palpha->cval)), 
              (doublecomplex*)COMPLEXP(x), &incx,
              (doublecomplex*)COMPLEXP(y), &incy);
     }
   Py_RETURN_NONE;
 }

In ``c/_gpaw.c``, we find::

 static PyMethodDef functions[] = {
   {"axpy", axpy, METH_VARARGS, 0},
   {0, 0, 0, 0}
 };

 DL_EXPORT(void) init_gpaw(void)
 {
   PyObject* m = Py_InitModule3("_gpaw", functions, doc);
   if (m == NULL)
     return;
   import_array();
 }

We could use the C extension code directly as::

  import numpy as np
  import _gpaw
  a = 2.7
  x = np.array([1.1, 1.2, 1.3])
  y = np.zeros(3)
  _gpaw.axpy(a, x, y)
 
Instead, we wrap the code in a Python function ``axpy``
in the file ``gpaw/utilities/blas.py``::

 def axpy(alpha, x, y):
     assert x.shape == y.shape
     assert x.flags.contiguous and y.flags.contiguous
     assert x.dtype == y.dtype
     if isinstance(alpha, complex):
         assert x.dtype == complex
     else:
         assert isinstance(alpha, float)
     _gpaw.axpy(alpha, x, y)

 if not debug:
     axpy = _gpaw.axpy

The Python ``axpy`` function takes care of all value and type checking
of the arguments to the function.  There is therefore no need to do
those checks in the C code (where it would be much more cumbersome to
code).  If the code is run in production mode (``debug == False``,
default), then the Python wrapper function is bypassed for calls to
the function.


.. _Python documentation: http://docs.python.org/ext/ext.html
