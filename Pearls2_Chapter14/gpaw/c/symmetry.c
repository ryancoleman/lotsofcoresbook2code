/*  Copyright (C) 2010-2011 CAMd
 *  Please see the accompanying LICENSE file for further information. */
#include "extensions.h"

//
// Apply symmetry operation op_cc to a and add result to b:
// 
//     =T_       _
//   b(U g) += a(g),
// 
// where:
// 
//   =                         _T
//   U     = op_cc[c1, c2] and g = (g0, g1, g2).
//    c1,c2
//
PyObject* symmetrize(PyObject *self, PyObject *args)
{
    PyArrayObject* a_g_obj;
    PyArrayObject* b_g_obj;
    PyArrayObject* op_cc_obj;
    if (!PyArg_ParseTuple(args, "OOO", &a_g_obj, &b_g_obj, &op_cc_obj)) 
        return NULL;

    const long* C = (const long*)PyArray_DATA(op_cc_obj);
    int ng0 = PyArray_DIMS(a_g_obj)[0];
    int ng1 = PyArray_DIMS(a_g_obj)[1];
    int ng2 = PyArray_DIMS(a_g_obj)[2];

    const double* a_g = (const double*)PyArray_DATA(a_g_obj);
    double* b_g = (double*)PyArray_DATA(b_g_obj);
    for (int g0 = 0; g0 < ng0; g0++)
        for (int g1 = 0; g1 < ng1; g1++)
	    for (int g2 = 0; g2 < ng2; g2++) {
	      int p0 = ((C[0] * g0 + C[3] * g1 + C[6] * g2) % ng0 + ng0) % ng0;
	      int p1 = ((C[1] * g0 + C[4] * g1 + C[7] * g2) % ng1 + ng1) % ng1;
	      int p2 = ((C[2] * g0 + C[5] * g1 + C[8] * g2) % ng2 + ng2) % ng2;
              b_g[(p0 * ng1 + p1) * ng2 + p2] += *a_g++;
	    }
    
    Py_RETURN_NONE;
}

PyObject* symmetrize_ft(PyObject *self, PyObject *args)
{
    PyArrayObject* a_g_obj;
    PyArrayObject* b_g_obj;
    PyArrayObject* op_cc_obj;
    PyArrayObject* ft_c_obj;
   
    if (!PyArg_ParseTuple(args, "OOOO", &a_g_obj, &b_g_obj, &op_cc_obj, &ft_c_obj)) 
        return NULL;

    const double* ft = (const double*)PyArray_DATA(ft_c_obj);
    const long* C = (const long*)PyArray_DATA(op_cc_obj);
    int ng0 = PyArray_DIMS(a_g_obj)[0];
    int ng1 = PyArray_DIMS(a_g_obj)[1];
    int ng2 = PyArray_DIMS(a_g_obj)[2];

    int ft0 = (int)(ft[0]*ng0);
    int ft1 = (int)(ft[1]*ng1);
    int ft2 = (int)(ft[2]*ng2);   
    
    const double* a_g = (const double*)PyArray_DATA(a_g_obj);
    double* b_g = (double*)PyArray_DATA(b_g_obj);
    for (int g0 = 0; g0 < ng0; g0++)
        
        for (int g1 = 0; g1 < ng1; g1++)
            for (int g2 = 0; g2 < ng2; g2++) {
              int p0 = ((C[0] * g0 + C[3] * g1 + C[6] * g2 - ft0) % ng0 + ng0) % ng0;
              int p1 = ((C[1] * g0 + C[4] * g1 + C[7] * g2 - ft1) % ng1 + ng1) % ng1;
              int p2 = ((C[2] * g0 + C[5] * g1 + C[8] * g2 - ft2) % ng2 + ng2) % ng2;
              b_g[(p0 * ng1 + p1) * ng2 + p2] += *a_g++;
            }
    
    Py_RETURN_NONE;
}

PyObject* symmetrize_wavefunction(PyObject *self, PyObject *args)
{
    PyArrayObject* a_g_obj;
    PyArrayObject* b_g_obj;
    PyArrayObject* op_cc_obj;
    PyArrayObject* kpt0_obj;
    PyArrayObject* kpt1_obj;
    
    if (!PyArg_ParseTuple(args, "OOOOO", &a_g_obj, &b_g_obj, &op_cc_obj, &kpt0_obj, &kpt1_obj)) 
        return NULL;

    const long* C = (const long*)PyArray_DATA(op_cc_obj);
    const double* kpt0 = (const double*) PyArray_DATA(kpt0_obj);
    const double* kpt1 = (const double*) PyArray_DATA(kpt1_obj);
    int ng0 = PyArray_DIMS(a_g_obj)[0];
    int ng1 = PyArray_DIMS(a_g_obj)[1];
    int ng2 = PyArray_DIMS(a_g_obj)[2];
    
    const double complex* a_g = (const double complex*)PyArray_DATA(a_g_obj);
    double complex* b_g = (double complex*)PyArray_DATA(b_g_obj);
    
    for (int g0 = 0; g0 < ng0; g0++)
        for (int g1 = 0; g1 < ng1; g1++)
	    for (int g2 = 0; g2 < ng2; g2++) {
	      int p0 = ((C[0] * g0 + C[3] * g1 + C[6] * g2) % ng0 + ng0) % ng0;
	      int p1 = ((C[1] * g0 + C[4] * g1 + C[7] * g2) % ng1 + ng1) % ng1;
	      int p2 = ((C[2] * g0 + C[5] * g1 + C[8] * g2) % ng2 + ng2) % ng2;

	      double complex phase = cexp(I * 2. * M_PI * 
					  (kpt1[0]/ng0*p0 +
					   kpt1[1]/ng1*p1 +
					   kpt1[2]/ng2*p2 -
					   kpt0[0]/ng0*g0 -
					   kpt0[1]/ng1*g1 -
					   kpt0[2]/ng2*g2));
	      b_g[(p0 * ng1 + p1) * ng2 + p2] += (*a_g * phase);
	      a_g++;
	    }
    
    Py_RETURN_NONE;
}

PyObject* symmetrize_return_index(PyObject *self, PyObject *args)
{
    PyArrayObject* a_g_obj;
    PyArrayObject* b_g_obj;
    PyArrayObject* op_cc_obj;
    PyArrayObject* kpt0_obj;
    PyArrayObject* kpt1_obj;
    
    if (!PyArg_ParseTuple(args, "OOOOO", &a_g_obj, &b_g_obj, &op_cc_obj, &kpt0_obj, &kpt1_obj)) 
        return NULL;

    const long* C = (const long*)PyArray_DATA(op_cc_obj);
    const double* kpt0 = (const double*) PyArray_DATA(kpt0_obj);
    const double* kpt1 = (const double*) PyArray_DATA(kpt1_obj);

    int ng0 = PyArray_DIMS(a_g_obj)[0];
    int ng1 = PyArray_DIMS(a_g_obj)[1];
    int ng2 = PyArray_DIMS(a_g_obj)[2];
    
    unsigned long* a_g = (unsigned long*)PyArray_DATA(a_g_obj);
    double complex* b_g = (double complex*)PyArray_DATA(b_g_obj);
    
    for (int g0 = 0; g0 < ng0; g0++)
        for (int g1 = 0; g1 < ng1; g1++)
	    for (int g2 = 0; g2 < ng2; g2++) {
	      int p0 = ((C[0] * g0 + C[3] * g1 + C[6] * g2) % ng0 + ng0) % ng0;
	      int p1 = ((C[1] * g0 + C[4] * g1 + C[7] * g2) % ng1 + ng1) % ng1;
	      int p2 = ((C[2] * g0 + C[5] * g1 + C[8] * g2) % ng2 + ng2) % ng2;

	      double complex phase = cexp(I * 2. * M_PI * 
					  (kpt1[0]/ng0*p0 +
					   kpt1[1]/ng1*p1 +
					   kpt1[2]/ng2*p2 -
					   kpt0[0]/ng0*g0 -
					   kpt0[1]/ng1*g1 -
					   kpt0[2]/ng2*g2));
	      *a_g++ = (p0 * ng1 + p1) * ng2 + p2;
	      *b_g++ = phase;
	    }
    
    Py_RETURN_NONE;
}

PyObject* symmetrize_with_index(PyObject *self, PyObject *args)
{
    PyArrayObject* a_g_obj;
    PyArrayObject* b_g_obj;
    PyArrayObject* index_g_obj;
    PyArrayObject* phase_g_obj;
    
    if (!PyArg_ParseTuple(args, "OOOO", &a_g_obj, &b_g_obj, &index_g_obj, &phase_g_obj)) 
        return NULL;

    int ng0 = PyArray_DIMS(a_g_obj)[0];
    int ng1 = PyArray_DIMS(a_g_obj)[1];
    int ng2 = PyArray_DIMS(a_g_obj)[2];
    
    const unsigned long* index_g = (const unsigned long*)PyArray_DATA(index_g_obj);
    const double complex* phase_g = (const double complex*)PyArray_DATA(phase_g_obj);
    const double complex* a_g = (const double complex*)PyArray_DATA(a_g_obj);
    double complex* b_g = (double complex*)PyArray_DATA(b_g_obj);
 
   
    for (int g0 = 0; g0 < ng0; g0++)
        for (int g1 = 0; g1 < ng1; g1++)
	    for (int g2 = 0; g2 < ng2; g2++) {
	      b_g[*index_g] += (*a_g * *phase_g);
	      a_g++;
	      phase_g++;
	      index_g++;
	    }
    
    Py_RETURN_NONE;
}

PyObject* map_k_points(PyObject *self, PyObject *args)
{
    PyArrayObject* bzk_kc_obj;
    PyArrayObject* U_scc_obj;
    double tol;
    PyArrayObject* bz2bz_ks_obj;
    int ka, kb;

    if (!PyArg_ParseTuple(args, "OOdOii", &bzk_kc_obj, &U_scc_obj,
			   &tol, &bz2bz_ks_obj, &ka, &kb)) 
        return NULL;

    const long* U_scc = (const long*)PyArray_DATA(U_scc_obj);
    const double* bzk_kc = (const double*)PyArray_DATA(bzk_kc_obj);
    long* bz2bz_ks = (long*)PyArray_DATA(bz2bz_ks_obj);

    int nbzkpts = PyArray_DIMS(bzk_kc_obj)[0];
    int nsym = PyArray_DIMS(U_scc_obj)[0];

    for (int k1 = ka; k1 < kb; k1++) {
        const double* q = bzk_kc + k1 * 3;
	 for (int s = 0; s < nsym; s++) {
	     const long* U = U_scc + s * 9;
	     double q0 = U[0] * q[0] + U[1] * q[1] + U[2] * q[2];
	     double q1 = U[3] * q[0] + U[4] * q[1] + U[5] * q[2];
	     double q2 = U[6] * q[0] + U[7] * q[1] + U[8] * q[2];
	     for (int k2 = 0; k2 < nbzkpts; k2++) {
	         double p0 = q0 - bzk_kc[k2 * 3];
		 if (fabs(p0 - round(p0)) > tol)
		     continue;
	         double p1 = q1 - bzk_kc[k2 * 3 + 1];
		 if (fabs(p1 - round(p1)) > tol)
		     continue;
	         double p2 = q2 - bzk_kc[k2 * 3 + 2];
		 if (fabs(p2 - round(p2)) > tol)
		     continue;
		 bz2bz_ks[k1 * nsym + s] = k2;
	     }
	 }
    }
    Py_RETURN_NONE;
}
