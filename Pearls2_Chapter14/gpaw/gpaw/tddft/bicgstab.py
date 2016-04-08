# Written by Lauri Lehtovaara, 2008

"""This module defines BiCGStab-class, which implements biconjugate
gradient stabilized method. Requires Numpy and GPAW's own BLAS."""

import numpy as np

from gpaw.utilities.blas import axpy
from gpaw.utilities.blas import dotc
from gpaw.mpi import rank

class BiCGStab:
    """Biconjugate gradient stabilized method
    
    This class solves a set of linear equations A.x = b using biconjugate 
    gradient stabilized method (BiCGStab). The matrix A is a general, 
    non-singular matrix, e.g., it can be nonsymmetric, complex, and 
    indefinite. The method requires only access to matrix-vector product 
    A.x = b, which is called A.dot(x). Thus A must provide the member 
    function dot(self,x,b), where x and b are complex arrays 
    (numpy.array([], complex), and x is the known vector, and 
    b is the result.

    Now x and b are multivectors, i.e., list of vectors.
    """ 
    
    def __init__( self, gd, timer = None,
                  tolerance = 1e-15, max_iterations = 1000, eps=1e-15 ):
        """Create the BiCGStab-object.
        
        Tolerance should not be smaller than attainable accuracy, which is 
        order of kappa(A) * eps, where kappa(A) is the (spectral) condition 
        number of the matrix. The maximum number of iterations should be 
        significantly less than matrix size, approximately 
        .5 sqrt(kappa) ln(2/tolerance). A small number is treated as zero
        if it's magnitude is smaller than argument eps.
        
        Parameters
        ----------
        gd: GridDescriptor
            grid descriptor for coarse (pseudowavefunction) grid
        timer: Timer
            timer
        tolerance: float
            tolerance for the norm of the residual ||b - A.x||^2
        max_iterations: integer
            maximum number of iterations
        eps: float
            if abs(rho) or (omega) < eps, it's regarded as zero 
            and the method breaks down

        """
        
        self.tol = tolerance
        self.max_iter = max_iterations
        if ( eps <= tolerance ):
            self.eps = eps
        else:
            raise RuntimeError("BiCGStab method got invalid tolerance (tol = %le < eps = %le)." % (tolerance,eps))

        self.iterations = -1
        
        self.gd = gd
        self.timer = timer
        

    def solve(self, A, x, b):
        """Solve a set of linear equations A.x = b.
        
        Parameters:
        A           matrix A
        x           initial guess x_0 (on entry) and the result (on exit)
        b           right-hand side (multi)vector

        """
        if self.timer is not None:
            self.timer.start('BiCGStab')

        # number of vectors
        nvec = len(x)

        # r_0 = b - A x_0
        r = self.gd.zeros(nvec, dtype=complex)
        A.dot(-x,r)
        r += b

        q = self.gd.zeros(nvec, dtype=complex)
        q[:] = r

        p = self.gd.zeros(nvec, dtype=complex)
        v = self.gd.zeros(nvec, dtype=complex)
        t = self.gd.zeros(nvec, dtype=complex)
        m = self.gd.zeros(nvec, dtype=complex)
        alpha = np.zeros((nvec,), dtype=complex) 
        rho  = np.zeros((nvec,), dtype=complex) 
        rhop  = np.zeros((nvec,), dtype=complex) 
        omega = np.zeros((nvec,), dtype=complex) 
        scale = np.zeros((nvec,), dtype=complex) 
        tmp = np.zeros((nvec,), dtype=complex) 

        rhop[:] = 1.
        omega[:] = 1.

        # Multivector dot product, a^H b, where ^H is conjugate transpose
        def multi_zdotc(s, x,y, nvec):
            for i in range(nvec):
                s[i] = dotc(x[i],y[i])
            self.gd.comm.sum(s)
            return s
        # Multivector ZAXPY: a x + y => y
        def multi_zaxpy(a,x,y, nvec):
            for i in range(nvec):
                axpy(a[i]*(1+0J), x[i], y[i])
        # Multiscale: a x => x
        def multi_scale(a,x, nvec):
            for i in range(nvec):
                x[i] *= a[i]

        # scale = square of the norm of b
        multi_zdotc(scale, b,b, nvec)
        scale = np.abs( scale )

        # if scale < eps, then convergence check breaks down
        if (scale < self.eps).any():
            raise RuntimeError("BigCGStab method detected underflow for squared norm of right-hand side (scale = %le < eps = %le)." % (scale,eps))

        #print 'Scale = ', scale

        slow_convergence_iters = 50

        for i in range(self.max_iter):
            # rho_i-1 = q^H r_i-1
            multi_zdotc(rho, q, r, nvec)

            #print 'Rho = ', rho

            # if i=1, p_i = r_i-1
            # else beta = (rho_i-1 / rho_i-2) (alpha_i-1 / omega_i-1)
            #      p_i = r_i-1 + b_i-1 (p_i-1 - omega_i-1 v_i-1)
            beta = (rho / rhop) * (alpha / omega)

            #print 'Beta = ', beta

            # if abs(beta) / scale < eps, then BiCGStab breaks down
            if ( (i > 0) and
                 ((np.abs(beta) / scale) < self.eps).any() ):
                raise RuntimeError("Biconjugate gradient stabilized method failed (abs(beta)=%le < eps = %le)." % (np.min(np.abs(beta)),self.eps))


            # p = r + beta * (p - omega * v)
            # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            multi_zaxpy(-omega, v, p, nvec)
            multi_scale(beta, p, nvec)
            p += r
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            # v_i = A.(M^-1.p)
            A.apply_preconditioner(p,m)
            A.dot(m,v)
            # alpha_i = rho_i-1 / (q^H v_i)
            multi_zdotc(alpha, q,v, nvec)
            alpha = rho / alpha
            # s = r_i-1 - alpha_i v_i
            multi_zaxpy(-alpha, v, r, nvec)
            # s is denoted by r

            #print 'Alpha = ', alpha

            # x_i = x_i-1 + alpha_i (M^-1.p_i) + omega_i (M^-1.s)
            # next line is x_i = x_i-1 + alpha (M^-1.p_i)
            multi_zaxpy(alpha, m, x, nvec)

            # if ( |s|^2 < tol^2 ) done
            multi_zdotc(tmp, r,r, nvec)
            if ( (np.abs(tmp) / scale) < self.tol*self.tol ).all():
                #print 'R2 of proc #', rank, '  = ' , tmp, \
                #    ' after ', i+1, ' iterations'
                break

            # print if slow convergence
            if ((i+1) % slow_convergence_iters) == 0:
                print('Log10 S2 of proc #', rank, '  = ' , np.round(np.log10(np.abs(tmp)),1), \
                      ' after ', i+1, ' iterations')

            # t = A.(M^-1.s), M = 1
            A.apply_preconditioner(r,m)
            A.dot(m,t)
            # omega_i = t^H s / (t^H t)
            multi_zdotc(omega, t,r, nvec)
            multi_zdotc(tmp, t,t, nvec)
            omega = omega / tmp

            #print 'Omega = ', omega

            # x_i = x_i-1 + alpha_i (M^-1.p_i) + omega_i (M^-1.s)
            # next line is x_i = ... + omega_i (M^-1.s)
            multi_zaxpy(omega, m, x, nvec)
            # r_i = s - omega_i * t
            multi_zaxpy(-omega, t, r, nvec)
            # s is no longer denoted by r

            # if ( |r|^2 < tol^2 ) done
            multi_zdotc(tmp, r,r, nvec)
            if ( (np.abs(tmp) / scale) < self.tol*self.tol ).all():
                #print 'R2 of proc #', rank, '  = ' , tmp, \
                #    ' after ', i+1, ' iterations'
                break

            # print if slow convergence
            if ((i+1) % slow_convergence_iters) == 0:
                print('Log10 R2 of proc #', rank, '  = ' , np.round(np.log10(np.abs(tmp)),1), \
                      ' after ', i+1, ' iterations')

            # if abs(omega) < eps, then BiCGStab breaks down
            if ( (np.abs(omega) / scale) < self.eps ).any():
                raise RuntimeError("Biconjugate gradient stabilized method failed (abs(omega)/scale=%le < eps = %le)." % (np.min(np.abs(omega)) / scale, self.eps))
            # finally update rho
            rhop[:] = rho


        # if max iters reached, raise error
        if (i >= self.max_iter-1):
            raise RuntimeError("Biconjugate gradient stabilized method failed to converged within given number of iterations (= %d)." % self.max_iter)


        # done
        self.iterations = i+1
        #print 'BiCGStab iterations = ', self.iterations

        if self.timer is not None:
            self.timer.stop('BiCGStab')

        return self.iterations
        #print self.iterations

