#ifndef QPHIX_INVBICGSTAB_H
#define QPHIX_INVBICGSTAB_H

///// USER SERVICEABLE OPTIONS -- SHOULD BE MOVED TO AUTOCONF CONTROL
#define QPHIX_VERBOSE_BICGSTAB
#define QPHIX_TIMING_BICGSTAB




#include "qphix/linearOp.h"
#include "qphix/blas_new_c.h"
#include "qphix/print_utils.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"


namespace QPhiX
{

  template<typename FT, int V, int S, bool compress12>
  class InvBiCGStab : public AbstractSolver<FT,V,S,compress12> {
  public:
    typedef typename Geometry<FT,V,S,compress12>::FourSpinorBlock Spinor;
    InvBiCGStab(EvenOddLinearOperator<FT,V,S,compress12>& M_,
		int MaxIters_) : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_)
    {
      r=geom.allocCBFourSpinor();
      r0=geom.allocCBFourSpinor();
      p=geom.allocCBFourSpinor();
      v=geom.allocCBFourSpinor();
      t=geom.allocCBFourSpinor();

      norm2Threads= geom.getNSIMT();
      xmyThreads  = geom.getNSIMT();
      copyThreads = geom.getNSIMT();
      zeroThreads = geom.getNSIMT();
      innerProductThreads = geom.getNSIMT();
      pUpdateThreads = geom.getNSIMT();
      sUpdateThreads = geom.getNSIMT();
      rxUpdateThreads = geom.getNSIMT();

    }
    
    ~InvBiCGStab() {
      geom.free(r);
      geom.free(r0);
      geom.free(p);
      geom.free(v);
      geom.free(t);
    }


    void operator()(Spinor *x, 
		    const Spinor *rhs, 
		    const double RsdTarget,
		    int& n_iters, 
		    double& rsd_sq_final, 
		    unsigned long& site_flops,
		    unsigned long& mv_apps, 
		    int isign,
		    bool verbose)
    {
      site_flops = 0;
      mv_apps = 0;
      double rhs_sq;
      double r_norm;
      // Double chi_sq = norm2(chi,s)
      norm2Spinor<FT,V,S,compress12>(rhs_sq,rhs,geom,norm2Threads);
      site_flops +=4*12;

      
      double rsd_sq = rhs_sq * RsdTarget*RsdTarget;
      if( verbose ) masterPrintf("BICGSTAB: ||b||^2 = %e Target Rsd= %e\n", rhs_sq, rsd_sq);

      // Compute r=r0=rhs - A x
      // A(r0,psi,isign)
      M(r0,x,isign);
      norm2Spinor<FT,V,S,compress12>(r_norm,r0,geom,norm2Threads);
      mv_apps++;

      // r = chi - r0
      bicgstab_xmy<FT,V,S,compress12>(rhs,r0,geom, xmyThreads);
      site_flops += 24;

      // r = r0
      copySpinor<FT,V,S,compress12>(r,r0,geom, copyThreads);

      // Now intitialize p and v to zero.
      zeroSpinor<FT,V,S,compress12>(p,geom,zeroThreads);
      zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);

      // rho_prev = 1
      double rho_prev_c[2] = { (double)1, (double)0 };

      // alpha = 1
      double alpha_c[2] = { (double)1, (double)0 };
      
      // omega = 1
      double omega_c[2] = { (double)1, (double)0 };
      
      double rho_c[2] = { (double)0, (double)0 };

      bool notConvP = true;
      int k=1;
      // Now iterate
      for(k=1; (k <= MaxIters) && notConvP; k++) { 
	
	// rho_{k+1} = < r_0 | r > 
	innerProduct<FT,V,S,compress12>(rho_c, r0, r, geom, innerProductThreads);
	site_flops += 8*12;

	if ( (rho_c[0] == 0) && ( rho_c[1] == 0 ) ) {
	  masterPrintf( "BICGSTAB: Breakdown in iteration %d. rho = 0\n", k);

	  rsd_sq_final = r_norm;
	  n_iters = k;
	  return;
	}

	double beta_c[2];
	// beta = ( rho_{k+1}/ rho_{k} ) ( alpha / omega )
	double tmp_c[2];
	double tmp2_c[2];
	complex_div( tmp_c, rho_c, rho_prev_c);
	complex_div( tmp2_c, alpha_c, omega_c );
	complex_mul( beta_c, tmp_c, tmp2_c );

	// p = r + beta(p - omega v)
	// p = r + betap - beta*omega v;
	//   y = x + a(y - bz) 
	
	bicgstab_p_update<FT,V,S,compress12>(r,p,v,beta_c,omega_c,geom, pUpdateThreads);
	site_flops += 16*12;

	// v = Ap
	M(v,p,isign);
	mv_apps++;

	innerProduct<FT,V,S,compress12>(tmp_c,r0,v,geom, innerProductThreads);
	site_flops += 8*12;

	if(  (tmp_c[0] == 0) && (tmp_c[1] == 0) ) { 
	  masterPrintf( "BICGSTAB: Breakdown in iteration %d. <r_0|v> = 0\n", k);
	  n_iters = k;
	  return;
	}
	
	// Alpha = rho/<r0,v>
	complex_div(alpha_c, rho_c, tmp_c);

	// Preserve rho as rho_prev
	rho_prev_c[0] = rho_c[0];
	rho_prev_c[1] = rho_c[1];

	// s = r - alpha v
	// I can overlap s with r because I recompute r at the end
		
	// r = s = r - alpha v:   complex y=ax+y with a=-alpha x=v, y=r
	//	FT alpha_cr[2] = { (FT)(alpha_c[0]), (FT)(alpha_c[1]) };		

	bicgstab_s_update<FT,V,S,compress12>(alpha_c, r, v, geom, sUpdateThreads);
	site_flops += 8*12;

	// t = As
	M(t,r,isign); 
	mv_apps++;

	double t_norm=0;
	norm2Spinor<FT,V,S,compress12>(t_norm, t,geom, norm2Threads);
	site_flops += 4*12;

	if( t_norm == 0 ) { 
	  masterPrintf( "BICGSTAB: Breakdown in iteration %d. ||t|| = 0\n", k);
	  rsd_sq_final = r_norm;
	  n_iters = k;
	  return;
	}
	
	innerProduct<FT,V,S,compress12>(omega_c, t, r, geom, innerProductThreads);
	site_flops += 8*12;

	omega_c[0] /= t_norm;
	omega_c[1] /= t_norm;

	// x = omega r + x +  alpha p 
	// r = r - omega t            
	// r_norm = norm2(r);          

	bicgstab_rxupdate<FT,V,S,compress12>(x, r, t, p, omega_c,alpha_c, r_norm, geom, rxUpdateThreads);
	site_flops += 28*12;

	if( verbose ) masterPrintf("BICGSTAB: iter %d r_norm = %e  target = %e \n" , k,r_norm, rsd_sq);
	if( r_norm < rsd_sq ) { 
	  notConvP = false; // Converged
	}

      } // Loop over iterations
      
      rsd_sq_final=r_norm;
      n_iters = k;
      
      if( notConvP == true ) { 
	masterPrintf("Solver did not converge in %d iterations\n", k);
      }
      return;
      
    }


    void tune(void) 
    {
      int iters=100;
      tuneZeroThreads(iters); 
      tuneNorm2Threads(iters);
      tuneXMYThreads(iters);
      tuneCopyThreads(iters);
      tuneInnerProductThreads(iters);
      tunePUpdateThreads(iters);
      tuneSUpdateThreads(iters);
      tuneRXUpdateThreads(iters);
      reportTuning();
    }

    void reportTuning()
    {
      masterPrintf("TuningResults: \n");
      masterPrintf("\t zeroThreads=%d threads\n", zeroThreads);
      masterPrintf("\t copyThreads=%d threads\n", copyThreads);
      masterPrintf("\t xmyThreads=%d threads\n", xmyThreads);
      masterPrintf("\t norm2Threads=%d threads\n", norm2Threads);
      masterPrintf("\t innerProductThreads=%d threads\n", innerProductThreads);
      masterPrintf("\t pUpdateThreads=%d threads\n", pUpdateThreads);
      masterPrintf("\t sUpdateThreads=%d threads\n", sUpdateThreads);
      masterPrintf("\t rxUpdateThreads=%d threads\n", rxUpdateThreads);
    }

    Geometry<FT,V,S,compress12>& getGeometry(){
      return geom;
    }


    private:

    EvenOddLinearOperator<FT, V,S,compress12>& M;
    Geometry<FT,V,S,compress12>& geom;
    int MaxIters;

    inline
    void complex_div(double res[2], double l[2], double r[2])
    {
      double tmp = (double)1/(r[0]*r[0] + r[1]*r[1]);

      res[0] = (l[0]*r[0] + l[1]*r[1])*tmp;
      res[1] = (l[1]*r[0] - l[0]*r[1])*tmp;

    }

    inline
    void complex_mul(double res[2], double mul1[2], double mul2[2]) 
    {
      res[0]=mul1[0]*mul2[0] - mul1[1]*mul2[1];
      res[1]=mul1[0]*mul2[1] + mul1[1]*mul2[0];
    }


    Spinor *r;
    Spinor *r0;
    Spinor *p;
    Spinor *v;
    Spinor *t;


    int norm2Threads;
    int xmyThreads;
    int copyThreads;
    int zeroThreads;
    int innerProductThreads;
    int pUpdateThreads;
    int sUpdateThreads;
    int rxUpdateThreads;

    void tuneNorm2Threads(int iters) {
      if( r != 0x0 ) { 
	// Do first with 1 thread
	double rnorm;
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	norm2Threads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) { 
	  norm2Spinor(rnorm,r,geom,norm2Threads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneNorm2Threads: threads = %d, current_time=%g (s)\n", norm2Threads,  best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    norm2Spinor(rnorm,r,geom,threads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneNorm2Threads: threads = %d, current_time = %g (s), best=%g (s)\n", threads, current_time, best_time);

	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    norm2Threads=threads;
	  }
	}
      }
    }
 
    void  tuneXMYThreads(int iters)     
    {
      if( r != 0x0 && v != 0 ) { 
	// Do first with 1 thread

	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	xmyThreads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  bicgstab_xmy<FT,V,S,compress12>(r,v,geom, xmyThreads); 
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneXMYThreads: threads = %d, current_time=%g (s)\n", xmyThreads, best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    bicgstab_xmy<FT,V,S,compress12>(r,v,geom, threads); 
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneXMYThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    xmyThreads=threads;
	  }
	}
      }
    }

    void  tuneCopyThreads(int iters)
    {
      if( r != 0x0 && v != 0 ) { 
	// Do first with 1 thread
	
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	copyThreads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  copySpinor<FT,V,S,compress12>(r,v,geom, copyThreads); 
	}
	
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneCopyThreads: threads = %d, current_time=%g (s)\n", copyThreads, best_time);

	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    copySpinor<FT,V,S,compress12>(r,v,geom, threads); 
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneCopyThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    copyThreads=threads;
	  }
	}
      }
    }

    void  tuneZeroThreads(int iters) 
    {
      if( r != 0x0) { 

	zeroThreads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  zeroSpinor<FT,V,S,compress12>(r,geom, zeroThreads); 
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneZeroThreads: threads = %d, current_time=%g (s)\n", zeroThreads, best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    zeroSpinor<FT,V,S,compress12>(r,geom, threads); 
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneZeroThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    zeroThreads=threads;
	  }
	}
      }
    
    }

    void  tuneInnerProductThreads(int iters) {
     if( r != 0x0 && v != 0 ) { 
	// Do first with 1 thread
	
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	innerProductThreads=1;
	double iprod_c[2];
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  innerProduct<FT,V,S,compress12>(iprod_c, r,v,geom, innerProductThreads); 
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneInnerProductThreads: threads = %d, current_time=%g (s)\n", innerProductThreads, best_time);

	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    innerProduct<FT,V,S,compress12>(iprod_c, r,v,geom, threads); 
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneInnerProductThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    innerProductThreads=threads;
	  }
	}
      }
    }

    void  tunePUpdateThreads(int iters) 
    {
      if( r != 0x0 && v != 0 && p !=0 ) { 
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(p,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	pUpdateThreads=1;
	double  beta_cr[2]={(double)1.0, (double)0.5};
	double omega_cr[2]={(double)2.0, (double)-0.5};
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  bicgstab_p_update<FT,V,S,compress12>(r,p,v,beta_cr,omega_cr,geom, pUpdateThreads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;

	masterPrintf("tunePUpdateThreads: threads = %d, current_time=%g (s)\n", pUpdateThreads, best_time);
	
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    bicgstab_p_update<FT,V,S,compress12>(r,p,v,beta_cr,omega_cr,geom, threads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  

	  masterPrintf("tunePUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    pUpdateThreads=threads;
	  }
	}
      }
    }
    
    void  tuneSUpdateThreads(int iters) 
    {
     if( r != 0x0 && v != 0 ) { 
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	sUpdateThreads=1;
	double alpha_cr[2]={(double)1.0, (double)0.5};
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  bicgstab_s_update<FT,V,S,compress12>(alpha_cr, r, v, geom, sUpdateThreads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;

	masterPrintf("tuneSUpdateThreads: threads = %d, current_time=%g (s)\n", sUpdateThreads, best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    bicgstab_s_update<FT,V,S,compress12>(alpha_cr, r, v, geom, threads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneSUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    sUpdateThreads=threads;
	  }
	}
      }
    }


    void tuneRXUpdateThreads(int iters) 
    {
      if( r != 0x0 && v != 0 && p!=0 && t !=0 ) { 
	zeroSpinor<FT,V,S,compress12>(r,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(v,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(t,geom,zeroThreads);
	zeroSpinor<FT,V,S,compress12>(p,geom,zeroThreads);
	double r_norm=0;
	rxUpdateThreads=1;
	double alpha_cr[2]={(double)1.0, (double)0.5};
	double omega_cr[2]={(double)-1.0, (double)-0.25};
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  bicgstab_rxupdate<FT,V,S,compress12>(r, v, t, p, omega_cr,alpha_cr, r_norm, geom, rxUpdateThreads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;

	masterPrintf("tuneRXUpdateThreads: threads = %d, current_time=%g (s)\n", rxUpdateThreads, best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    bicgstab_rxupdate<FT,V,S,compress12>(r, v, t, p, omega_cr,alpha_cr, r_norm, geom, threads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneRXUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    rxUpdateThreads=threads;
	  }
	}
      }
    }

  };





};



#endif
