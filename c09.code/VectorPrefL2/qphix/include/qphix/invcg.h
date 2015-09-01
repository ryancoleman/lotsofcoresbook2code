#ifndef QPHIX_INVCG2_H
#define QPHIX_INVCG2_H

///// USER SERVICEABLE OPTIONS -- SHOULD BE MOVED TO AUTOCONF CONTROL
#define QPHIX_TIMING_CG




#include "qphix/linearOp.h"
#include "qphix/print_utils.h"
#include "qphix/blas_new_c.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

namespace  QPhiX
{

  // FIXME!!!: passing compress12 as a template to solver breaks 'delegation of responsibilities rule.
  // Solution may be to take it out of geometry and create a fields class. Geometry can then deal
  // exclusively with the blocking and stuff...
  //
  // That will be a second refactor step when everything works
  template<typename FT, int veclen, int soalen, bool compress12>
    class InvCG : public AbstractSolver<FT,veclen,soalen,compress12> {
  public: 
    typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock Spinor;
    InvCG(EvenOddLinearOperator<FT,veclen,soalen,compress12>& M_,
	  int MaxIters_) : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_)
    {
      
      // Length of vectors in floats, including PADs
      masterPrintf("Initializing CG Solver: Nvec=%d Ny=%d Nz=%d Nt=%d\n", geom.nVecs(), geom.Ny(), geom.Nz(), geom.Nt());
      
      
      mp = (Spinor *)geom.allocCBFourSpinor();
      mmp = (Spinor *)geom.allocCBFourSpinor();
      p = (Spinor *)geom.allocCBFourSpinor();
      r = (Spinor *)geom.allocCBFourSpinor();
      

      copyThreads=geom.getNSIMT();
      aypxThreads=geom.getNSIMT();
      xmyNormThreads=geom.getNSIMT();
      rmammpNorm2rxpapThreads=geom.getNSIMT(); 
      norm2Threads=geom.getNSIMT(); 
      
      
    }

    ~InvCG() {
      geom.free(mp);
      geom.free(mmp);
      geom.free(p);
      geom.free(r);
    }
    

    void operator()(Spinor *x, 
		    const Spinor *rhs, 
		    const double RsdTarget,
		    int& n_iters, 
		    double& rsd_sq_final, 
		    unsigned long& site_flops,
		    unsigned long& mv_apps, 
		    int isign,
		    bool verboseP) 
    {
#ifdef QPHIX_TIMING_CG
      TSC_tick cg_start;
      TSC_tick cg_end;
      TSC_tick cg_time=0;

      // LINOP_TIMINGS
      TSC_tick lstart;
      TSC_tick lend;
      TSC_tick ltime=0;
      
      // BLAS TIMINGS: norm2Spinor
      TSC_tick norm2_start;
      TSC_tick norm2_end;
      TSC_tick norm2_time=0;
      
      TSC_tick blas_time = 0;
      // BLAS TIMINGS: xmyNorm 
      TSC_tick xmyNorm_start;
      TSC_tick xmyNorm_end;
      TSC_tick xmyNorm_time=0;

     // BLAS_TIMINGS: aypx
      TSC_tick aypx_start;
      TSC_tick aypx_end;
      TSC_tick aypx_time=0;

      // BLAS_TIMINGS: copy
      TSC_tick copy_start;
      TSC_tick copy_end;
      TSC_tick copy_time = 0;

      // BLAS_TIMINGS: rMammpNorm2Rpmax
      TSC_tick rmammp_start;
      TSC_tick rmammp_end;
      TSC_tick rmammp_time =0;



      TSC_tick missing_time;
#endif

      int n_cores;
      int n_simt;
      double cp;
      
      
#ifdef TIMING_CG
      CLOCK_NOW(cg_start);
#endif
      
      site_flops=0;
      mv_apps=0;
      n_cores = geom.getNumCores();
      n_simt = geom.getNSIMT();
      
#ifdef TIMING_CG
      CLOCK_NOW(norm2_start);
#endif
      double chi_sq;


      norm2Spinor<FT,veclen,soalen,compress12>(chi_sq, rhs, geom, norm2Threads);        
      masterPrintf("Chi_sq=%g RsdTarget=%g\n", chi_sq, RsdTarget);

#ifdef CGDEBUG
      norm2Spinor<FT,veclen,soalen,compress12>(chi_sq, x, geom, norm2Threads);
      masterPrintf("||x||=%lf\n", chi_sq);

      norm2Spinor<FT,veclen,soalen,compress12>(chi_sq, rhs, geom, norm2Threads);        
      masterPrintf("Chi_sq=%g RsdTarget=%g\n", chi_sq, RsdTarget);
#endif      

#ifdef TIMING_CG
      CLOCK_NOW(norm2_end);
      norm2_time += (norm2_end-norm2_start);
#endif
      
      site_flops += 4*12;
      double rsd_sq = (RsdTarget*RsdTarget)*chi_sq;
      
      double tmp_d;
      // M^\dagger M psi
#ifdef TIMING_CG
      CLOCK_NOW(lstart);
#endif
      M(mp,x,isign);

#ifdef CGDEBUG
      norm2Spinor<FT,veclen,soalen,compress12>(tmp_d,mp,geom,norm2Threads);
      masterPrintf("M p = %lf \n", tmp_d);
#endif

      M(mmp,mp,-isign);

#ifdef CGDEBUG
      norm2Spinor<FT,veclen,soalen,compress12>(tmp_d,mmp,geom,norm2Threads);
      masterPrintf("MM p = %lf \n", tmp_d);
#endif

#ifdef TIMING_CG 
      CLOCK_NOW(lend);
      ltime += (lend - lstart);
#endif
      mv_apps += 2;
      
      
      // r = chi_internal - mmp,   cp=norm2(r)
#ifdef TIMING_CG
      CLOCK_NOW(xmyNorm_start);
#endif
      xmyNorm2Spinor<FT,veclen,soalen,compress12>(r,rhs,mmp, cp, geom, xmyNormThreads);
      
#ifdef TIMING_CG
      CLOCK_NOW(xmyNorm_end);
      xmyNorm_time += (xmyNorm_end - xmyNorm_start);
#endif
      site_flops += 6*12;
      
      
      if( verboseP ) masterPrintf("CG: r0 = %g   target = %g \n", cp, rsd_sq);

      if( cp <= rsd_sq ) {
	n_iters=0;
	rsd_sq_final = cp / chi_sq;

#ifdef TIMING_CG
	CLOCK_NOW(cg_end);
	cg_time = cg_end - cg_start;

	masterPrintf("CG Total time: %ld ticks", cg_time );
	masterPrintf("Linop Total time: %ld ticks %u applications, %g percent of CG\n", ltime, mv_apps, 100.0*((double)ltime)/((double)cg_time));
	blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
	masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",  blas_time, 100.0*((double)blas_time)/((double)cg_time));
	masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n", norm2_time, 100.0*((double)norm2_time)/((double)cg_time));
	masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n", xmyNorm_time, 100.0*((double)xmyNorm_time)/((double)cg_time));
	masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n", aypx_time, 100.0*((double)aypx_time)/((double)cg_time));
	masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n", rmammp_time, 100.0*((double)rmammp_time)/((double)cg_time));
	masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n", copy_time, 100.0*((double)copy_time)/((double)cg_time));
	missing_time = cg_time - ltime -blas_time;
	masterPrintf("Missing time = %ld ticks = %g percent of CG\n", missing_time, 100.0* (double)(missing_time)/ ((double)cg_time));
#endif
	return;

      }

#ifdef TIMING_CG
      CLOCK_NOW(copy_start);
#endif

      copySpinor<FT,veclen,soalen,compress12>(p,r, geom, copyThreads);

#ifdef TIMING_CG
      CLOCK_NOW(copy_end);
      copy_time += (copy_end - copy_start);
#endif

      double a,b,c,d;
      for(int k=0; k <= MaxIters; ++k) {
	c = cp;

#ifdef TIMING_CG 
	CLOCK_NOW(lstart);
#endif
	M(mp,p, isign);
	M(mmp,mp,-isign);

#ifdef TIMING_CG
	CLOCK_NOW(lend);
	ltime += (lend - lstart);
#endif
	mv_apps +=2;
	

#ifdef TIMING_CG
	CLOCK_NOW(norm2_start);
#endif

	norm2Spinor<FT,veclen,soalen,compress12>(d,mp,geom,norm2Threads);

#ifdef TIMING_CG
	CLOCK_NOW(norm2_end);
	norm2_time += (norm2_end - norm2_start);
#endif

	site_flops+=4*12;

	a = c/d;

#ifdef CGDEBUG
	masterPrintf("CG:   iter %d c=%lf d=%lf a=%lf \n", k,c,d,a);
#endif

#ifdef TIMING_CG
	CLOCK_NOW(rmammp_start);
#endif

	rmammpNorm2rxpap<FT,veclen,soalen,compress12>(r,a,mmp,cp,x,p,geom, rmammpNorm2rxpapThreads);

#ifdef CGDEBUG
	norm2Spinor<FT,veclen,soalen,compress12>(tmp_d, r, geom, norm2Threads);        
	masterPrintf("CG:   iter %d: r2 = %lf \n", k, tmp_d);

	norm2Spinor<FT,veclen,soalen,compress12>(tmp_d, x, geom, norm2Threads);        
	masterPrintf("CG:   iter %d: x2 = %lf \n", k, tmp_d);
#endif

#ifdef TIMING_CG
	CLOCK_NOW(rmammp_end);
	rmammp_time += ( rmammp_end - rmammp_start);
#endif
	site_flops += 12*12;

	if (verboseP) masterPrintf("CG: iter %d:  r2 = %g   target = %g \n", k, cp, rsd_sq);	  

	if( cp < rsd_sq ) { 
	  n_iters = k;

#ifdef TIMING_CG
	  CLOCK_NOW(lstart);
	  M(mp,x,isign);
	  M(mmp,mp,-isign);
	  CLOCK_NOW(lend);
	  ltime += (lend -lstart);
	  mv_apps += 2;

	  CLOCK_NOW(xmyNorm_start);
	  xmyNorm2Spinor<FT,veclen,soalen,compress12>(r, rhs, mmp, cp, geom, xmyNormThreads);
	  CLOCK_NOW(xmyNorm_end);
	  xmyNorm_time += (xmyNorm_end - xmyNorm_start);

	  site_flops += 6*12;
	  rsd_sq_final = cp ;
	  CLOCK_NOW(cg_end);
	  cg_time = cg_end - cg_start;

	  masterPrintf("Iters=%d Final Rsd = %e Final RsdSq=%e Final ResidRelSq=%e Final ResidRel=%e\n", n_iters, sqrt(rsd_sq_final), rsd_sq_final, rsd_sq_final/chi_sq, sqrt(rsd_sq_final/chi_sq));
	  masterPrintf("CG Total time: %ld ticks\n", cg_time);
	  masterPrintf("Linop Total time: %ld ticks %u applications, %g percent of CG\n", ltime, mv_apps, 100.0*((double)ltime)/((double)cg_time));
	  blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
	  masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",  blas_time, 100.0*((double)blas_time)/((double)cg_time));
	  masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n", norm2_time, 100.0*((double)norm2_time)/((double)cg_time));
	  masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n", xmyNorm_time, 100.0*((double)xmyNorm_time)/((double)cg_time));
	  masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n", aypx_time, 100.0*((double)aypx_time)/((double)cg_time));
	  masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n", rmammp_time, 100.0*((double)rmammp_time)/((double)cg_time));
	  masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n", copy_time, 100.0*((double)copy_time)/((double)cg_time));
	  missing_time = cg_time - ltime -blas_time;
	  masterPrintf("Missing time = %ld ticks = %g percent of CG \n", missing_time, 100.0*(double)(missing_time)/ ((double)cg_time));
#endif

	  return;
	}

	b = cp/c;

#ifdef TIMING_CG
	CLOCK_NOW(aypx_start);
#endif
	aypx<FT,veclen,soalen,compress12>(b,r,p, geom, aypxThreads);

#ifdef TIMING_CG
	CLOCK_NOW(aypx_end);
	aypx_time += (aypx_end - aypx_start);
#endif
	site_flops += 4*12;
      }
      

      n_iters = MaxIters;
      rsd_sq_final = cp;
 
#ifdef TIMING_CG
      CLOCK_NOW(cg_end);
      cg_time=cg_end - cg_start;

      masterPrintf("CG Total time: %ld ticks\n", cg_time);
      masterPrintf("Linop Total time: %ld ticks %u applications, %g percent of CG\n", ltime, mv_apps, 100.0*((double)ltime)/((double)cg_time));
      blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
      masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",  blas_time, 100.0*((double)blas_time)/((double)cg_time));
      masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n", norm2_time, 100.0*((double)norm2_time)/((double)cg_time));
      masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n", xmyNorm_time, 100.0*((double)xmyNorm_time)/((double)cg_time));
      masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n", aypx_time, 100.0*((double)aypx_time)/((double)cg_time));
      masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n", rmammp_time, 100.0*((double)rmammp_time)/((double)cg_time));
      masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n", copy_time, 100.0*((double)copy_time)/((double)cg_time));
      missing_time = cg_time - ltime -blas_time;
      masterPrintf("Missing time = %ld ticks = %g percent of CG\n", missing_time, 100.0* (double)(missing_time)/ ((double)cg_time));
#endif
      return;
    }

    
    void setCopyThreads(int c) { copyThreads = c; }
    void setAypxThreads(int c) { aypxThreads = c; }
    void setXmyNormThreads(int c) { xmyNormThreads = c; }
    void setRmammpNorm2rxpapThreads(int c )  {rmammpNorm2rxpapThreads = c;}
    void setNorm2Threads(int c) { norm2Threads = c; }
    
    int getCopyThreads(void) { return copyThreads; }
    int getAypxThreads(void) { return aypxThreads; }
    int getXmyNormThreads(void) { return xmyNormThreads; }
    int getRmammpNorm2rxpapThreads(void) { return rmammpNorm2rxpapThreads; }
    int getNorm2Threads(void) { return norm2Threads; }


    void tune()
    {
      int iters=100;
#if 1
      tuneCopyThreads(iters);
      tuneAypxThreads(iters);
      tuneNorm2Threads(iters);
      tuneXMYNorm2Threads(iters);
      tuneRXUpdateThreads(iters);

      reportTuning();
#endif
    }

    Geometry<FT,veclen,soalen,compress12>& getGeometry() {
      return geom;
    }

    void reportTuning() {
	masterPrintf("TuningResults: \n");
	masterPrintf("\t copyThreads=%d threads\n", copyThreads);
	masterPrintf("\t aypxThreads=%d threads\n", aypxThreads);
	masterPrintf("\t xmyNormThreads=%d threads\n", xmyNormThreads);
	masterPrintf("\t rmammpNorm2rxpapThreads=%d threads\n", rmammpNorm2rxpapThreads);
	masterPrintf("\t norm2Threads=%d threads\n", norm2Threads);
    }
    private:

    EvenOddLinearOperator<FT, veclen,soalen,compress12>& M;
    Geometry<FT, veclen,soalen,compress12>& geom;
    int MaxIters;

    Spinor *mp;
    Spinor *mmp;
    Spinor *p;
    Spinor *r;


    int copyThreads;
    int aypxThreads;
    int xmyNormThreads;
    int rmammpNorm2rxpapThreads;
    int norm2Threads;



    void  tuneCopyThreads(int iters)
    {
      if( r != 0x0 && p != 0 ) { 
	// Do first with 1 thread
	
	zeroSpinor<FT,veclen,soalen,compress12>(r,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(p,geom,1);
	copyThreads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  copySpinor<FT,veclen,soalen,compress12>(r,p,geom, copyThreads); 
	}
	
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneCopyThreads: threads = %d, current_time=%g (s)\n", copyThreads, best_time);

	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    copySpinor<FT,veclen,soalen,compress12>(r,p,geom, threads); 
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

    void  tuneAypxThreads(int iters)
    {
      if( r != 0x0 && p != 0 ) { 
	// Do first with 1 thread
	
	zeroSpinor<FT,veclen,soalen,compress12>(r,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(p,geom,1);
	aypxThreads=1;
	double br=(double)0.5;

	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  aypx<FT,veclen,soalen,compress12>(br,r,p, geom, aypxThreads);
	}
	
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneAypxThreads: threads = %d, current_time=%g (s)\n", aypxThreads, best_time);

	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    aypx<FT,veclen,soalen,compress12>(br,r,p, geom, threads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneAypxThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    aypxThreads=threads;
	  }
	}
      }
    }

    void tuneNorm2Threads(int iters) {
      if( r != 0x0 ) { 
	// Do first with 1 thread
	double rnorm;
	zeroSpinor<FT,veclen,soalen,compress12>(r,geom,1);
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

    void  tuneXMYNorm2Threads(int iters)     
    {
      if( r != 0x0 && p != 0 && mmp != 0) { 
	// Do first with 1 thread

	zeroSpinor<FT,veclen,soalen,compress12>(r,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(p,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(mmp,geom,1);
	xmyNormThreads=1;
	double reduction=0;

	double start_time=omp_get_wtime();	
	for(int i=0; i < iters; i++) {
	  xmyNorm2Spinor<FT,veclen,soalen,compress12>(r,p,mmp, reduction, geom, xmyNormThreads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;
	masterPrintf("tuneXmyNorm2Threads: threads = %d, current_time=%g (s)\n", xmyNormThreads, best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) {
	    xmyNorm2Spinor<FT,veclen,soalen,compress12>(r,p,mmp, reduction, geom, threads); 
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneXMYThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    xmyNormThreads=threads;
	  }
	}
      }
    }

    void tuneRXUpdateThreads(int iters) 
    {
      if( r != 0x0 && p != 0 && mp!=0 && mmp !=0 ) { 
	zeroSpinor<FT,veclen,soalen,compress12>(r,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(p,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(mp,geom,1);
	zeroSpinor<FT,veclen,soalen,compress12>(mmp,geom,1);
	double r_norm=0;
	double a = (double)0.5;
	int rmammpNorm2rxpapThreads=1;
	double start_time=omp_get_wtime();
	for(int i=0; i < iters; i++) {
	  rmammpNorm2rxpap<FT,veclen,soalen,compress12>(r,a,mmp,r_norm,mp,p,geom, rmammpNorm2rxpapThreads);
	}
	double stop_time=omp_get_wtime();
	double best_time=stop_time - start_time;

	masterPrintf("tuneRXUpdateThreads: threads = %d, current_time=%g (s)\n",rmammpNorm2rxpapThreads , best_time);
	for(int threads = 2; threads <=geom.getNSIMT(); threads++) { 
	  start_time=omp_get_wtime();
	  for(int i=0; i < iters; i++) { 
	    rmammpNorm2rxpap<FT,veclen,soalen,compress12>(r,a,mmp,r_norm,mp,p,geom, rmammpNorm2rxpapThreads);
	  }
	  stop_time=omp_get_wtime();
	  double current_time=stop_time-start_time;
	  
	  masterPrintf("tuneRXUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n", threads, current_time, best_time);
	  if ( current_time < best_time ) { 
	    best_time = current_time;
	    rmammpNorm2rxpapThreads=threads;
	  }
	}
      }
    }


  };





};



#endif
