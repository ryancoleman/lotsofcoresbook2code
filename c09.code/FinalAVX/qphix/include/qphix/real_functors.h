#ifndef QPHIX_REAL_FUNCTORS_H
#define QPHIX_REAL_FUNCTORS_H

#include "qphix/qphix_config.h"
#include "qphix/blas_utils.h"
#include "qphix/print_utils.h"

namespace QPhiX
{

  template<typename FT, int V, int S, bool compress>
  class CopyFunctor {
  public:
    
  CopyFunctor( typename Geometry<FT,V,S,compress>::FourSpinorBlock* res_,
	       const typename Geometry<FT,V,S,compress>::FourSpinorBlock* src_ ) : res(res_), src(src_) {}
    
    ~CopyFunctor() {}
    
    inline void 
    func(int block)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* srcbase=&src[block][0][0][0][0];
      FT* resbase = &res[block][0][0][0][0];

#if defined(__MIC__)
      const int prefdist1 = 12;
      const char* pref1base = (const char *)srcbase+prefdist1*64;
      
      const int prefdist2 = 64;
      const char* pref2base = (const char *)srcbase+prefdist2*64;
#endif
      
      for(int numv=0; numv < nvec_in_spinor; numv++) { 
#if defined(__MIC__)
	_mm_prefetch(&pref1base[numv*V*sizeof(FT)], _MM_HINT_T0);
	_mm_prefetch(&pref2base[numv*V*sizeof(FT)], _MM_HINT_T1);
#endif
	
#pragma simd
#pragma vector aligned(srcbase)
#pragma vector nontemporal(resbase)
	for(int s=0; s < V; s++) {
	  resbase[numv*V + s ] = srcbase[ numv*V + s ];
	}
      }
    }
    
  private: 
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict src;
  };
  
  
  template<typename FT, int V, int S, bool compress>
  class ZeroFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;
    
    ZeroFunctor( typename Geometry<FT,V,S,compress>::FourSpinorBlock* res_) : res(res_) {}
    ~ZeroFunctor() {}
    
    inline void 
    func(int block)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      FT* resbase = &res[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock res_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock res_spinor;
#endif

       // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN
            for(int i=0; i < S; i++) {
              res_spinor[col][spin][reim][i] = rep<AT,double>((double)0);
            }
#else
	   res_spinor[col][spin][reim][:] = rep<AT,double>((double)0);
#endif
	  }
	}
      }
      
      BLASUtils::streamOutSpinor<FT,V>(resbase, (const AT *)res_spinor, nvec_in_spinor);

    }
    
  private: 
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res;
  };
  
  template<typename FT, int V, int S, bool compress>
  class AXPYFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;

    AXPYFunctor(double a_, 
		const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_,
		typename Geometry<FT,V,S,compress>::FourSpinorBlock* y_) : a(rep<AT,double>(a_)), x(x_), y(y_) {}
    ~AXPYFunctor() {}
    
    inline void 
    func(int block)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* xbase=&x[block][0][0][0][0];
      FT* ybase = &y[block][0][0][0][0];
      
      // Temporary storage to stream into and out of
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor;
#endif

      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)y_spinor, ybase, nvec_in_spinor);
      
      // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN
            for(int i=0; i < S; i++) {
              y_spinor[col][spin][reim][i] = a*y_spinor[col][spin][reim][i] + x_spinor[col][spin][reim][i];
            }
#else
            y_spinor[col][spin][reim][:] = a*y_spinor[col][spin][reim][:] + x_spinor[col][spin][reim][:];
#endif
	  }
	}
      }
      
      BLASUtils::streamOutSpinor<FT,V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
    }
    
  private: 
    AT a;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y;
  };


  template<typename FT, int V, int S, bool compress>
  class Norm2Functor {
  public:
    typedef typename ArithType<FT>::Type AT;

    Norm2Functor(const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_) : x(x_) {}
    ~Norm2Functor() {}
    
    inline void 
      func(int block, double* reduction)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* xbase=&x[block][0][0][0][0];
      
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
#endif

      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
	    for(int s=0; s < S; s++) {
              double xfoo = rep<double,AT>(x_spinor[col][spin][reim][s]);
	      reduction[s] += xfoo*xfoo;
	    }
	  }
	}
     }

 
    }
    
  private: 
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x;
  };
  

  template<typename FT, int V, int S, bool compress>
  class XMYNorm2Functor {
  public:
    typedef typename ArithType<FT>::Type AT;
    
    XMYNorm2Functor(typename Geometry<FT,V,S,compress>::FourSpinorBlock* res_,
		    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_,
		    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* y_
		    ) : res(res_), x(x_), y(y_) {}
    
    ~XMYNorm2Functor() {}
    
    inline 
    void 
    func(int block, double* reduction) {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* xbase=&x[block][0][0][0][0];
      const FT* ybase=&y[block][0][0][0][0];
      FT* resbase=&res[block][0][0][0][0];
      
      // Temporary storage to stream into and out of
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock res_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock res_spinor;
#endif      
      
      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)y_spinor, ybase, nvec_in_spinor);
      
      // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN 
            for(int i = 0; i < S; i++)
              res_spinor[col][spin][reim][i] = x_spinor[col][spin][reim][i] -  y_spinor[col][spin][reim][i];
#else
            res_spinor[col][spin][reim][:] = x_spinor[col][spin][reim][:] -  y_spinor[col][spin][reim][:];
#endif
	    for(int s=0; s < S; s++) {
	      reduction[s] += (double)res_spinor[col][spin][reim][s]*(double)res_spinor[col][spin][reim][s];
	    }
	  }
	}
      }
      BLASUtils::streamOutSpinor<FT,V>(resbase, (const AT *)res_spinor, nvec_in_spinor);
    }
    
  private: 
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y;
  };
  
  
  template<typename FT, int V, int S, bool compress>
  class XMYFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;
    
    XMYFunctor(const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_,
	       typename Geometry<FT,V,S,compress>::FourSpinorBlock* y_
	       ) :  x(x_), y(y_) {}
    
    ~XMYFunctor() {}
    
    inline 
    void 
    func(int block) {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* xbase=&x[block][0][0][0][0];
      FT* ybase=&y[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock y_spinor;
#endif
 
      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)y_spinor, ybase, nvec_in_spinor);
      // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN
            for(int i = 0; i < S; i++)
              y_spinor[col][spin][reim][i] = x_spinor[col][spin][reim][i] -  y_spinor[col][spin][reim][i];
#else
            y_spinor[col][spin][reim][:] = x_spinor[col][spin][reim][:] -  y_spinor[col][spin][reim][:];
#endif
	  }
	}
      }
      BLASUtils::streamOutSpinor<FT,V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
    }
   
 private: 
   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
   typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y;
 };
  

  template<typename FT, int V, int S, bool compress>
  class   RmammpNorm2rxpapFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;

    RmammpNorm2rxpapFunctor(double a_,
			    typename Geometry<FT,V,S,compress>::FourSpinorBlock* r_,
			    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* mmp_,
			    typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_,
			    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* p_
			    ) : a(rep<AT,double>(a_)), r(r_), mmp(mmp_), x(x_), p(p_) {}
    
    ~RmammpNorm2rxpapFunctor() {}
    
    inline 
    void 
    func(int block, double*  reduction) {
      int nvec_in_spinor = (3*4*2*S)/V;
      FT* rbase=&r[block][0][0][0][0];
      const FT* mmpbase=&mmp[block][0][0][0][0];
      FT* xbase=&x[block][0][0][0][0];
      const FT* pbase=&p[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock mmp_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN))); 
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock mmp_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor;
#endif
 
      BLASUtils::streamInSpinor<FT,V>((AT *)r_spinor, rbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)mmp_spinor, mmpbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)p_spinor, pbase, nvec_in_spinor);
      
      // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN
            for(int i = 0; i < S; i++)
              r_spinor[col][spin][reim][i] = r_spinor[col][spin][reim][i] -  a * mmp_spinor[col][spin][reim][i];
#else
            r_spinor[col][spin][reim][:] = r_spinor[col][spin][reim][:] -  a * mmp_spinor[col][spin][reim][:];
#endif
	    for(int s =0 ; s < S; s++) { 
	      reduction[s] += (double)r_spinor[col][spin][reim][s]*(double)r_spinor[col][spin][reim][s];
	    }

#ifndef QPHIX_USE_CEAN
            for(int i = 0; i < S; i++)
              x_spinor[col][spin][reim][i] = x_spinor[col][spin][reim][i] + a * p_spinor[col][spin][reim][i];
#else
            x_spinor[col][spin][reim][:] = x_spinor[col][spin][reim][:] + a * p_spinor[col][spin][reim][:];
#endif
	  }
	}
      }
	
      BLASUtils::streamOutSpinor<FT,V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
      BLASUtils::streamOutSpinor<FT,V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
    }


private: 
  AT a;
  typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r;
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict mmp;
  typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p;
};


  template<typename FT, int V, int S, bool compress>
    class   RichardsonRXUpdateNormRFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;

    RichardsonRXUpdateNormRFunctor(
				   typename Geometry<FT,V,S,compress>::FourSpinorBlock* x_,
				   typename Geometry<FT,V,S,compress>::FourSpinorBlock* r_,

				   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* delta_x_,
				   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* delta_r_
				   ) : x(x_), r(r_), delta_x(delta_x_), delta_r(delta_r_) {}
    
    ~RichardsonRXUpdateNormRFunctor() {}
    
    inline 
    void 
    func(int block, double*  reduction) {
      FT* xbase=&x[block][0][0][0][0];
      FT* rbase=&r[block][0][0][0][0];
      const FT* delta_xbase=&delta_x[block][0][0][0][0];
      const FT* delta_rbase=&delta_r[block][0][0][0][0];
      int nvec_in_spinor = (3*4*2*S)/V;

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock delta_x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock delta_r_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock delta_x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock delta_r_spinor;
#endif
 
      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)r_spinor, rbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)delta_x_spinor, delta_xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)delta_r_spinor, delta_rbase, nvec_in_spinor);
      
      // Now we are hopefully both in L1 and in the right layout so
      for(int col=0; col < 3; col++) { 
	for(int spin=0; spin < 4; spin ++) { 
	  for(int reim=0; reim < 2; reim++) { 
#ifndef QPHIX_USE_CEAN
            for(int i = 0; i < S; i++) {
              x_spinor[col][spin][reim][i] += delta_x_spinor[col][spin][reim][i];
              r_spinor[col][spin][reim][i] -= delta_r_spinor[col][spin][reim][i];
            }
#else
            x_spinor[col][spin][reim][:] += delta_x_spinor[col][spin][reim][:];
            r_spinor[col][spin][reim][:] -= delta_r_spinor[col][spin][reim][:];
#endif
	    for(int s =0 ; s < S; s++) { 
	      reduction[s] += (double)r_spinor[col][spin][reim][s]*(double)r_spinor[col][spin][reim][s];
	    }
	  }
	}
      }

      BLASUtils::writeSpinor<FT,V>(xbase, (const AT *)x_spinor, nvec_in_spinor);	
      BLASUtils::writeSpinor<FT,V>(rbase, (const AT *)r_spinor, nvec_in_spinor);

    }


private: 
  typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
  typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r;
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict delta_x;
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict delta_r;
};


}; // Namespace

#endif
