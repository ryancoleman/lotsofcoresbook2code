#ifndef QPHIX_COMPLEX_FUNCTORS_H
#define QPHIX_COMPLEX_FUNCTORS_H

#include "qphix/blas_utils.h"

namespace QPhiX  { 


  template<typename FT, int V, int S, bool compress>
  class InnerProductFunctor {
  public:
    
    typedef typename ArithType<FT>::Type AT;

    InnerProductFunctor( const typename Geometry<FT,V,S,compress>::FourSpinorBlock* l_,
			 const typename Geometry<FT,V,S,compress>::FourSpinorBlock* r_ ) : l(l_), r(r_) {}
    
    ~InnerProductFunctor() {}
    
    inline void 
    func(int block, double* reduction_re, double *reduction_im)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* lbase=&l[block][0][0][0][0];
      const FT* rbase=&r[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock l_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock l_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor;
#endif

      BLASUtils::streamInSpinor<FT,V>((AT *)l_spinor, lbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)r_spinor, rbase, nvec_in_spinor);

      // Accumulate the inner product from this spinor

      for(int col=0; col < 3; col++) {
	for(int spin=0; spin < 4; spin++) { 
	  for(int s = 0; s < S; s++) { 
	    reduction_re[s] += l_spinor[col][spin][0][s]*r_spinor[col][spin][0][s]
	      + l_spinor[col][spin][1][s]*r_spinor[col][spin][1][s];
      	    reduction_im[s] += l_spinor[col][spin][0][s]*r_spinor[col][spin][1][s]
	      - l_spinor[col][spin][1][s]*r_spinor[col][spin][0][s];
	  }
	}
      }
    }
    
  private: 
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict l;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r;
  };
  


  template<typename FT, int V, int S, bool compress>
  class BiCGStabPUpdateFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;

    BiCGStabPUpdateFunctor(const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r_,
			   typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p_,
			   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict v_,
			   double beta_[2], double omega_[2]) : r(r_), p(p_), v(v_)
    {
      beta[0]=rep<AT,double>(beta_[0]);
      beta[1]=rep<AT,double>(beta_[1]);
      omega[0]=rep<AT,double>(omega_[0]);
      omega[1]=rep<AT,double>(omega_[1]);
    }


		     
    ~BiCGStabPUpdateFunctor() {}
    
    inline void 
    func(int block)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      const FT* rbase=&r[block][0][0][0][0];
      FT* pbase=&p[block][0][0][0][0];
      const FT* vbase=&v[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor  __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor  __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock v_spinor  __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock v_spinor;
#endif
      
      BLASUtils::streamInSpinor<FT,V>((AT *)r_spinor, rbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)p_spinor, pbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)v_spinor, vbase, nvec_in_spinor);
      
      for(int col=0; col < 3; col++){ 
	for(int spin=0; spin < 4; spin++) { 

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
          AT tmp_cmpx[2][S] __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
          __declspec(align(QPHIX_LLC_CACHE_ALIGN)) AT tmp_cmpx[2][S];
#endif

	  // Evaluate: p = r + beta (p - omega v)
	  //           p = r + beta tmp
	  // with tmp = p - omega
	  
	  // tmp = -omega v + p = p - omega v
	  BLASUtils::cnmadd<AT,S>(tmp_cmpx, omega, v_spinor[col][spin], p_spinor[col][spin]);
	  
	  // p = r + beta tmp 
	  BLASUtils::cmadd<AT,S>(p_spinor[col][spin], beta, tmp_cmpx, r_spinor[col][spin]);
	}
      }

      // Should this be a stream out spinor?
      BLASUtils::writeSpinor<FT,V>(pbase, (const AT *)p_spinor, nvec_in_spinor);

    }
  
private: 
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r;
  typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p;
  const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict v;
  AT beta[2];
  AT omega[2];
  
};

  template<typename FT, int V, int S, bool compress>
  class BiCGStabSUpdateFunctor {
  public:

    typedef typename ArithType<FT>::Type AT;

    BiCGStabSUpdateFunctor(double alpha_[2],
			   typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict s_,
			   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict  v_)
			   : s(s_), v(v_)
    {
      alpha[0]=rep<AT,double>(alpha_[0]);
      alpha[1]=rep<AT,double>(alpha_[1]);
    }
		     
    ~BiCGStabSUpdateFunctor() {}
    
    inline void 
    func(int block)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      FT* sbase=&s[block][0][0][0][0];
      const FT* vbase=&v[block][0][0][0][0];

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock s_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock v_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock s_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock v_spinor;
#endif
      
      BLASUtils::streamInSpinor<FT,V>((AT *)s_spinor, sbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)v_spinor, vbase, nvec_in_spinor);
      
      for(int col=0; col < 3; col++){ 
	for(int spin=0; spin < 4; spin++) { 
	  // s = s - alpha v
	  //   = -alpha v + s
	  // 
	  // 's' is not involved in the complex mul, so its components
	  // are not mixed, so I will leave it as an alias to the cnmadd
	  BLASUtils::cnmadd<AT,S>( s_spinor[col][spin],
				   alpha,
				   v_spinor[col][spin], 
				   s_spinor[col][spin]);


	}
      }
      BLASUtils::writeSpinor<FT,V>(sbase, (const AT *)s_spinor, nvec_in_spinor);
    }
    
  private: 
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict s;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict v;
    AT alpha[2];
  };


  template<typename FT, int V, int S, bool compress>
  class BiCGStabRXUpdateFunctor {
  public:
    typedef typename ArithType<FT>::Type AT;
    BiCGStabRXUpdateFunctor(
			    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x_,
			    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r_,
			    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict t_,
			    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p_,
			    double omega_[2],
			    double alpha_[2]
			    ) : x(x_),r(r_),t(t_),p(p_)
    {
      omega[0]=rep<AT,double>(omega_[0]);
      omega[1]=rep<AT,double>(omega_[1]);
      
      alpha[0]=rep<AT,double>(alpha_[0]);
      alpha[1]=rep<AT,double>(alpha_[1]);
    }
		     
    ~BiCGStabRXUpdateFunctor() {}
    
    inline void 
    func(int block, double* reduction)
    {
      int nvec_in_spinor = (3*4*2*S)/V;
      FT* xbase=&x[block][0][0][0][0];
      FT* rbase=&r[block][0][0][0][0];
      const FT* tbase=&t[block][0][0][0][0];
      const FT* pbase=&p[block][0][0][0][0];
      
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock t_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
      typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock x_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock r_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock t_spinor;
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) typename Geometry<AT,V,S,compress>::FourSpinorBlock p_spinor;
#endif
      
      BLASUtils::streamInSpinor<FT,V>((AT *)x_spinor, xbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)r_spinor, rbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)t_spinor, tbase, nvec_in_spinor);
      BLASUtils::streamInSpinor<FT,V>((AT *)p_spinor, pbase, nvec_in_spinor);
      
      

      for(int col=0; col < 3; col++){ 
	for(int spin=0; spin < 4; spin++) { 
	  AT tmp_cmpx[2][S];

	  /* tmp = alpha p + x */
	  BLASUtils::cmadd<AT,S>(tmp_cmpx,alpha,p_spinor[col][spin],x_spinor[col][spin]);
	  /* x = omega r + tmp = omega r + alpha p + x */ 
	  BLASUtils::cmadd<AT,S>(x_spinor[col][spin], omega, r_spinor[col][spin],tmp_cmpx);
	  
	  /* r = -omega t + r = r - omega t */
	  BLASUtils::cnmadd<AT,S>(r_spinor[col][spin],omega,t_spinor[col][spin],r_spinor[col][spin]);

	  /* accumulate new r_norm into reduction */
	  for(int cmpx=0; cmpx<2; cmpx++) { 
	    for(int s=0; s < S; s++) { 
	      reduction[s] += (double)r_spinor[col][spin][cmpx][s]*(double)r_spinor[col][spin][cmpx][s];
	    }
	  }
	  

	}
      }
      /* Write back x and r */
      /* Should these be streamouts */
      BLASUtils::writeSpinor<FT,V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
      BLASUtils::writeSpinor<FT,V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
    }
			   
  private:
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x;
    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict t;
    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p;
    AT alpha[2];
    AT omega[2];
  };
  

}; // Namespace

#endif
