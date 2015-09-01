// -*- C++ -*-
#ifndef __clover_term_qdp_w_h__
#define __clover_term_qdp_w_h__

#include "clover_fermact_params_w.h"
#include "mesfield.h"

//! Special structure used for triangular objects
template<typename R>
struct PrimitiveClovTriang
{
  RScalar<R>   diag[2][2*Nc];
  RComplex<R>  offd[2][2*Nc*Nc-Nc];
};


//! Clover term
/*!
 * \ingroup linop
 *
 */
template<typename T, typename U>
class QDPCloverTermT  {
public:

  // Typedefs to save typing
  typedef typename WordType<T>::Type_t REALT;
  
  typedef OLattice< PScalar< PScalar< RScalar< typename WordType<T>::Type_t> > > > LatticeREAL;
  typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;
  
  //! Empty constructor. Must use create later
  QDPCloverTermT();

  //! No real need for cleanup here
  ~QDPCloverTermT() {}
  
    //! Creation routine
  void create(const multi1d<U>& u_,
	      const CloverFermActParams& param_);

  void create(const multi1d<U>& u_,
	      const CloverFermActParams& param_,
	      const QDPCloverTermT<T,U>& from_);

  //! Computes the inverse of the term on cb using Cholesky
  /*!
   * \param cb   checkerboard of work (Read)
   */
  void choles(int cb);

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  void apply (T& chi, const T& psi, int isign, int cb) const;


  void applySite(T& chi, const T& psi, int isign, int site) const;
  const multi1d<PrimitiveClovTriang<REALT> >& getTriBuffer() const {
    return tri;
  }

protected:
  //! Create the clover term on cb
  /*!
   *  \param f         field strength tensor F(mu,nu)        (Read)
   *  \param cb        checkerboard                          (Read)
   */
  void makeClov(const multi1d<U>& f, const RealT& diag_mass);
  
  //! Invert the clover term on cb
  void ldagdlinv(LatticeREAL& tr_log_diag, int cb);
  
  //! Get the u field
  const multi1d<U>& getU() const {return u;}

  //! Calculates Tr_D ( Gamma_mat L )
  Real getCloverCoeff(int mu, int nu) const;


private:
  multi1d<U>  u;
  CloverFermActParams          param;
  LatticeREAL                  tr_log_diag_; // Fill this out during create
                                                  // but save the global sum until needed.
  multi1d<bool> choles_done;   // Keep note of whether the decomposition has been done
                                 // on a particular checkerboard. 

  multi1d<PrimitiveClovTriang<REALT> >  tri;
    
};


// Empty constructor. Must use create later
template<typename T, typename U>
QDPCloverTermT<T,U>::QDPCloverTermT() {}

// Now copy
template<typename T, typename U>
void QDPCloverTermT<T,U>::create(const multi1d<U>& u_,
				 const CloverFermActParams& param_,
				 const QDPCloverTermT<T,U>& from)
{
    
    u.resize(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      u[mu] = u_[mu];
    }
    param = param_;
    
    // Apply anisotropy 
    {
      RealT ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= Real(0.5) * ff;
      param.clovCoeffT *= Real(0.5);
    }
    
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }
    
        
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++) {
      choles_done[i] = from.choles_done[i];
    }
    
    tr_log_diag_ = from.tr_log_diag_;
    
    tri = from.tri;
      
  }


  //! Creation routine
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::create(const multi1d<U>& u_,
				   const CloverFermActParams& param_)
  {
    
   
    u.resize(Nd);
    u=u_;
    param = param_;

    {
      RealT ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= RealT(0.5) * ff;
      param.clovCoeffT *= RealT(0.5);
    }
    
    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }
    
    /* Calculate F(mu,nu) */
    multi1d<U> f;
    mesField(f, u);
    makeClov(f, diag_mass);
    
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++) {
      choles_done[i] = false;
    }
    
  }


  /*
   * MAKCLOV 
   *
   *  In this routine, MAKCLOV calculates

   *    1 - (1/4)*sigma(mu,nu) F(mu,nu)

   *  using F from mesfield

   *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]

   *  using basis of SPPROD and stores in a lower triangular matrix
   *  (no diagonal) plus real diagonal

   *  where
   *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
   *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
   *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
   *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

   *  and

   *         | sigF(1)   sigF(3)     0         0     |
   *  sigF = | sigF(5)  -sigF(1)     0         0     |
   *         |   0         0      -sigF(0)  -sigF(2) |
   *         |   0         0      -sigF(4)   sigF(0) |
   *  where
   *    sigF(i) is a color matrix

   *  sigF(0) = i*(ClovT*E_z + ClovR*B_z)
   *          = i*(ClovT*F(3,2) + ClovR*F(1,0))
   *  sigF(1) = i*(ClovT*E_z - ClovR*B_z)
   *          = i*(ClovT*F(3,2) - ClovR*F(1,0))
   *  sigF(2) = i*(E_+ + B_+)
   *  sigF(3) = i*(E_+ - B_+)
   *  sigF(4) = i*(E_- + B_-)
   *  sigF(5) = i*(E_- - B_-)
   *  i*E_+ = (i*ClovT*E_x - ClovT*E_y)
   *        = (i*ClovT*F(3,0) - ClovT*F(3,1))
   *  i*E_- = (i*ClovT*E_x + ClovT*E_y)
   *        = (i*ClovT*F(3,0) + ClovT*F(3,1))
   *  i*B_+ = (i*ClovR*B_x - ClovR*B_y)
   *        = (i*ClovR*F(2,1) + ClovR*F(2,0))
   *  i*B_- = (i*ClovR*B_x + ClovR*B_y)
   *        = (i*ClovR*F(2,1) - ClovR*F(2,0))

   *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

   *  NOTE: the above definitions assume that the time direction, t_dir,
   *        is 3. In general F(k,j) is multiplied with ClovT if either
   *        k=t_dir or j=t_dir, and with ClovR otherwise.

   *+++
   *  Here are some notes on the origin of this routine. NOTE, ClovCoeff or u0
   *  are not actually used in MAKCLOV.
   *
   *  The clover mass term is suppose to act on a vector like
   *
   *  chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi

   *  Definitions used here (NOTE: no "i")
   *   sigma(mu,nu) = gamma(mu)*gamma(nu) - gamma(nu)*gamma(mu)
   *                = 2*gamma(mu)*gamma(nu)   for mu != nu
   *       
   *   chi = sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu < nu
   *       = (1/2) * sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu != nu
   *       = (1/4) * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)*psi
   *
   *
   * chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi
   *     = psi - (ClovCoeff/u0^3) * kappa * chi
   *     == psi - kappa * chi
   *
   *  We have absorbed ClovCoeff/u0^3 into kappa. A u0 was previously absorbed into kappa
   *  for compatibility to ancient conventions. 
   *---

   * Arguments:
   *  \param f         field strength tensor F(cb,mu,nu)        (Read)
   *  \param diag_mass effective mass term                      (Read)
   */

  /*Threaded this. It needs a QMT arg struct and I've extracted the site loop */
  namespace QDPCloverEnv { 
    
    template<typename U>
    struct QDPCloverMakeClovArg {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;
      const RealT& diag_mass;
      const U& f0;
      const U& f1;
      const U& f2;
      const U& f3;
      const U& f4;
      const U& f5;
      multi1d< PrimitiveClovTriang < REALT > >& tri;
    };
    
    
    /* This is the extracted site loop for makeClover */
    template<typename U>
    inline 
    void makeClovSiteLoop(int lo, int hi, int myId, QDPCloverMakeClovArg<U> *a)
    {
      typedef typename QDPCloverMakeClovArg<U>::RealT RealT;
      typedef typename QDPCloverMakeClovArg<U>::REALT REALT;
      
      const RealT& diag_mass = a->diag_mass;
      const U& f0=a->f0;
      const U& f1=a->f1;
      const U& f2=a->f2;
      const U& f3=a->f3;
      const U& f4=a->f4;
      const U& f5=a->f5;
      multi1d<PrimitiveClovTriang < REALT > >& tri=a->tri;

      // SITE LOOP STARTS HERE
      for(int site = lo; site < hi; ++site)  {
	/*# Construct diagonal */
	
	for(int jj = 0; jj < 2; jj++) {
	  
	  for(int ii = 0; ii < 2*Nc; ii++) {
	    
	    tri[site].diag[jj][ii] = diag_mass.elem().elem().elem();
	  }
	}
	
       

	RComplex<REALT> E_minus;
	RComplex<REALT> B_minus;
	RComplex<REALT> ctmp_0;
	RComplex<REALT> ctmp_1;
	RScalar<REALT> rtmp_0;
	RScalar<REALT> rtmp_1;
	
	for(int i = 0; i < Nc; ++i) {
	  
	  /*# diag_L(i,0) = 1 - i*diag(E_z - B_z) */
	  /*#             = 1 - i*diag(F(3,2) - F(1,0)) */
	  ctmp_0 = f5.elem(site).elem().elem(i,i);
	  ctmp_0 -= f0.elem(site).elem().elem(i,i);
	  rtmp_0 = imag(ctmp_0);
	  tri[site].diag[0][i] += rtmp_0;
	  
	  /*# diag_L(i+Nc,0) = 1 + i*diag(E_z - B_z) */
	  /*#                = 1 + i*diag(F(3,2) - F(1,0)) */
	  tri[site].diag[0][i+Nc] -= rtmp_0;
	  
	  /*# diag_L(i,1) = 1 + i*diag(E_z + B_z) */
	  /*#             = 1 + i*diag(F(3,2) + F(1,0)) */
	  ctmp_1 = f5.elem(site).elem().elem(i,i);
	  ctmp_1 += f0.elem(site).elem().elem(i,i);
	  rtmp_1 = imag(ctmp_1);
	  tri[site].diag[1][i] -= rtmp_1;
	  
	  /*# diag_L(i+Nc,1) = 1 - i*diag(E_z + B_z) */
	  /*#                = 1 - i*diag(F(3,2) + F(1,0)) */
	  tri[site].diag[1][i+Nc] += rtmp_1;
	}
	
	/*# Construct lower triangular portion */
	/*# Block diagonal terms */
	for(int i = 1; i < Nc; ++i) {
	  
	  for(int j = 0; j < i; ++j) {
	    
	    int elem_ij  = i*(i-1)/2 + j;
	    int elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;
	    
	    /*# L(i,j,0) = -i*(E_z - B_z)[i,j] */
	    /*#          = -i*(F(3,2) - F(1,0)) */
	    ctmp_0 = f0.elem(site).elem().elem(i,j);
	    ctmp_0 -= f5.elem(site).elem().elem(i,j);
	    tri[site].offd[0][elem_ij] = timesI(ctmp_0);
	    
	    /*# L(i+Nc,j+Nc,0) = +i*(E_z - B_z)[i,j] */
	    /*#                = +i*(F(3,2) - F(1,0)) */
	    tri[site].offd[0][elem_tmp] = -tri[site].offd[0][elem_ij];
	    
	    /*# L(i,j,1) = i*(E_z + B_z)[i,j] */
	    /*#          = i*(F(3,2) + F(1,0)) */
	    ctmp_1 = f5.elem(site).elem().elem(i,j);
	    ctmp_1 += f0.elem(site).elem().elem(i,j);
	    tri[site].offd[1][elem_ij] = timesI(ctmp_1);
	    
	    /*# L(i+Nc,j+Nc,1) = -i*(E_z + B_z)[i,j] */
	    /*#                = -i*(F(3,2) + F(1,0)) */
	    tri[site].offd[1][elem_tmp] = -tri[site].offd[1][elem_ij];
	  }
	}
	
	/*# Off-diagonal */
	for(int i = 0; i < Nc; ++i) {
	  
	  for(int j = 0; j < Nc; ++j) {
	    
	    // Flipped index
	    // by swapping i <-> j. In the past i would run slow
	    // and now j runs slow
	    int elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;
	    
	    /*# i*E_- = (i*E_x + E_y) */
	    /*#       = (i*F(3,0) + F(3,1)) */
	    E_minus = timesI(f2.elem(site).elem().elem(i,j));
	    E_minus += f4.elem(site).elem().elem(i,j);
	    
	    /*# i*B_- = (i*B_x + B_y) */
	    /*#       = (i*F(2,1) - F(2,0)) */
	    B_minus = timesI(f3.elem(site).elem().elem(i,j));
	    B_minus -= f1.elem(site).elem().elem(i,j);
	    
	    /*# L(i+Nc,j,0) = -i*(E_- - B_-)  */
	    tri[site].offd[0][elem_ij] = B_minus - E_minus;
	    
	    /*# L(i+Nc,j,1) = +i*(E_- + B_-)  */
	    tri[site].offd[1][elem_ij] = E_minus + B_minus;
	  }
	}
      } /* End Site loop */
    } /* Function */
  }
  
  /* This now just sets up and dispatches... */
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::makeClov(const multi1d<U>& f, const RealT& diag_mass)
  {
    
    
    if ( Nd != 4 ){
      QDPIO::cerr << __func__ << ": expecting Nd==4" << endl;
      QDP_abort(1);
    }
    
    if ( Ns != 4 ){
      QDPIO::cerr << __func__ << ": expecting Ns==4" << endl;
      QDP_abort(1);
    }
  
    U f0 = f[0] * getCloverCoeff(0,1);
    U f1 = f[1] * getCloverCoeff(0,2);
    U f2 = f[2] * getCloverCoeff(0,3);
    U f3 = f[3] * getCloverCoeff(1,2);
    U f4 = f[4] * getCloverCoeff(1,3);
    U f5 = f[5] * getCloverCoeff(2,3);    

    const int nodeSites = QDP::Layout::sitesOnNode();

    tri.resize(nodeSites);  // hold local lattice

    QDPCloverEnv::QDPCloverMakeClovArg<U> arg = {diag_mass, f0,f1,f2,f3,f4,f5,tri };
    dispatch_to_threads(nodeSites, arg, QDPCloverEnv::makeClovSiteLoop<U>);
              

    
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::choles(int cb)
  {
    

    // When you are doing the cholesky - also fill out the trace_log_diag piece)
    // chlclovms(tr_log_diag_, cb);
    // Switch to LDL^\dag inversion
    ldagdlinv(tr_log_diag_,cb);

    
  }


  namespace QDPCloverEnv { 
    template<typename U>
    struct LDagDLInvArgs { 
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;
      typedef OLattice< PScalar< PScalar< RScalar<REALT> > > > LatticeRealT;
      LatticeRealT& tr_log_diag;
      multi1d<PrimitiveClovTriang<REALT> >& tri;
      int cb;
    };

    template<typename U>
    inline 
    void LDagDLInvSiteLoop(int lo, int hi, int myId, LDagDLInvArgs<U>* a) 
    {
      typedef typename LDagDLInvArgs<U>::REALT REALT;
      typedef typename LDagDLInvArgs<U>::RealT RealT;
      typedef typename LDagDLInvArgs<U>::LatticeRealT LatticeRealT;

      LatticeRealT& tr_log_diag = a->tr_log_diag;
      multi1d<PrimitiveClovTriang < REALT> >& tri = a->tri;
      int cb = a->cb;
      
      
      RScalar<REALT> zip=0;
      int N = 2*Nc;
      
      // Loop through the sites.
      for(int ssite=lo; ssite < hi; ++ssite)  {

	int site = rb[cb].siteTable()[ssite];
	
	int site_neg_logdet=0;
	// Loop through the blocks on the site.
	for(int block=0; block < 2; block++) { 
	  
	  // Triangular storage 
	  RScalar<REALT> inv_d[6] QDP_ALIGN16;
	  RComplex<REALT> inv_offd[15] QDP_ALIGN16;
	  RComplex<REALT> v[6] QDP_ALIGN16;
	  RScalar<REALT>  diag_g[6] QDP_ALIGN16;
	  // Algorithm 4.1.2 LDL^\dagger Decomposition
	  // From Golub, van Loan 3rd ed, page 139
	  for(int i=0; i < N; i++) { 
	    inv_d[i] = tri[site].diag[block][i];
	  }
	  
	  for(int i=0; i < 15; i++) { 
	    inv_offd[i]  =tri[site].offd[block][i];
	  }
	  
	  for(int j=0; j < N; ++j) { 
	    
	    // Compute v(0:j-1)
	    //
	    // for i=0:j-2
	    //   v(i) = A(j,i) A(i,i)
	    // end
	    
	    
	    for(int i=0; i < j; i++) { 
	      int elem_ji = j*(j-1)/2 + i;
	      
	      RComplex<REALT> A_ii = cmplx( inv_d[i], zip );
	      v[i] = A_ii*adj(inv_offd[elem_ji]);
	    }
	    
	    // v(j) = A(j,j) - A(j, 0:j-2) v(0:j-2)
	    //                 ^ This is done with a loop over k ie:
	    //
	    // v(j) = A(j,j) - sum_k A*(j,k) v(k)     k=0...j-2
	    //
	    //      = A(j,j) - sum_k A*(j,k) A(j,k) A(k,k)
	    //      = A(j,j) - sum_k | A(j,k) |^2 A(k,k)
	    
	    v[j] = cmplx(inv_d[j],zip);
	    
	    for(int k=0; k < j; k++) { 
	      int elem_jk = j*(j-1)/2 + k;
	      v[j] -= inv_offd[elem_jk]*v[k];
	    }
	    
	    
	    // At this point in time v[j] has to be real, since
	    // A(j,j) is from diag ie real and all | A(j,k) |^2 is real
	    // as is A(k,k)
	    
	    // A(j,j) is the diagonal element - so store it.
	    inv_d[j] = real( v[j] );
	    
	    // Last line of algorithm:
	    // A( j+1 : n, j) = ( A(j+1:n, j) - A(j+1:n, 1:j-1)v(1:k-1) ) / v(j)
	    //
	    // use k as first colon notation and l as second so
	    // 
	    // for k=j+1 < n-1
	    //      A(k,j) = A(k,j) ;
	    //      for l=0 < j-1
	    //         A(k,j) -= A(k, l) v(l)
	    //      end
	    //      A(k,j) /= v(j);
	    //
	    for(int k=j+1; k < N; k++) { 
	      int elem_kj = k*(k-1)/2 + j;
	      for(int l=0; l < j; l++) { 
		int elem_kl = k*(k-1)/2 + l;
		inv_offd[elem_kj] -= inv_offd[elem_kl] * v[l];
	      }
	      inv_offd[elem_kj] /= v[j];
	    }
	  }
	  
	  // Now fix up the inverse
	  RScalar<REALT> one;
	  one.elem() = (REALT)1;
	  
	  for(int i=0; i < N; i++) { 
	    diag_g[i] = one/inv_d[i];
	    
	    // Compute the trace log
	    // NB we are always doing trace log | A | 
	    // (because we are always working with actually A^\dagger A
	    //  even in one flavour case where we square root)
	    tr_log_diag.elem(site).elem().elem().elem() += log(fabs(inv_d[i].elem()));
	    // However, it is worth counting just the no of negative logdets
	    // on site
	    if( inv_d[i].elem() < 0 ) { 
	      site_neg_logdet++;
	    }
	  }
	  // Now we need to invert the L D L^\dagger 
	  // We can do this by solving:
	  //
	  //  L D L^\dagger M^{-1} = 1   
	  //
	  // This can be done by solving L D X = 1  (X = L^\dagger M^{-1})
	  //
	  // Then solving L^\dagger M^{-1} = X
	  //
	  // LD is lower diagonal and so X will also be lower diagonal.
	  // LD X = 1 can be solved by forward substitution.
	  //
	  // Likewise L^\dagger is strictly upper triagonal and so
	  // L^\dagger M^{-1} = X can be solved by forward substitution.
	  RComplex<REALT> sum;
	  for(int k = 0; k < N; ++k) {

	    for(int i = 0; i < k; ++i) {
	      zero_rep(v[i]);
	    }
	    
	    /*# Forward substitution */
	    
	    // The first element is the inverse of the diagonal
	    v[k] = cmplx(diag_g[k],zip);
	    
	    for(int i = k+1; i < N; ++i) {
	      zero_rep(v[i]);
	      
	      for(int j = k; j < i; ++j) {
		int elem_ij = i*(i-1)/2+j;	
		
		// subtract l_ij*d_j*x_{kj}
		v[i] -= inv_offd[elem_ij] *inv_d[j]*v[j];
		
	      }
	      
	      // scale out by 1/d_i
	      v[i] *= diag_g[i];
	    }
	    
	    /*# Backward substitution */
	    // V[N-1] remains unchanged
	    // Start from V[N-2]
	    
	    for(int i = N-2; (int)i >= (int)k; --i) {
	      for(int j = i+1; j < N; ++j) {
		int elem_ji = j*(j-1)/2 + i;
		// Subtract terms of typ (l_ji)*x_kj
		v[i] -= adj(inv_offd[elem_ji]) * v[j];
	      }
	    }
	    
	    /*# Overwrite column k of invcl.offd */
	    inv_d[k] = real(v[k]);
	    for(int i = k+1; i < N; ++i) {

	      int elem_ik = i*(i-1)/2+k;
	      inv_offd[elem_ik] = v[i];
	    }
	  }
	  
	  
	  // Overwrite original data
	  for(int i=0; i < N; i++) { 
	    tri[site].diag[block][i] = inv_d[i];
	  }
	  for(int i=0; i < 15; i++) { 
	    tri[site].offd[block][i] = inv_offd[i];
	  }
	}
	
	if( site_neg_logdet != 0 ) { 
	  // Report if site has any negative terms. (-ve def)
	  std::cout << "WARNING: found " << site_neg_logdet
		    << " negative eigenvalues in Clover DET at site: " << site << endl;
	}
      }/* End Site Loop */
    } /* End Function */
  } /* End Namespace */


  /*! An LDL^\dag decomposition and inversion? */
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::ldagdlinv(LatticeREAL& tr_log_diag, int cb)
  {
    

    if ( 2*Nc < 3 )
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << endl;
      QDP_abort(1);
    }

    // Zero trace log
    tr_log_diag = zero;

    QDPCloverEnv::LDagDLInvArgs<U> a = { tr_log_diag, tri, cb };
    dispatch_to_threads(rb[cb].numSiteTable(), a, QDPCloverEnv::LDagDLInvSiteLoop<U>);

    
    // This comes from the days when we used to do Cholesky
    choles_done[cb] = true;
    
  }
 



  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::applySite(T& chi, const T& psi, 
			    int isign, int site) const
  {
    

    if ( Ns != 4 )
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << endl;
      QDP_abort(1);
    }

    int n = 2*Nc;

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));


    cchi[ 0] = tri[site].diag[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 0]) * ppsi[ 1]
      +   conj(tri[site].offd[0][ 1]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 3]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 6]) * ppsi[ 4]
      +   conj(tri[site].offd[0][10]) * ppsi[ 5];
    
    cchi[ 1] = tri[site].diag[0][ 1]  * ppsi[ 1]
      +        tri[site].offd[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 2]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 4]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 7]) * ppsi[ 4]
      +   conj(tri[site].offd[0][11]) * ppsi[ 5];
    
    cchi[ 2] = tri[site].diag[0][ 2]  * ppsi[ 2]
      +        tri[site].offd[0][ 1]  * ppsi[ 0]
      +        tri[site].offd[0][ 2]  * ppsi[ 1]
      +   conj(tri[site].offd[0][ 5]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 8]) * ppsi[ 4]
      +   conj(tri[site].offd[0][12]) * ppsi[ 5];
    
    cchi[ 3] = tri[site].diag[0][ 3]  * ppsi[ 3]
      +        tri[site].offd[0][ 3]  * ppsi[ 0]
      +        tri[site].offd[0][ 4]  * ppsi[ 1]
      +        tri[site].offd[0][ 5]  * ppsi[ 2]
      +   conj(tri[site].offd[0][ 9]) * ppsi[ 4]
      +   conj(tri[site].offd[0][13]) * ppsi[ 5];
    
    cchi[ 4] = tri[site].diag[0][ 4]  * ppsi[ 4]
      +        tri[site].offd[0][ 6]  * ppsi[ 0]
      +        tri[site].offd[0][ 7]  * ppsi[ 1]
      +        tri[site].offd[0][ 8]  * ppsi[ 2]
      +        tri[site].offd[0][ 9]  * ppsi[ 3]
      +   conj(tri[site].offd[0][14]) * ppsi[ 5];
    
    cchi[ 5] = tri[site].diag[0][ 5]  * ppsi[ 5]
      +        tri[site].offd[0][10]  * ppsi[ 0]
      +        tri[site].offd[0][11]  * ppsi[ 1]
      +        tri[site].offd[0][12]  * ppsi[ 2]
      +        tri[site].offd[0][13]  * ppsi[ 3]
      +        tri[site].offd[0][14]  * ppsi[ 4];
    
    cchi[ 6] = tri[site].diag[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 0]) * ppsi[ 7]
      +   conj(tri[site].offd[1][ 1]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 3]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 6]) * ppsi[10]
      +   conj(tri[site].offd[1][10]) * ppsi[11];
    
    cchi[ 7] = tri[site].diag[1][ 1]  * ppsi[ 7]
      +        tri[site].offd[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 2]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 4]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 7]) * ppsi[10]
      +   conj(tri[site].offd[1][11]) * ppsi[11];
    
    cchi[ 8] = tri[site].diag[1][ 2]  * ppsi[ 8]
      +        tri[site].offd[1][ 1]  * ppsi[ 6]
      +        tri[site].offd[1][ 2]  * ppsi[ 7]
      +   conj(tri[site].offd[1][ 5]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 8]) * ppsi[10]
      +   conj(tri[site].offd[1][12]) * ppsi[11];
    
    cchi[ 9] = tri[site].diag[1][ 3]  * ppsi[ 9]
      +        tri[site].offd[1][ 3]  * ppsi[ 6]
      +        tri[site].offd[1][ 4]  * ppsi[ 7]
      +        tri[site].offd[1][ 5]  * ppsi[ 8]
      +   conj(tri[site].offd[1][ 9]) * ppsi[10]
      +   conj(tri[site].offd[1][13]) * ppsi[11];
    
    cchi[10] = tri[site].diag[1][ 4]  * ppsi[10]
      +        tri[site].offd[1][ 6]  * ppsi[ 6]
      +        tri[site].offd[1][ 7]  * ppsi[ 7]
      +        tri[site].offd[1][ 8]  * ppsi[ 8]
      +        tri[site].offd[1][ 9]  * ppsi[ 9]
      +   conj(tri[site].offd[1][14]) * ppsi[11];
    
    cchi[11] = tri[site].diag[1][ 5]  * ppsi[11]
      +        tri[site].offd[1][10]  * ppsi[ 6]
      +        tri[site].offd[1][11]  * ppsi[ 7]
      +        tri[site].offd[1][12]  * ppsi[ 8]
      +        tri[site].offd[1][13]  * ppsi[ 9]
      +        tri[site].offd[1][14]  * ppsi[10];


    
  }





  //! Returns the appropriate clover coefficient for indices mu and nu
  template<typename T, typename U>
  Real
  QDPCloverTermT<T,U>::getCloverCoeff(int mu, int nu) const 
  { 
    

    if( param.anisoParam.anisoP )  {
      if (mu==param.anisoParam.t_dir || nu == param.anisoParam.t_dir) { 
	return param.clovCoeffT;
      }
      else { 
	// Otherwise return the spatial coeff
	return param.clovCoeffR;
      }
    }
    else { 
      // If there is no anisotropy just return the spatial one, it will
      // be the same as the temporal one
      return param.clovCoeffR; 
    } 
    
    
  }


  namespace QDPCloverEnv { 

    template<typename T>
    struct ApplyArgs {
       typedef typename WordType<T>::Type_t REALT;
      T& chi;
      const T& psi;
      const multi1d<PrimitiveClovTriang<REALT> >& tri;
      int cb;
    };


    template<typename T>
    void applySiteLoop(int lo, int hi, int MyId,
		       ApplyArgs<T>* arg)
      
    {
      
      // This is essentially the body of the previous "Apply"
      // but now the args are handed in through user arg struct...
      
      

      typedef typename WordType<T>::Type_t REALT;
      // Unwrap the args...
      T& chi=arg->chi;
      const T& psi=arg->psi;
      const multi1d<PrimitiveClovTriang<REALT> >& tri = arg->tri;
      int cb = arg->cb;
      

      int n = 2*Nc;
      
      const multi1d<int>& tab = rb[cb].siteTable();
      
      // Now just loop from low to high sites...
      for(int ssite=lo; ssite < hi; ++ssite)  {
	
	int site = tab[ssite];
	
	RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
	const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));
#if 0
#warning "Using unrolled clover term"
      // Rolled version
      for(int i = 0; i < n; ++i)
      {
	cchi[0*n+i] = tri[site].diag[0][i] * ppsi[0*n+i];
	cchi[1*n+i] = tri[site].diag[1][i] * ppsi[1*n+i];
      }

      int kij = 0;  
      for(int i = 0; i < n; ++i)
      {
	for(int j = 0; j < i; j++)
	{
	  cchi[0*n+i] += tri[site].offd[0][kij] * ppsi[0*n+j];
	  cchi[0*n+j] += conj(tri[site].offd[0][kij]) * ppsi[0*n+i];
	  cchi[1*n+i] += tri[site].offd[1][kij] * ppsi[1*n+j];
	  cchi[1*n+j] += conj(tri[site].offd[1][kij]) * ppsi[1*n+i];
	  kij++;
	}
      }
#elif 0
#warning "Using unrolled clover term - version 1"
      // Unrolled version - basically copying the loop structure
      cchi[ 0] = tri[site].diag[0][0] * ppsi[ 0];
      cchi[ 1] = tri[site].diag[0][1] * ppsi[ 1];
      cchi[ 2] = tri[site].diag[0][2] * ppsi[ 2];
      cchi[ 3] = tri[site].diag[0][3] * ppsi[ 3];
      cchi[ 4] = tri[site].diag[0][4] * ppsi[ 4];
      cchi[ 5] = tri[site].diag[0][5] * ppsi[ 5];
      cchi[ 6] = tri[site].diag[1][0] * ppsi[ 6];
      cchi[ 7] = tri[site].diag[1][1] * ppsi[ 7];
      cchi[ 8] = tri[site].diag[1][2] * ppsi[ 8];
      cchi[ 9] = tri[site].diag[1][3] * ppsi[ 9];
      cchi[10] = tri[site].diag[1][4] * ppsi[10];
      cchi[11] = tri[site].diag[1][5] * ppsi[11];

      // cchi[0*n+i] += tri[site].offd[0][kij] * ppsi[0*n+j];
      cchi[ 1] += tri[site].offd[0][ 0] * ppsi[ 0];
      cchi[ 2] += tri[site].offd[0][ 1] * ppsi[ 0];
      cchi[ 2] += tri[site].offd[0][ 2] * ppsi[ 1];
      cchi[ 3] += tri[site].offd[0][ 3] * ppsi[ 0];
      cchi[ 3] += tri[site].offd[0][ 4] * ppsi[ 1];
      cchi[ 3] += tri[site].offd[0][ 5] * ppsi[ 2];
      cchi[ 4] += tri[site].offd[0][ 6] * ppsi[ 0];
      cchi[ 4] += tri[site].offd[0][ 7] * ppsi[ 1];
      cchi[ 4] += tri[site].offd[0][ 8] * ppsi[ 2];
      cchi[ 4] += tri[site].offd[0][ 9] * ppsi[ 3];
      cchi[ 5] += tri[site].offd[0][10] * ppsi[ 0];
      cchi[ 5] += tri[site].offd[0][11] * ppsi[ 1];
      cchi[ 5] += tri[site].offd[0][12] * ppsi[ 2];
      cchi[ 5] += tri[site].offd[0][13] * ppsi[ 3];
      cchi[ 5] += tri[site].offd[0][14] * ppsi[ 4];

      // cchi[0*n+j] += conj(tri[site].offd[0][kij]) * ppsi[0*n+i];
      cchi[ 0] += conj(tri[site].offd[0][ 0]) * ppsi[ 1];
      cchi[ 0] += conj(tri[site].offd[0][ 1]) * ppsi[ 2];
      cchi[ 1] += conj(tri[site].offd[0][ 2]) * ppsi[ 2];
      cchi[ 0] += conj(tri[site].offd[0][ 3]) * ppsi[ 3];
      cchi[ 1] += conj(tri[site].offd[0][ 4]) * ppsi[ 3];
      cchi[ 2] += conj(tri[site].offd[0][ 5]) * ppsi[ 3];
      cchi[ 0] += conj(tri[site].offd[0][ 6]) * ppsi[ 4];
      cchi[ 1] += conj(tri[site].offd[0][ 7]) * ppsi[ 4];
      cchi[ 2] += conj(tri[site].offd[0][ 8]) * ppsi[ 4];
      cchi[ 3] += conj(tri[site].offd[0][ 9]) * ppsi[ 4];
      cchi[ 0] += conj(tri[site].offd[0][10]) * ppsi[ 5];
      cchi[ 1] += conj(tri[site].offd[0][11]) * ppsi[ 5];
      cchi[ 2] += conj(tri[site].offd[0][12]) * ppsi[ 5];
      cchi[ 3] += conj(tri[site].offd[0][13]) * ppsi[ 5];
      cchi[ 4] += conj(tri[site].offd[0][14]) * ppsi[ 5];

      // cchi[1*n+i] += tri[site].offd[1][kij] * ppsi[1*n+j];
      cchi[ 7] += tri[site].offd[1][ 0] * ppsi[ 6];
      cchi[ 8] += tri[site].offd[1][ 1] * ppsi[ 6];
      cchi[ 8] += tri[site].offd[1][ 2] * ppsi[ 7];
      cchi[ 9] += tri[site].offd[1][ 3] * ppsi[ 6];
      cchi[ 9] += tri[site].offd[1][ 4] * ppsi[ 7];
      cchi[ 9] += tri[site].offd[1][ 5] * ppsi[ 8];
      cchi[10] += tri[site].offd[1][ 6] * ppsi[ 6];
      cchi[10] += tri[site].offd[1][ 7] * ppsi[ 7];
      cchi[10] += tri[site].offd[1][ 8] * ppsi[ 8];
      cchi[10] += tri[site].offd[1][ 9] * ppsi[ 9];
      cchi[11] += tri[site].offd[1][10] * ppsi[ 6];
      cchi[11] += tri[site].offd[1][11] * ppsi[ 7];
      cchi[11] += tri[site].offd[1][12] * ppsi[ 8];
      cchi[11] += tri[site].offd[1][13] * ppsi[ 9];
      cchi[11] += tri[site].offd[1][14] * ppsi[10];

      // cchi[1*n+j] += conj(tri[site].offd[1][kij]) * ppsi[1*n+i];
      cchi[ 6] += conj(tri[site].offd[1][ 0]) * ppsi[ 7];
      cchi[ 6] += conj(tri[site].offd[1][ 1]) * ppsi[ 8];
      cchi[ 7] += conj(tri[site].offd[1][ 2]) * ppsi[ 8];
      cchi[ 6] += conj(tri[site].offd[1][ 3]) * ppsi[ 9];
      cchi[ 7] += conj(tri[site].offd[1][ 4]) * ppsi[ 9];
      cchi[ 8] += conj(tri[site].offd[1][ 5]) * ppsi[ 9];
      cchi[ 6] += conj(tri[site].offd[1][ 6]) * ppsi[10];
      cchi[ 7] += conj(tri[site].offd[1][ 7]) * ppsi[10];
      cchi[ 8] += conj(tri[site].offd[1][ 8]) * ppsi[10];
      cchi[ 9] += conj(tri[site].offd[1][ 9]) * ppsi[10];
      cchi[ 6] += conj(tri[site].offd[1][10]) * ppsi[11];
      cchi[ 7] += conj(tri[site].offd[1][11]) * ppsi[11];
      cchi[ 8] += conj(tri[site].offd[1][12]) * ppsi[11];
      cchi[ 9] += conj(tri[site].offd[1][13]) * ppsi[11];
      cchi[10] += conj(tri[site].offd[1][14]) * ppsi[11];
#elif 0
#warning "Using unrolled clover term - version 2"
      // Unrolled version - collect all LHS terms into 1 expression
      cchi[ 0]  =      tri[site].diag[0][ 0]  * ppsi[ 0];
      cchi[ 0] += conj(tri[site].offd[0][ 0]) * ppsi[ 1];
      cchi[ 0] += conj(tri[site].offd[0][ 1]) * ppsi[ 2];
      cchi[ 0] += conj(tri[site].offd[0][ 3]) * ppsi[ 3];
      cchi[ 0] += conj(tri[site].offd[0][ 6]) * ppsi[ 4];
      cchi[ 0] += conj(tri[site].offd[0][10]) * ppsi[ 5];

      cchi[ 1]  = tri[site].diag[0][ 1]  * ppsi[ 1];
      cchi[ 1] += tri[site].offd[0][ 0]  * ppsi[ 0];
      cchi[ 1] += conj(tri[site].offd[0][ 2]) * ppsi[ 2];
      cchi[ 1] += conj(tri[site].offd[0][ 4]) * ppsi[ 3];
      cchi[ 1] += conj(tri[site].offd[0][ 7]) * ppsi[ 4];
      cchi[ 1] += conj(tri[site].offd[0][11]) * ppsi[ 5];

      cchi[ 2]  = tri[site].diag[0][ 2]  * ppsi[ 2];
      cchi[ 2] += tri[site].offd[0][ 1]  * ppsi[ 0];
      cchi[ 2] += tri[site].offd[0][ 2]  * ppsi[ 1];
      cchi[ 2] += conj(tri[site].offd[0][ 5]) * ppsi[ 3];
      cchi[ 2] += conj(tri[site].offd[0][ 8]) * ppsi[ 4];
      cchi[ 2] += conj(tri[site].offd[0][12]) * ppsi[ 5];

      cchi[ 3]  = tri[site].diag[0][ 3]  * ppsi[ 3];
      cchi[ 3] += tri[site].offd[0][ 3]  * ppsi[ 0];
      cchi[ 3] += tri[site].offd[0][ 4]  * ppsi[ 1];
      cchi[ 3] += tri[site].offd[0][ 5]  * ppsi[ 2];
      cchi[ 3] += conj(tri[site].offd[0][ 9]) * ppsi[ 4];
      cchi[ 3] += conj(tri[site].offd[0][13]) * ppsi[ 5];

      cchi[ 4]  = tri[site].diag[0][ 4]  * ppsi[ 4];
      cchi[ 4] += tri[site].offd[0][ 6]  * ppsi[ 0];
      cchi[ 4] += tri[site].offd[0][ 7]  * ppsi[ 1];
      cchi[ 4] += tri[site].offd[0][ 8]  * ppsi[ 2];
      cchi[ 4] += tri[site].offd[0][ 9]  * ppsi[ 3];
      cchi[ 4] += conj(tri[site].offd[0][14]) * ppsi[ 5];

      cchi[ 5]  = tri[site].diag[0][ 5]  * ppsi[ 5];
      cchi[ 5] += tri[site].offd[0][10]  * ppsi[ 0];
      cchi[ 5] += tri[site].offd[0][11]  * ppsi[ 1];
      cchi[ 5] += tri[site].offd[0][12]  * ppsi[ 2];
      cchi[ 5] += tri[site].offd[0][13]  * ppsi[ 3];
      cchi[ 5] += tri[site].offd[0][14]  * ppsi[ 4];

      cchi[ 6]  = tri[site].diag[1][ 0]  * ppsi[ 6];
      cchi[ 6] += conj(tri[site].offd[1][ 0]) * ppsi[ 7];
      cchi[ 6] += conj(tri[site].offd[1][ 1]) * ppsi[ 8];
      cchi[ 6] += conj(tri[site].offd[1][ 3]) * ppsi[ 9];
      cchi[ 6] += conj(tri[site].offd[1][ 6]) * ppsi[10];
      cchi[ 6] += conj(tri[site].offd[1][10]) * ppsi[11];

      cchi[ 7]  = tri[site].diag[1][ 1]  * ppsi[ 7];
      cchi[ 7] += tri[site].offd[1][ 0]  * ppsi[ 6];
      cchi[ 7] += conj(tri[site].offd[1][ 2]) * ppsi[ 8];
      cchi[ 7] += conj(tri[site].offd[1][ 4]) * ppsi[ 9];
      cchi[ 7] += conj(tri[site].offd[1][ 7]) * ppsi[10];
      cchi[ 7] += conj(tri[site].offd[1][11]) * ppsi[11];

      cchi[ 8]  = tri[site].diag[1][ 2]  * ppsi[ 8];
      cchi[ 8] += tri[site].offd[1][ 1]  * ppsi[ 6];
      cchi[ 8] += tri[site].offd[1][ 2]  * ppsi[ 7];
      cchi[ 8] += conj(tri[site].offd[1][ 5]) * ppsi[ 9];
      cchi[ 8] += conj(tri[site].offd[1][ 8]) * ppsi[10];
      cchi[ 8] += conj(tri[site].offd[1][12]) * ppsi[11];

      cchi[ 9]  = tri[site].diag[1][ 3]  * ppsi[ 9];
      cchi[ 9] += tri[site].offd[1][ 3]  * ppsi[ 6];
      cchi[ 9] += tri[site].offd[1][ 4]  * ppsi[ 7];
      cchi[ 9] += tri[site].offd[1][ 5]  * ppsi[ 8];
      cchi[ 9] += conj(tri[site].offd[1][ 9]) * ppsi[10];
      cchi[ 9] += conj(tri[site].offd[1][13]) * ppsi[11];

      cchi[10]  = tri[site].diag[1][ 4]  * ppsi[10];
      cchi[10] += tri[site].offd[1][ 6]  * ppsi[ 6];
      cchi[10] += tri[site].offd[1][ 7]  * ppsi[ 7];
      cchi[10] += tri[site].offd[1][ 8]  * ppsi[ 8];
      cchi[10] += tri[site].offd[1][ 9]  * ppsi[ 9];
      cchi[10] += conj(tri[site].offd[1][14]) * ppsi[11];

      cchi[11]  = tri[site].diag[1][ 5]  * ppsi[11];
      cchi[11] += tri[site].offd[1][10]  * ppsi[ 6];
      cchi[11] += tri[site].offd[1][11]  * ppsi[ 7];
      cchi[11] += tri[site].offd[1][12]  * ppsi[ 8];
      cchi[11] += tri[site].offd[1][13]  * ppsi[ 9];
      cchi[11] += tri[site].offd[1][14]  * ppsi[10];
#elif 1
	// Unrolled version 3. 
	// Took unrolled version 2 and wrote out in real() and imag() 
	// parts. Rearranged so that all the reals follow each other
	//  in the output so that we can write linearly


	cchi[ 0].real()  = tri[site].diag[0][0].elem()  * ppsi[0].real();
	cchi[ 0].real() += tri[site].offd[0][0].real()  * ppsi[1].real();
	cchi[ 0].real() += tri[site].offd[0][0].imag()  * ppsi[1].imag();
	cchi[ 0].real() += tri[site].offd[0][1].real()  * ppsi[2].real();
	cchi[ 0].real() += tri[site].offd[0][1].imag()  * ppsi[2].imag();
	cchi[ 0].real() += tri[site].offd[0][3].real()  * ppsi[3].real();
	cchi[ 0].real() += tri[site].offd[0][3].imag()  * ppsi[3].imag();
	cchi[ 0].real() += tri[site].offd[0][6].real()  * ppsi[4].real();
	cchi[ 0].real() += tri[site].offd[0][6].imag()  * ppsi[4].imag();
	cchi[ 0].real() += tri[site].offd[0][10].real() * ppsi[5].real();
	cchi[ 0].real() += tri[site].offd[0][10].imag() * ppsi[5].imag();
	
	cchi[ 0].imag()  = tri[site].diag[0][0].elem()  * ppsi[ 0].imag();
	cchi[ 0].imag() += tri[site].offd[0][0].real()  * ppsi[1].imag();
	cchi[ 0].imag() -= tri[site].offd[0][0].imag()  * ppsi[1].real();
	cchi[ 0].imag() += tri[site].offd[0][3].real()  * ppsi[3].imag();
	cchi[ 0].imag() -= tri[site].offd[0][3].imag()  * ppsi[3].real();
	cchi[ 0].imag() += tri[site].offd[0][1].real()  * ppsi[2].imag();
	cchi[ 0].imag() -= tri[site].offd[0][1].imag()  * ppsi[2].real();
	cchi[ 0].imag() += tri[site].offd[0][6].real()  * ppsi[4].imag();
	cchi[ 0].imag() -= tri[site].offd[0][6].imag()  * ppsi[4].real();
	cchi[ 0].imag() += tri[site].offd[0][10].real() * ppsi[5].imag();
	cchi[ 0].imag() -= tri[site].offd[0][10].imag() * ppsi[5].real();
	
	
	cchi[ 1].real()  = tri[site].diag[0][ 1].elem() * ppsi[ 1].real();
	cchi[ 1].real() += tri[site].offd[0][ 0].real() * ppsi[ 0].real();
	cchi[ 1].real() -= tri[site].offd[0][ 0].imag() * ppsi[ 0].imag();
	cchi[ 1].real() += tri[site].offd[0][ 2].real() * ppsi[ 2].real();
	cchi[ 1].real() += tri[site].offd[0][ 2].imag() * ppsi[ 2].imag();
	cchi[ 1].real() += tri[site].offd[0][ 4].real() * ppsi[ 3].real();
	cchi[ 1].real() += tri[site].offd[0][ 4].imag() * ppsi[ 3].imag();
	cchi[ 1].real() += tri[site].offd[0][ 7].real() * ppsi[ 4].real();
	cchi[ 1].real() += tri[site].offd[0][ 7].imag() * ppsi[ 4].imag();
	cchi[ 1].real() += tri[site].offd[0][11].real() * ppsi[ 5].real();
	cchi[ 1].real() += tri[site].offd[0][11].imag() * ppsi[ 5].imag();
	
	
	cchi[ 1].imag()  = tri[site].diag[0][ 1].elem() * ppsi[ 1].imag();
	cchi[ 1].imag() += tri[site].offd[0][ 0].real() * ppsi[ 0].imag();
	cchi[ 1].imag() += tri[site].offd[0][ 0].imag() * ppsi[ 0].real();
	cchi[ 1].imag() += tri[site].offd[0][ 2].real() * ppsi[ 2].imag();
	cchi[ 1].imag() -= tri[site].offd[0][ 2].imag() * ppsi[ 2].real();
	cchi[ 1].imag() += tri[site].offd[0][ 4].real() * ppsi[ 3].imag();
	cchi[ 1].imag() -= tri[site].offd[0][ 4].imag() * ppsi[ 3].real();
	cchi[ 1].imag() += tri[site].offd[0][ 7].real() * ppsi[ 4].imag();
	cchi[ 1].imag() -= tri[site].offd[0][ 7].imag() * ppsi[ 4].real();
	cchi[ 1].imag() += tri[site].offd[0][11].real() * ppsi[ 5].imag();
	cchi[ 1].imag() -= tri[site].offd[0][11].imag() * ppsi[ 5].real();
	
	
	cchi[ 2].real() = tri[site].diag[0][ 2].elem()  * ppsi[ 2].real();
	cchi[ 2].real() += tri[site].offd[0][ 1].real() * ppsi[ 0].real();
	cchi[ 2].real() -= tri[site].offd[0][ 1].imag() * ppsi[ 0].imag();
	cchi[ 2].real() += tri[site].offd[0][ 2].real() * ppsi[ 1].real();
	cchi[ 2].real() -= tri[site].offd[0][ 2].imag() * ppsi[ 1].imag();
	cchi[ 2].real() += tri[site].offd[0][5].real()  * ppsi[ 3].real();
	cchi[ 2].real() += tri[site].offd[0][5].imag()  * ppsi[ 3].imag();
	cchi[ 2].real() += tri[site].offd[0][8].real()  * ppsi[ 4].real();
	cchi[ 2].real() += tri[site].offd[0][8].imag()  * ppsi[ 4].imag();
	cchi[ 2].real() += tri[site].offd[0][12].real() * ppsi[ 5].real();
	cchi[ 2].real() += tri[site].offd[0][12].imag() * ppsi[ 5].imag();
	
	
	cchi[ 2].imag() = tri[site].diag[0][ 2].elem()  * ppsi[ 2].imag();
	cchi[ 2].imag() += tri[site].offd[0][ 1].real() * ppsi[ 0].imag();
	cchi[ 2].imag() += tri[site].offd[0][ 1].imag() * ppsi[ 0].real();
	cchi[ 2].imag() += tri[site].offd[0][ 2].real() * ppsi[ 1].imag();
	cchi[ 2].imag() += tri[site].offd[0][ 2].imag() * ppsi[ 1].real();
	cchi[ 2].imag() += tri[site].offd[0][5].real()  * ppsi[ 3].imag();
	cchi[ 2].imag() -= tri[site].offd[0][5].imag()  * ppsi[ 3].real();
	cchi[ 2].imag() += tri[site].offd[0][8].real()  * ppsi[ 4].imag();
	cchi[ 2].imag() -= tri[site].offd[0][8].imag()  * ppsi[ 4].real();
	cchi[ 2].imag() += tri[site].offd[0][12].real() * ppsi[ 5].imag();
	cchi[ 2].imag() -= tri[site].offd[0][12].imag() * ppsi[ 5].real();
	
	
	cchi[ 3].real()  = tri[site].diag[0][ 3].elem() * ppsi[ 3].real();
	cchi[ 3].real() += tri[site].offd[0][ 3].real() * ppsi[ 0].real();
	cchi[ 3].real() -= tri[site].offd[0][ 3].imag() * ppsi[ 0].imag();
	cchi[ 3].real() += tri[site].offd[0][ 4].real() * ppsi[ 1].real();
	cchi[ 3].real() -= tri[site].offd[0][ 4].imag() * ppsi[ 1].imag();
	cchi[ 3].real() += tri[site].offd[0][ 5].real() * ppsi[ 2].real();
	cchi[ 3].real() -= tri[site].offd[0][ 5].imag() * ppsi[ 2].imag();
	cchi[ 3].real() += tri[site].offd[0][ 9].real() * ppsi[ 4].real();
	cchi[ 3].real() += tri[site].offd[0][ 9].imag() * ppsi[ 4].imag();
	cchi[ 3].real() += tri[site].offd[0][13].real() * ppsi[ 5].real();
	cchi[ 3].real() += tri[site].offd[0][13].imag() * ppsi[ 5].imag();
	
	
	cchi[ 3].imag()  = tri[site].diag[0][ 3].elem() * ppsi[ 3].imag();
	cchi[ 3].imag() += tri[site].offd[0][ 3].real() * ppsi[ 0].imag();
	cchi[ 3].imag() += tri[site].offd[0][ 3].imag() * ppsi[ 0].real();
	cchi[ 3].imag() += tri[site].offd[0][ 4].real() * ppsi[ 1].imag();
	cchi[ 3].imag() += tri[site].offd[0][ 4].imag() * ppsi[ 1].real();
	cchi[ 3].imag() += tri[site].offd[0][ 5].real() * ppsi[ 2].imag();
	cchi[ 3].imag() += tri[site].offd[0][ 5].imag() * ppsi[ 2].real();
	cchi[ 3].imag() += tri[site].offd[0][ 9].real() * ppsi[ 4].imag();
	cchi[ 3].imag() -= tri[site].offd[0][ 9].imag() * ppsi[ 4].real();
	cchi[ 3].imag() += tri[site].offd[0][13].real() * ppsi[ 5].imag();
	cchi[ 3].imag() -= tri[site].offd[0][13].imag() * ppsi[ 5].real();
	
	
	cchi[ 4].real()  = tri[site].diag[0][ 4].elem() * ppsi[ 4].real();
	cchi[ 4].real() += tri[site].offd[0][ 6].real() * ppsi[ 0].real();
	cchi[ 4].real() -= tri[site].offd[0][ 6].imag() * ppsi[ 0].imag();
	cchi[ 4].real() += tri[site].offd[0][ 7].real() * ppsi[ 1].real();
	cchi[ 4].real() -= tri[site].offd[0][ 7].imag() * ppsi[ 1].imag();
	cchi[ 4].real() += tri[site].offd[0][ 8].real() * ppsi[ 2].real();
	cchi[ 4].real() -= tri[site].offd[0][ 8].imag() * ppsi[ 2].imag();
	cchi[ 4].real() += tri[site].offd[0][ 9].real() * ppsi[ 3].real();
	cchi[ 4].real() -= tri[site].offd[0][ 9].imag() * ppsi[ 3].imag();
	cchi[ 4].real() += tri[site].offd[0][14].real() * ppsi[ 5].real();
	cchi[ 4].real() += tri[site].offd[0][14].imag() * ppsi[ 5].imag();
	
	
	cchi[ 4].imag()  = tri[site].diag[0][ 4].elem() * ppsi[ 4].imag();
	cchi[ 4].imag() += tri[site].offd[0][ 6].real() * ppsi[ 0].imag();
	cchi[ 4].imag() += tri[site].offd[0][ 6].imag() * ppsi[ 0].real();
	cchi[ 4].imag() += tri[site].offd[0][ 7].real() * ppsi[ 1].imag();
	cchi[ 4].imag() += tri[site].offd[0][ 7].imag() * ppsi[ 1].real();
	cchi[ 4].imag() += tri[site].offd[0][ 8].real() * ppsi[ 2].imag();
	cchi[ 4].imag() += tri[site].offd[0][ 8].imag() * ppsi[ 2].real();
	cchi[ 4].imag() += tri[site].offd[0][ 9].real() * ppsi[ 3].imag();
	cchi[ 4].imag() += tri[site].offd[0][ 9].imag() * ppsi[ 3].real();
	cchi[ 4].imag() += tri[site].offd[0][14].real() * ppsi[ 5].imag();
	cchi[ 4].imag() -= tri[site].offd[0][14].imag() * ppsi[ 5].real();
	
	
	cchi[ 5].real()  = tri[site].diag[0][ 5].elem() * ppsi[ 5].real();
	cchi[ 5].real() += tri[site].offd[0][10].real() * ppsi[ 0].real();
	cchi[ 5].real() -= tri[site].offd[0][10].imag() * ppsi[ 0].imag();
	cchi[ 5].real() += tri[site].offd[0][11].real() * ppsi[ 1].real();
	cchi[ 5].real() -= tri[site].offd[0][11].imag() * ppsi[ 1].imag();
	cchi[ 5].real() += tri[site].offd[0][12].real() * ppsi[ 2].real();
	cchi[ 5].real() -= tri[site].offd[0][12].imag() * ppsi[ 2].imag();
	cchi[ 5].real() += tri[site].offd[0][13].real() * ppsi[ 3].real();
	cchi[ 5].real() -= tri[site].offd[0][13].imag() * ppsi[ 3].imag();
	cchi[ 5].real() += tri[site].offd[0][14].real() * ppsi[ 4].real();
	cchi[ 5].real() -= tri[site].offd[0][14].imag() * ppsi[ 4].imag();
	
	
	cchi[ 5].imag()  = tri[site].diag[0][ 5].elem() * ppsi[ 5].imag();
	cchi[ 5].imag() += tri[site].offd[0][10].real() * ppsi[ 0].imag();
	cchi[ 5].imag() += tri[site].offd[0][10].imag() * ppsi[ 0].real();
	cchi[ 5].imag() += tri[site].offd[0][11].real() * ppsi[ 1].imag();
	cchi[ 5].imag() += tri[site].offd[0][11].imag() * ppsi[ 1].real();
	cchi[ 5].imag() += tri[site].offd[0][12].real() * ppsi[ 2].imag();
	cchi[ 5].imag() += tri[site].offd[0][12].imag() * ppsi[ 2].real();
	cchi[ 5].imag() += tri[site].offd[0][13].real() * ppsi[ 3].imag();
	cchi[ 5].imag() += tri[site].offd[0][13].imag() * ppsi[ 3].real();
	cchi[ 5].imag() += tri[site].offd[0][14].real() * ppsi[ 4].imag();
	cchi[ 5].imag() += tri[site].offd[0][14].imag() * ppsi[ 4].real();
	
	
	cchi[ 6].real()  = tri[site].diag[1][0].elem()  * ppsi[6].real();
	cchi[ 6].real() += tri[site].offd[1][0].real()  * ppsi[7].real();
	cchi[ 6].real() += tri[site].offd[1][0].imag()  * ppsi[7].imag();
	cchi[ 6].real() += tri[site].offd[1][1].real()  * ppsi[8].real();
	cchi[ 6].real() += tri[site].offd[1][1].imag()  * ppsi[8].imag();
	cchi[ 6].real() += tri[site].offd[1][3].real()  * ppsi[9].real();
	cchi[ 6].real() += tri[site].offd[1][3].imag()  * ppsi[9].imag();
	cchi[ 6].real() += tri[site].offd[1][6].real()  * ppsi[10].real();
	cchi[ 6].real() += tri[site].offd[1][6].imag()  * ppsi[10].imag();
	cchi[ 6].real() += tri[site].offd[1][10].real() * ppsi[11].real();
	cchi[ 6].real() += tri[site].offd[1][10].imag() * ppsi[11].imag();
	
	cchi[ 6].imag()  = tri[site].diag[1][0].elem()  * ppsi[6].imag();
	cchi[ 6].imag() += tri[site].offd[1][0].real()  * ppsi[7].imag();
	cchi[ 6].imag() -= tri[site].offd[1][0].imag()  * ppsi[7].real();
	cchi[ 6].imag() += tri[site].offd[1][1].real()  * ppsi[8].imag();
	cchi[ 6].imag() -= tri[site].offd[1][1].imag()  * ppsi[8].real();
	cchi[ 6].imag() += tri[site].offd[1][3].real()  * ppsi[9].imag();
	cchi[ 6].imag() -= tri[site].offd[1][3].imag()  * ppsi[9].real();
	cchi[ 6].imag() += tri[site].offd[1][6].real()  * ppsi[10].imag();
	cchi[ 6].imag() -= tri[site].offd[1][6].imag()  * ppsi[10].real();
	cchi[ 6].imag() += tri[site].offd[1][10].real() * ppsi[11].imag();
	cchi[ 6].imag() -= tri[site].offd[1][10].imag() * ppsi[11].real();
	
	
	cchi[ 7].real()  = tri[site].diag[1][ 1].elem()  * ppsi[ 7].real();
	cchi[ 7].real() += tri[site].offd[1][ 0].real()  * ppsi[ 6].real();
	cchi[ 7].real() -= tri[site].offd[1][ 0].imag()  * ppsi[ 6].imag();
	cchi[ 7].real() += tri[site].offd[1][ 2].real()  * ppsi[ 8].real();
	cchi[ 7].real() += tri[site].offd[1][ 2].imag()  * ppsi[ 8].imag();
	cchi[ 7].real() += tri[site].offd[1][ 4].real()  * ppsi[ 9].real();
	cchi[ 7].real() += tri[site].offd[1][ 4].imag()  * ppsi[ 9].imag();
	cchi[ 7].real() += tri[site].offd[1][ 7].real()  * ppsi[10].real();
	cchi[ 7].real() += tri[site].offd[1][ 7].imag()  * ppsi[10].imag();
	cchi[ 7].real() += tri[site].offd[1][11].real()  * ppsi[11].real();
	cchi[ 7].real() += tri[site].offd[1][11].imag()  * ppsi[11].imag();
	
	cchi[ 7].imag()  = tri[site].diag[1][ 1].elem()  * ppsi[ 7].imag();
	cchi[ 7].imag() += tri[site].offd[1][ 0].real()  * ppsi[ 6].imag();
	cchi[ 7].imag() += tri[site].offd[1][ 0].imag()  * ppsi[ 6].real();
	cchi[ 7].imag() += tri[site].offd[1][ 2].real()  * ppsi[ 8].imag();
	cchi[ 7].imag() -= tri[site].offd[1][ 2].imag()  * ppsi[ 8].real();
	cchi[ 7].imag() += tri[site].offd[1][ 4].real()  * ppsi[ 9].imag();
	cchi[ 7].imag() -= tri[site].offd[1][ 4].imag()  * ppsi[ 9].real();
	cchi[ 7].imag() += tri[site].offd[1][ 7].real()  * ppsi[10].imag();
	cchi[ 7].imag() -= tri[site].offd[1][ 7].imag()  * ppsi[10].real();
	cchi[ 7].imag() += tri[site].offd[1][11].real()  * ppsi[11].imag();
	cchi[ 7].imag() -= tri[site].offd[1][11].imag()  * ppsi[11].real();
	
	
	cchi[ 8].real()  = tri[site].diag[1][ 2].elem()  * ppsi[ 8].real();
	cchi[ 8].real() += tri[site].offd[1][ 1].real()  * ppsi[ 6].real();
	cchi[ 8].real() -= tri[site].offd[1][ 1].imag()  * ppsi[ 6].imag();
	cchi[ 8].real() += tri[site].offd[1][ 2].real()  * ppsi[ 7].real();
	cchi[ 8].real() -= tri[site].offd[1][ 2].imag()  * ppsi[ 7].imag();
	cchi[ 8].real() += tri[site].offd[1][5].real()   * ppsi[ 9].real();
	cchi[ 8].real() += tri[site].offd[1][5].imag()   * ppsi[ 9].imag();
	cchi[ 8].real() += tri[site].offd[1][8].real()   * ppsi[10].real();
	cchi[ 8].real() += tri[site].offd[1][8].imag()   * ppsi[10].imag();
	cchi[ 8].real() += tri[site].offd[1][12].real()  * ppsi[11].real();
	cchi[ 8].real() += tri[site].offd[1][12].imag()  * ppsi[11].imag();
	
	cchi[ 8].imag() = tri[site].diag[1][ 2].elem()   * ppsi[ 8].imag();
	cchi[ 8].imag() += tri[site].offd[1][ 1].real()  * ppsi[ 6].imag();
	cchi[ 8].imag() += tri[site].offd[1][ 1].imag()  * ppsi[ 6].real();
	cchi[ 8].imag() += tri[site].offd[1][ 2].real()  * ppsi[ 7].imag();
	cchi[ 8].imag() += tri[site].offd[1][ 2].imag()  * ppsi[ 7].real();
	cchi[ 8].imag() += tri[site].offd[1][5].real()   * ppsi[ 9].imag();
	cchi[ 8].imag() -= tri[site].offd[1][5].imag()   * ppsi[ 9].real();
	cchi[ 8].imag() += tri[site].offd[1][8].real()   * ppsi[10].imag();
	cchi[ 8].imag() -= tri[site].offd[1][8].imag()   * ppsi[10].real();
	cchi[ 8].imag() += tri[site].offd[1][12].real()  * ppsi[11].imag();
	cchi[ 8].imag() -= tri[site].offd[1][12].imag()  * ppsi[11].real();
	
	
	cchi[ 9].real()  = tri[site].diag[1][ 3].elem()  * ppsi[ 9].real();
	cchi[ 9].real() += tri[site].offd[1][ 3].real()  * ppsi[ 6].real();
	cchi[ 9].real() -= tri[site].offd[1][ 3].imag()  * ppsi[ 6].imag();
	cchi[ 9].real() += tri[site].offd[1][ 4].real()  * ppsi[ 7].real();
	cchi[ 9].real() -= tri[site].offd[1][ 4].imag()  * ppsi[ 7].imag();
	cchi[ 9].real() += tri[site].offd[1][ 5].real()  * ppsi[ 8].real();
	cchi[ 9].real() -= tri[site].offd[1][ 5].imag()  * ppsi[ 8].imag();
	cchi[ 9].real() += tri[site].offd[1][ 9].real()  * ppsi[10].real();
	cchi[ 9].real() += tri[site].offd[1][ 9].imag()  * ppsi[10].imag();
	cchi[ 9].real() += tri[site].offd[1][13].real()  * ppsi[11].real();
	cchi[ 9].real() += tri[site].offd[1][13].imag()  * ppsi[11].imag();
	
	cchi[ 9].imag()  = tri[site].diag[1][ 3].elem()  * ppsi[ 9].imag();
	cchi[ 9].imag() += tri[site].offd[1][ 3].real()  * ppsi[ 6].imag();
	cchi[ 9].imag() += tri[site].offd[1][ 3].imag()  * ppsi[ 6].real();
	cchi[ 9].imag() += tri[site].offd[1][ 4].real()  * ppsi[ 7].imag();
	cchi[ 9].imag() += tri[site].offd[1][ 4].imag()  * ppsi[ 7].real();
	cchi[ 9].imag() += tri[site].offd[1][ 5].real()  * ppsi[ 8].imag();
	cchi[ 9].imag() += tri[site].offd[1][ 5].imag()  * ppsi[ 8].real();
	cchi[ 9].imag() += tri[site].offd[1][ 9].real()  * ppsi[10].imag();
	cchi[ 9].imag() -= tri[site].offd[1][ 9].imag()  * ppsi[10].real();
	cchi[ 9].imag() += tri[site].offd[1][13].real()  * ppsi[11].imag();
	cchi[ 9].imag() -= tri[site].offd[1][13].imag()  * ppsi[11].real();
	
	
	cchi[10].real()  = tri[site].diag[1][ 4].elem()  * ppsi[10].real();
	cchi[10].real() += tri[site].offd[1][ 6].real()  * ppsi[ 6].real();
	cchi[10].real() -= tri[site].offd[1][ 6].imag()  * ppsi[ 6].imag();
	cchi[10].real() += tri[site].offd[1][ 7].real()  * ppsi[ 7].real();
	cchi[10].real() -= tri[site].offd[1][ 7].imag()  * ppsi[ 7].imag();
	cchi[10].real() += tri[site].offd[1][ 8].real()  * ppsi[ 8].real();
	cchi[10].real() -= tri[site].offd[1][ 8].imag()  * ppsi[ 8].imag();
	cchi[10].real() += tri[site].offd[1][ 9].real()  * ppsi[ 9].real();
	cchi[10].real() -= tri[site].offd[1][ 9].imag()  * ppsi[ 9].imag();
	cchi[10].real() += tri[site].offd[1][14].real()  * ppsi[11].real();
	cchi[10].real() += tri[site].offd[1][14].imag()  * ppsi[11].imag();
	
	cchi[10].imag()  = tri[site].diag[1][ 4].elem()  * ppsi[10].imag();
	cchi[10].imag() += tri[site].offd[1][ 6].real()  * ppsi[ 6].imag();
	cchi[10].imag() += tri[site].offd[1][ 6].imag()  * ppsi[ 6].real();
	cchi[10].imag() += tri[site].offd[1][ 7].real()  * ppsi[ 7].imag();
	cchi[10].imag() += tri[site].offd[1][ 7].imag()  * ppsi[ 7].real();
	cchi[10].imag() += tri[site].offd[1][ 8].real()  * ppsi[ 8].imag();
	cchi[10].imag() += tri[site].offd[1][ 8].imag()  * ppsi[ 8].real();
	cchi[10].imag() += tri[site].offd[1][ 9].real()  * ppsi[ 9].imag();
	cchi[10].imag() += tri[site].offd[1][ 9].imag()  * ppsi[ 9].real();
	cchi[10].imag() += tri[site].offd[1][14].real()  * ppsi[11].imag();
	cchi[10].imag() -= tri[site].offd[1][14].imag()  * ppsi[11].real();
	
	
	cchi[11].real()  = tri[site].diag[1][ 5].elem()  * ppsi[11].real();
	cchi[11].real() += tri[site].offd[1][10].real()  * ppsi[ 6].real();
	cchi[11].real() -= tri[site].offd[1][10].imag()  * ppsi[ 6].imag();
	cchi[11].real() += tri[site].offd[1][11].real()  * ppsi[ 7].real();
	cchi[11].real() -= tri[site].offd[1][11].imag()  * ppsi[ 7].imag();
	cchi[11].real() += tri[site].offd[1][12].real()  * ppsi[ 8].real();
	cchi[11].real() -= tri[site].offd[1][12].imag()  * ppsi[ 8].imag();
	cchi[11].real() += tri[site].offd[1][13].real()  * ppsi[ 9].real();
	cchi[11].real() -= tri[site].offd[1][13].imag()  * ppsi[ 9].imag();
	cchi[11].real() += tri[site].offd[1][14].real()  * ppsi[10].real();
	cchi[11].real() -= tri[site].offd[1][14].imag()  * ppsi[10].imag();
	
	cchi[11].imag()  = tri[site].diag[1][ 5].elem()  * ppsi[11].imag();
	cchi[11].imag() += tri[site].offd[1][10].real()  * ppsi[ 6].imag();
	cchi[11].imag() += tri[site].offd[1][10].imag()  * ppsi[ 6].real();
	cchi[11].imag() += tri[site].offd[1][11].real()  * ppsi[ 7].imag();
	cchi[11].imag() += tri[site].offd[1][11].imag()  * ppsi[ 7].real();
	cchi[11].imag() += tri[site].offd[1][12].real()  * ppsi[ 8].imag();
	cchi[11].imag() += tri[site].offd[1][12].imag()  * ppsi[ 8].real();
	cchi[11].imag() += tri[site].offd[1][13].real()  * ppsi[ 9].imag();
	cchi[11].imag() += tri[site].offd[1][13].imag()  * ppsi[ 9].real();
	cchi[11].imag() += tri[site].offd[1][14].real()  * ppsi[10].imag();
	cchi[11].imag() += tri[site].offd[1][14].imag()  * ppsi[10].real();
#endif
      }

      
    }// Function
  } // Namespace 

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  template<typename T, typename U>
  void QDPCloverTermT<T,U>::apply(T& chi, const T& psi, 
			    int isign, int cb) const
  {
    
    
    if ( Ns != 4 ) {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << endl;
      QDP_abort(1);
    }

    QDPCloverEnv::ApplyArgs<T> arg = { chi,psi,tri,cb };
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPCloverEnv::applySiteLoop<T>);
    // (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    
  }



  typedef QDPCloverTermT<LatticeFermion, LatticeColorMatrix> QDPCloverTerm;
  typedef QDPCloverTermT<LatticeFermionF, LatticeColorMatrixF> QDPCloverTermF;
  typedef QDPCloverTermT<LatticeFermionD, LatticeColorMatrixD> QDPCloverTermD;



#endif
