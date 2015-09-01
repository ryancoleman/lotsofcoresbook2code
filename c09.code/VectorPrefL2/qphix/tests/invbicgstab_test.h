#ifndef BICGSTAB_DEBUG_H
#define BICGSTAB_DEBUG_H


#undef DEBUG

template<typename T, typename U>
  void AMat( T& res,
	     const multi1d<U>& u,
	     const T& src, 
	     int isign,
	     int targetCB,
	     const Real& massFactor, 
	     const Real& betaFactor) 
{
  int otherCB = (targetCB == 0) ? 1 : 0;
  T tmp1;
  T tmp2;

  dslash(tmp1,u, src,isign,otherCB);
  dslash(tmp2,u, tmp1,isign,targetCB);
  res[rb[targetCB]] = massFactor*src - betaFactor*tmp2;
}

template<typename T, typename U, typename CR>
void
InvBiCGStabSolve(const Real& massFactor,
	    const Real& betaFactor,
	    multi1d<U> u,
	    const T& chi,
	    T& psi,
	    const Real& RsdBiCGStab,
	    int MaxBiCGStab, 
	    int isign,
	    int targetCB)

{
  int niters =0 ; 
  bool convP = false;	
  Double chi_sq =  norm2(chi,rb[targetCB]);
  Double rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;

  QDPIO::cout << " || chi ||^2 || =  " << chi_sq <<endl;
  QDPIO::cout << " Target is = " << rsd_sq << endl;

  // First get r = r0 = chi - A psi
  T r;
  T r0;

  // Get A psi, use r0 as a temporary
  //A(r0, psi, isign);
  AMat(r0,u,psi,isign,targetCB,massFactor,betaFactor);

  // now work out r= chi - Apsi = chi - r0
  r[rb[targetCB]] = chi - r0;
  // Also copy back to r0. We are no longer in need of the
  // nth component
  r0[rb[targetCB]] = r;
  
  // Now we have r = r0 = chi - Mpsi
 
  // Now initialise v = p = 0
  T p;
  T v;

  p[rb[targetCB]] = zero;
  v[rb[targetCB]] = zero;

  T tmp;
  T t;

  ComplexD rho, rho_prev, alpha, omega;

  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho_prev = Double(1);
  alpha = Double(1);
  omega = Double(1);

  // The iterations 
  for(int k = 1; k <= MaxBiCGStab && !convP ; k++) { 
    
    // rho_{k+1} = < r_0 | r >
    rho = innerProduct(r0,r,rb[targetCB]);
#ifdef DEBUG
    QDPIO::cout << "k=" << k << " rho = " << rho << endl;
#endif

    if( toBool( real(rho) == 0 ) && toBool( imag(rho) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: rho = 0" << endl;
      QDP_abort(1);
    }

    // beta = ( rho_{k+1}/rho_{k})(alpha/omega)
    ComplexD beta;
    beta = ( rho / rho_prev ) * (alpha/omega);
    
    // p = r + beta(p - omega v)

    // first work out p - omega v 
    // into tmp
    // then do p = r + beta tmp
    CR omega_r = omega;
    CR beta_r = beta;
#ifdef DEBUG
    QDPIO::cout << "k=" << k << " beta=" << beta << "  omega=" << omega << endl;
#endif

 #ifdef DEBUG 
    QDPIO::cout << "k=" << k << " p_norm before p_update = " << norm2(p, rb[targetCB]) << endl;
#endif 

    tmp[rb[targetCB]] = p - omega_r*v;
    p[rb[targetCB]] = r + beta_r*tmp;

 #ifdef DEBUG 
    QDPIO::cout << "k=" << k << " p_norm after p_update =" << norm2(p, rb[targetCB]) << endl;
#endif 

    // v = Ap
    //    A(v,p,isign);
    AMat(v,u,p,isign,targetCB,massFactor,betaFactor);
    
    // alpha = rho_{k+1} / < r_0 | v >
    // put <r_0 | v > into tmp
#ifdef DEBUG

    QDPIO::cout << "k=" << k << "r0_norm=" << norm2(r0,rb[targetCB]) <<  "  vnorm="<< norm2(v,rb[targetCB]) << endl;
#endif

    DComplex ctmp = innerProduct(r0,v,rb[targetCB]);
#ifdef DEBUG
    QDPIO::cout << "k=" << k << "<r0|v>=" << ctmp << endl;
#endif

    if( toBool( real(ctmp) == 0 ) && toBool( imag(ctmp) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }

    alpha = rho / ctmp;
#ifdef DEBUG
    QDPIO::cout << "k=" << "  rho=" << rho << " ctmp=" << ctmp << " alpha=" << alpha << endl;
#endif
 
    // Done with rho now, so save it into rho_prev
    rho_prev = rho;

    // s = r - alpha v
    // I can overlap s with r, because I recompute it at the end.
    CR alpha_r = alpha;
    QDPIO::cout << "k=" << k << " alpha=" << alpha_r << endl; 
   r[rb[targetCB]]  -=  alpha_r*v;

#ifdef DEBUG
    QDPIO::cout << "k=" << k << " r_norm before tmult=" << norm2(r,rb[targetCB]) << endl;
#endif
    // t = As  = Ar 
    //    A(t,r,isign);
    AMat(t,u,r,isign,targetCB,massFactor,betaFactor);
    // omega = < t | s > / < t | t > = < t | r > / norm2(t);

    // This does the full 5D norm
    Double t_norm = norm2(t,rb[targetCB]);
#ifdef DEBUG
    QDPIO::cout << "k=" << k << " t_norm=" << t_norm << endl;
#endif

    if( toBool(t_norm == 0) ) { 
      QDPIO::cerr << "Breakdown || Ms || = || t || = 0 " << endl;
      QDP_abort(1);
    }

    // accumulate <t | s > = <t | r> into omega
    omega = innerProduct(t,r,rb[targetCB]);
    omega /= t_norm;

    // psi = psi + omega s + alpha p 
    //     = psi + omega r + alpha p
    //
    // use tmp to compute psi + omega r
    // then add in the alpha p
#ifdef DEBUG
    QDPIO::cout << "k=" << k << " alpha=" << alpha << " omega=" << omega << endl;
#endif

    omega_r = omega;
    alpha_r = alpha;
    tmp[rb[targetCB]] = psi + omega_r*r;   
    psi[rb[targetCB]] = tmp + alpha_r*p;



    // r = s - omega t = r - omega t1G

    
    r[rb[targetCB]] -= omega_r*t;


    Double r_norm = norm2(r,rb[targetCB]);


    QDPIO::cout << "Iteration " << k << " : r = " << r_norm << endl;

    if( toBool(r_norm < rsd_sq ) ) {
      convP = true;
      niters = k;
    }
    else { 
      convP = false;
      niters = k;
    }

    //-------BiCGStab Flopcounting --------------------------------------
    // flopcount.addSiteFlops(8*Nc*Ns,s);     // <r0|r>
    // flopcount.addSiteFlops(16*Nc*Ns,s);    // p = r + beta p - beta_omega v
    // flopcount.addSiteFlops(8*Nc*Ns,s);  //  <r0 | v>
    // flopcount.addSiteFlops(8*Nc*Ns,s);  //  r -= alpha v
    // flopcount.addSiteFlops(8*Nc*Ns, s); //  < t, r>
    // flopcount.addSiteFlops(4*Nc*Ns, s); //  < t, t>     
    // flopcount.addSiteFlops(16*Nc*Ns,s); // psi += omega r + alpha_p
    // flopcount.addSiteFlops(8*Nc*Ns,s); // r -=omega t
    // flopcount.addSiteFlops(4*Nc*Ns,s); // norm2(r)
    // flopcount.addFlops(2*A.nFlops());  // = 80*Nc*Ns cbsite flops + 2*A
    //----------------------------------------------------------------------
  }
  

  QDPIO::cout << "Finished in " << niters << " iterations" << endl;
  QDPIO::cout << "ConvP is " << convP << endl;
}

#endif
