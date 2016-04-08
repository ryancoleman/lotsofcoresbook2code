#ifndef QPHIX_BLAS_MIC_NEW_H
#define QPHIX_BLAS_MIC_NEW_H

#include <qphix/geometry.h>
#include <qphix/comm.h>
#include <omp.h>

namespace QPhiX { 

  // Typically the arithmetic type is the same as the type
  template<typename FT>
  struct ArithType { 
    typedef FT Type;
  };
  
  // But sometimes it is not.
  template<>
  struct ArithType<half> { 
    typedef float Type;
  };
};

#include "qphix/site_loops.h"
#include "qphix/real_functors.h"
#include "qphix/complex_functors.h"

namespace QPhiX { 

  template<typename FT, int V, int S, bool compress>
  void copySpinor( typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res,
		   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict src,
		   const Geometry<FT,V,S,compress>& geom, 
		   int n_blas_simt) 
  {
    CopyFunctor<FT,V,S,compress> f(res, src);
    siteLoopNoReduction<FT,V,S,compress,CopyFunctor<FT,V,S,compress> >(f,geom,n_blas_simt);
  }
  
  template<typename FT, int V, int S, bool compress>
  void zeroSpinor( typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res,
		   const Geometry<FT,V,S,compress>& geom, 
		   int n_blas_simt) 
  {
    ZeroFunctor<FT,V,S,compress> f(res);
    siteLoopNoReduction<FT,V,S,compress,ZeroFunctor<FT,V,S,compress> >(f,geom,n_blas_simt);
  }
  
  template<typename FT, int V, int S, bool compress>
  void aypx(const double alpha, 
	    const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
	    typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y,
	    const Geometry<FT,V,S,compress>& geom, 
	    int n_blas_simt) 
  {
    
    AXPYFunctor<FT,V,S,compress> f(alpha, x,y );
    siteLoopNoReduction<FT,V,S,compress, AXPYFunctor<FT,V,S,compress> >(f,geom,n_blas_simt);
  }
  

  template<typename FT, int V, int S, bool compress>
  void norm2Spinor(double& n2,
		   const typename Geometry<FT,V,S,compress>::FourSpinorBlock* x,
		   Geometry<FT,V,S,compress>& geom, 
		   int n_blas_simt) 
  {
    Norm2Functor<FT,V,S,compress> f(x);
    siteLoop1Reduction<FT,V,S,compress, Norm2Functor<FT,V,S,compress> >(f,n2,geom,n_blas_simt);
  }



  template<typename FT, int V, int S, bool compress>
  void xmyNorm2Spinor(typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict res,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
		      typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y,
		      double& n2res,
		      const Geometry<FT,V,S,compress>& geom, 
		      int n_blas_simt) 
  {
    XMYNorm2Functor<FT,V,S,compress> f(res,x,y);
    siteLoop1Reduction<FT,V,S,compress, XMYNorm2Functor<FT,V,S,compress> >(f,n2res,geom,n_blas_simt);
  } // End of Function. 
  
  
  
  template<typename FT, int V, int S, bool compress>
  void rmammpNorm2rxpap(
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r,
			const double& ar, 
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict mmp,
			double& cp, 
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p,
			const Geometry<FT,V,S,compress>& geom, 
			int n_blas_simt) 
  {
    RmammpNorm2rxpapFunctor<FT,V,S,compress> f(ar,r,mmp,x,p);
    siteLoop1Reduction<FT,V,S, compress, RmammpNorm2rxpapFunctor<FT,V,S,compress> >(f,cp,geom,n_blas_simt);
  } // End of Function. 
 
 template<typename FT, int V, int S, bool compress>
  void richardson_rxupdateNormR(
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
			typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r,
			const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict delta_x,
			const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict delta_r,
			double& cp, 
			const Geometry<FT,V,S,compress>& geom, 
			int n_blas_simt) 
  {
    RichardsonRXUpdateNormRFunctor<FT,V,S,compress> f(x,r,delta_x,delta_r);
    siteLoop1Reduction<FT,V,S,compress, RichardsonRXUpdateNormRFunctor<FT,V,S,compress> >(f,cp,geom,n_blas_simt);
  } // End of Function. 
 

  template<typename FT, int V, int S, bool compress>
    void bicgstab_xmy( const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
	       	       typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y,
		       const Geometry<FT,V,S,compress>& geom, 
		       int n_blas_simt) 
  {
    XMYFunctor<FT,V,S,compress> f(x,y);
    siteLoopNoReduction<FT,V,S,compress, XMYFunctor<FT,V,S,compress> >(f,geom,n_blas_simt);
  } // End of Function. 
  

  template<typename FT, int V, int S, bool compress> 
    void innerProduct(double results[2],
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict y,
		      const Geometry<FT,V,S,compress>& geom, 
		      int n_blas_simt) 
  {
    InnerProductFunctor<FT,V,S,compress> f(x,y);
    siteLoop2Reductions<FT,V,S,compress, InnerProductFunctor<FT, V, S,compress> >(f, results, geom, n_blas_simt);
  }

  template<typename FT, int V, int S, bool compress>
    void 
    bicgstab_p_update(		     
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r,
		      typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict v,
		      double beta[2],
		      double omega[2],
		      const Geometry<FT,V,S,compress>& geom, 
		      int n_blas_simt) {

    BiCGStabPUpdateFunctor<FT,V,S,compress> f(r, p, v, beta, omega);
    siteLoopNoReduction<FT,V,S,compress, BiCGStabPUpdateFunctor<FT,V,S,compress>  >(f,geom,n_blas_simt);
  }

  template<typename FT, int V, int S, bool compress>
    void 
    bicgstab_s_update(
		      double alpha[2],		     
		      typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict s,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict v,
		      
		      const Geometry<FT,V,S,compress>& geom, 
		      int n_blas_simt) {
    
    BiCGStabSUpdateFunctor<FT,V,S,compress> f(alpha,s,v);
    siteLoopNoReduction<FT,V,S,compress, BiCGStabSUpdateFunctor<FT,V,S,compress>  >(f,geom,n_blas_simt);
  }

  template<typename FT, int V, int S, bool compress>
    void 
    bicgstab_rxupdate(
		      typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict x,
		      typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict r,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict t,
		      const typename Geometry<FT,V,S,compress>::FourSpinorBlock* restrict p,
		      double omega[2],		     
		      double alpha[2],
		      double& r_norm,
		      const Geometry<FT,V,S,compress>& geom, 
		      int n_blas_simt) {
    
    BiCGStabRXUpdateFunctor<FT,V,S,compress> f(x,r,t,p,omega,alpha);
    siteLoop1Reduction<FT,V,S,compress, BiCGStabRXUpdateFunctor<FT,V,S,compress>  >(f,r_norm,geom,n_blas_simt);
  }
		      

  template<typename FTOut, int VOut, int SOut, bool CompressOut,
	   typename FTIn, int VIn, int SIn, bool CompressIn>
  void convert( 
	       typename Geometry<FTOut,VOut,SOut,CompressOut>::FourSpinorBlock* restrict spinor_out,
	       double scale_factor,
	       const typename Geometry<FTIn,VIn,SIn,CompressIn>::FourSpinorBlock* restrict spinor_in,
	       const Geometry<FTOut,VOut,SOut,CompressOut>& geom_out,
	       const Geometry<FTIn,VIn,SIn,CompressIn>& geom_in,
	       int n_blas_threads)
  {
    // Get the subgrid latt size.
    int Nt =  geom_out.Nt();
    int Nz =  geom_out.Nz();
    int Ny =  geom_out.Ny();
    int Nxh = geom_out.Nxh();
    int nvecs_out = geom_out.nVecs();
    int Pxy_out = geom_out.getPxy();
    int Pxyz_out = geom_out.getPxyz();

    int nvecs_in = geom_in.nVecs();
    int Pxy_in = geom_in.getPxy();
    int Pxyz_in = geom_in.getPxyz();

#pragma omp parallel for collapse(4)
    for(int t=0; t < Nt; t++) {
      for(int z=0; z < Nz; z++) {
	for(int y=0; y < Ny; y++) {
	  for(int s=0; s < nvecs_out; s++) { 
	    for(int col=0; col < 3; col++)  {
	      for(int spin=0; spin < 4; spin++) { 
		for(int x=0; x < SOut; x++) { 
		  
		  int ind_out = t*Pxyz_out+z*Pxy_out+y*nvecs_out+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		  int x_coord = s*SOut + x;
		  
		  int s_in  = x_coord / SIn;
		  int x_in  = x_coord - SIn*s_in;
		  
		  int ind_in = t*Pxyz_in + z*Pxy_in + y*nvecs_in + s_in;
		    
		  spinor_out[ind_out][col][spin][0][x] = 
		    rep<FTOut, typename ArithType<FTOut>::Type >(
			   rep<typename ArithType<FTOut>::Type,double>(scale_factor)
			   * rep<typename ArithType<FTOut>::Type, FTIn >( spinor_in[ind_in][col][spin][0][x_in] )
								 );

		  spinor_out[ind_out][col][spin][1][x] = 
		    rep<FTOut, typename ArithType<FTOut>::Type >(
			   rep<typename ArithType<FTOut>::Type,double>(scale_factor)
			   * rep<typename ArithType<FTOut>::Type, FTIn >( spinor_in[ind_in][col][spin][1][x_in] )
								 );



		  
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
    



}; // Namespace

#endif
