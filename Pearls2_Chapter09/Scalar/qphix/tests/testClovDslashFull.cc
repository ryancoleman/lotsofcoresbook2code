#include "unittest.h"
#include "testClovDslashFull.h"

// Stupid compiler thing
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END

#include "qphix/clover_dslash_def.h"
#include "qphix/clover_dslash_body.h"


#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif


#include "qphix/qdp_packer.h"
#include "clover_term_qdp_w.h"

#if 1
#include "qphix/clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#endif

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 4
#endif

#if defined(QPHIX_MIC_SOURCE)

#define VECLEN_SP 16 
#define VECLEN_HP 16 
#define VECLEN_DP 8

#elif defined(QPHIX_AVX_SOURCE) 

#define VECLEN_SP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SCALAR_SOURCE)
#define VECLEN_DP 1
#define VECLEN_SP 1
#define QPHIX_SOALEN 1

#elif defined(QPHIX_QPX_SOURCE)
#define VECLEN_DP 4
#define QPHIX_SOALEN 4

#endif

template<typename T>
struct tolerance { 
  static const Double small; // Always fail
};

template<>
const Double tolerance<half>::small = Double(5.0e-3);

template<>
 const Double tolerance<float>::small = Double(1.0e-6);


template<>
const Double tolerance<double>::small = Double(1.0e-13);

template<typename T>
struct rsdTarget { 
  static const double value;
};

template<>
const double rsdTarget<half>::value = (double)(1.0e-3);

template<>
const double rsdTarget<float>::value = (double)(1.0e-7);


template<>
const double rsdTarget<double>::value = (double)(1.0e-12);



template<typename FT, int V, int S, bool compress, typename U, typename Phi>
void
testClovDslashFull::runTest(void) 
{

  typedef typename Geometry<FT,V,S,compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT,V,S,compress>::CloverBlock Clover;

  bool verbose = true;
  QDPIO::cout << endl << "ENTERING CLOVER DSLASH TEST" << endl;

  // Diagnostic information:
  const multi1d<int>& lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Ny = lattSize[1];
  int Nz = lattSize[2];
  int Nt = lattSize[3];



  QDPIO::cout << "VECLEN=" << V<< "  SOALEN="<< S <<endl;
  QDPIO::cout << "Lattice Size: ";
  for(int mu=0; mu < lattSize.size(); mu++){ 
    QDPIO::cout << " " << lattSize[mu];
  }
  QDPIO::cout << endl;

  QDPIO::cout << "Block Sizes: By=" << By << " Bz=" << Bz << endl;
  QDPIO::cout << "N Cores" << NCores << endl;
  QDPIO::cout << "SMT Grid: Sy=" << Sy << " Sz=" << Sz << endl;
  QDPIO::cout << "Pad Factors: PadXY=" << PadXY << " PadXYZ=" << PadXYZ << endl;
  QDPIO::cout << "MinCt=" << MinCt << endl;
  QDPIO::cout << "Threads_per_core = " << N_simt << endl;

 

  // What we consider to be small enough...
  QDPIO::cout << "Inititalizing QDP++ gauge field"  << endl;


  // Make a random gauge field 
  // Start up the gauge field somehow
  // We choose: u = reunit(  1 + factor*gaussian(u) );
  // Adjust factor to get reasonable inversion times in invert test.
  // Bug gauge field is always unitary

  multi1d<U> u(4);
  {
    U g;
    U uf;
    for(int mu=0; mu < 4; mu++) { 
#if 1
      uf = 1;   // Unit gauge
      
      Real factor=Real(0.08);
      gaussian(g);
      u[mu] = uf + factor*g;
      reunit(u[mu]);
#else
      // Only Unit Gauge: testing....
      u[mu] = 1;
#endif
    }
  }


  // Set anisotropy parameters -- pick some random numbers
  double xi_0_f = 0.3;
  double nu_f = 1.4;

  // This is what makeFermCoeffs does under the hood.
  // Get the direct spatial and temporal anisotropy factors
  double aniso_fac_s=(double)nu_f / xi_0_f;
  double aniso_fac_t=(double)(1);

  QDPIO::cout << "Setting Clover Term Parameters" << endl;

  CloverFermActParams clparam;
  AnisoParam_t aniso;

  // Aniso prarams
  aniso.anisoP=true;
  aniso.xi_0 = xi_0_f;
  aniso.nu = nu_f;
  aniso.t_dir = 3;

  // Set up the Clover params
  clparam.anisoParam = aniso;

  // Some mass
  clparam.Mass = Real(0.01);
  
  // Some random clover coeffs
  clparam.clovCoeffR=Real(1);
  clparam.clovCoeffT=Real(1);

  // Set up the 'periodic BC dslash'
  QDPIO::cout << "Dslash will run with " << omp_get_max_threads() << " threads" << endl;

  double t_boundary = (double)(1);
  QDPIO::cout << "Instantiating ClovDslash<FT,"<<V<<","<<S<<">" << " with t_boundary = " << t_boundary <<  endl;

 
  Geometry<FT,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  ClovDslash<FT,V,S,compress> D32(&geom, t_boundary, aniso_fac_s, aniso_fac_t);

  // Make a random source
  QDPIO::cout << "Initializing QDP++ input spinor" << endl;
  Phi chi, clov_chi;
  Phi psi, chi2, clov_chi2;
  QDPIO::cout << "Filling psi with gaussian noise" << endl;
  gaussian(psi);
  
  
  QDPIO::cout << "Allocating packged gauge fields" << endl;
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();

  QDPIO::cout << "Allocating packed spinor fields" << endl;
  Spinor* psi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* psi_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_odd=(Spinor*)geom.allocCBFourSpinor();

  QDPIO::cout << "Allocate Packed Clover Term" << endl;
  Clover* A_cb0=(Clover*)geom.allocCBClov();
  Clover* A_cb1=(Clover*)geom.allocCBClov();
  Clover* A_inv_cb0=(Clover*)geom.allocCBClov();
  Clover* A_inv_cb1=(Clover*)geom.allocCBClov();
  Clover* invclov_packed[2] = { A_inv_cb0, A_inv_cb1 };
  Clover* clov_packed[2] = { A_cb0, A_cb1 };

  QDPIO::cout << "Fields allocated" << endl;
 
  
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  QDPIO::cout << "done" << endl;
  
  QDPIO::cout << " Packing fermions..." ;	
  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };

  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
  QDPIO::cout << "done" << endl; 

  QDPIO::cout << "Creating the Clover Term " << endl;

  // Clover term deals with anisotropy internally -- so use original u field.
  QDPCloverTermT<Phi, U> clov_qdp;
  clov_qdp.create(u, clparam);
  QDPIO::cout << "Inverting Clover Term" << endl;
  QDPCloverTermT<Phi, U> invclov_qdp(clov_qdp);
  for(int cb=0; cb < 2; cb++) { 
    invclov_qdp.choles(cb);
  }
  QDPIO::cout << "Done" << endl;

  QDPIO::cout << "Packing Clover term..." << endl;
  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(invclov_qdp.getTriBuffer(), invclov_packed[cb], geom, cb);
  }

  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(clov_qdp.getTriBuffer(), clov_packed[cb], geom, cb);
  }
  QDPIO::cout << "Done" << endl;

  QDPIO::cout << "Multiplying aniso coeffs into gauge field for testing" << endl;
  multi1d<U> u_test(4);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor;
    if (mu == 3) { 
      factor=Real(aniso_fac_t);
    }
    else { 
      factor=Real(aniso_fac_s);
    }
    u_test[mu] = factor*u[mu];
  }

#if 1
  // Test only Dslash operator.
  // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
  QDPIO::cout << "Testing Dslash \n" << endl;

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 

  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;

      clov_chi = zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);


      // Apply Optimized Dslash
      D32.dslash(chi_s[target_cb],	
		 psi_s[source_cb],
		 u_packed[target_cb],
		 invclov_packed[target_cb],
		 isign, 
		 target_cb);

      qdp_unpack_spinor<>(chi_even,chi_odd, clov_chi, geom);

      // Account for Clover term from QDP++
      // We should remove this once the actual clover implementation is ready.
      // invclov_qdp.apply(clov_chi,chi,isign, target_cb);
      // clov_chi = chi;

      
      // Apply QDP Dslash
      chi2 = zero;

      dslash(chi2,u_test,psi, isign, target_cb);
      invclov_qdp.apply(clov_chi2,chi2,isign, target_cb);

      // Check the difference per number in chi vector
      Phi diff = clov_chi2 -clov_chi;
  
      Double diff_norm = sqrt( norm2( diff, rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm > tolerance<FT>::small  ) ) {
	
	int Nxh=Nx/2;
	for(int t=0; t < Nt; t++){ 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++){ 
	      for(int x=0; x < Nxh; x++){ 
		
		// These are unpadded QDP++ indices...
		int ind = x + Nxh*(y + Ny*(z + Nz*t));
		for(int s =0 ; s < Ns; s++) { 
		  for(int c=0; c < Nc; c++) { 
		    REAL dr = diff.elem(rb[target_cb].start()+ind).elem(s).elem(c).real();
		    REAL di = diff.elem(rb[target_cb].start()+ind).elem(s).elem(c).imag();
		    if( toBool( fabs(dr) > tolerance<FT>::small ) || toBool ( fabs(di) > tolerance<FT>::small) ) {
		      QDPIO::cout <<"(x,y,z,t)=(" << x <<"," <<y<<","<<z<<","<<t<<") site=" << ind << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+ind).elem(s).elem(c) 
				  << "  chi = " << clov_chi.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << " qdp++ =" << clov_chi2.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << endl;
		      
		    }
		  }
		}
	      } // x 
	    } // y 
	  } // z 
	} // t
	assertion( toBool( diff_norm < tolerance<FT>::small ) );
      }


    } // cb
  } // isign


#endif



#if 1
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  // Test ax - bDslash y
  QDPIO::cout << "Testing dslashAchiMinusBDPsi" << endl;

  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;




      chi=zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      double beta = (double)(0.25); // Always 0.25

      // Apply Optimized Dslash
      D32.dslashAChiMinusBDPsi(chi_s[target_cb],	
			       psi_s[source_cb],
			       psi_s[target_cb],
			       u_packed[target_cb],
			       clov_packed[target_cb],
			       beta,
			       isign, 
			       target_cb);
      
      qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);


      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, target_cb);
      Phi res = zero;
      clov_qdp.apply(res,psi,isign,target_cb);
      res[rb[target_cb]] -= beta*chi2;

      // Check the difference per number in chi vector
      Phi diff = res-chi;

      Double diff_norm = sqrt( norm2( diff , rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm > tolerance<FT>::small ) ) {
	for(int i=0; i < rb[target_cb].siteTable().size(); i++){ 
	  for(int s =0 ; s < Ns; s++) { 
	    for(int c=0; c < Nc; c++) { 
	      QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+i).elem(s).elem(c) << endl;
	    }
	  }
	}
      }
      assertion( toBool( diff_norm < tolerance<FT>::small ) );

    }
   }
    
#endif


#if 1
  // Test only Dslash operator.
  // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
  QDPIO::cout << "Testing Dslash With antiperiodic BCs \n" << endl;
  t_boundary = (double)(-1);
  
  // Create Antiperiodic Dslash
  ClovDslash<FT,V,S,compress> D32_ap(&geom, t_boundary, aniso_fac_s, aniso_fac_t);
 
  
  // Step 1: Convert u_test into one with antiperiodic BCs.
  // NB: This alone does not need a re-pack for the gauge field. 
  QDPIO::cout << "Applying BCs to original gauge field" << endl;
  int mu_t = Nd -1;

  // Let us make u_test to be u with antiperiodic_BCs
  for(int mu=0; mu < Nd; mu++) { 
    u_test[mu] = u[mu];
  }

  u_test[mu_t] *= where(Layout::latticeCoordinate(mu_t) == (Layout::lattSize()[mu_t]-1),
  			Real(t_boundary), Real(1));

  QDPIO::cout << "Creating Clover term using Gauge field with antiperiodic BCs" << endl;
  QDPCloverTermT<Phi,U> clov_qdp_ap;
  clov_qdp_ap.create(u_test, clparam);
  QDPIO::cout << "Inverting Clover Term" << endl;
  QDPCloverTermT<Phi,U> invclov_qdp_ap(clov_qdp_ap);
  for(int cb=0; cb < 2; cb++) { 
    invclov_qdp_ap.choles(cb);
  }
  QDPIO::cout << "Done" << endl;

  // Now we need to repack this.
  QDPIO::cout << "Packing Clover term..." << endl;
  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(invclov_qdp_ap.getTriBuffer(), invclov_packed[cb], D32_ap.getGeometry(), cb);
  }

  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(clov_qdp_ap.getTriBuffer(), clov_packed[cb], D32_ap.getGeometry(), cb);
  }

  QDPIO::cout << "Folding aniso factors into gauge field for testing" << endl;
  for(int mu=0; mu < Nd; mu++) { 
    Real factor;
    if (mu == 3) { 
      factor=Real(aniso_fac_t);
    }
    else { 
      factor=Real(aniso_fac_s);
    }
    u_test[mu] *= factor;
  }

  // NB: Gauge field doesn't need to be repacked. It is the original 'u' without the aniso factors or boundaries
  // As these are now handled in the Dslash.
  
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;

      clov_chi = zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, D32_ap.getGeometry());


      // Apply Optimized Dslash
      D32_ap.dslash(chi_s[target_cb],	
		    psi_s[source_cb],
		    u_packed[target_cb],
		    invclov_packed[target_cb],
		    isign,
		    target_cb);

      qdp_unpack_spinor<>(chi_even,chi_odd, clov_chi, geom);

      // Account for Clover term from QDP++
      // We should remove this once the actual clover implementation is ready.
      // invclov_qdp.apply(clov_chi,chi,isign, target_cb);
      // clov_chi = chi;

      
      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, target_cb);
      invclov_qdp_ap.apply(clov_chi2,chi2,isign, target_cb);

      // Check the difference per number in chi vector
      Phi diff = clov_chi2 -clov_chi;
  
      Double diff_norm = sqrt( norm2( diff, rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm > tolerance<FT>::small ) ) {
	
	int Nxh=Nx/2;
	for(int t=0; t < Nt; t++){ 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++){ 
	      for(int x=0; x < Nxh; x++){ 
		
		// These are unpadded QDP++ indices...
		int ind = x + Nxh*(y + Ny*(z + Nz*t));
		for(int s =0 ; s < Ns; s++) { 
		  for(int c=0; c < Nc; c++) { 
		    REAL dr = diff.elem(rb[target_cb].start()+ind).elem(s).elem(c).real();
		    REAL di = diff.elem(rb[target_cb].start()+ind).elem(s).elem(c).imag();
		    if( toBool( fabs(dr) > tolerance<FT>::small ) || toBool ( fabs(di) > tolerance<FT>::small ) ) {
		      QDPIO::cout <<"(x,y,z,t)=(" << x <<"," <<y<<","<<z<<","<<t<<") site=" << ind << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+ind).elem(s).elem(c) 
				  << "  chi = " << clov_chi.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << " qdp++ =" << clov_chi2.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << endl;
		      
		    }
		  }
		}
	      } // x 
	    } // y 
	  } // z 
	} // t
	assertion( toBool( diff_norm < tolerance<FT>::small ) );
      }


    } // cb
  } // isign


#endif



#if 1
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  // Test ax - bDslash y
  QDPIO::cout << "Testing dslashAchiMinusBDPsi" << endl;

  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;




      chi=zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, D32_ap.getGeometry());

      double beta = (double)(0.25); // Always 0.25

      // Apply Optimized Dslash
      D32_ap.dslashAChiMinusBDPsi(chi_s[target_cb],	
				  psi_s[source_cb],
				  psi_s[target_cb],
				  u_packed[target_cb],
				  clov_packed[target_cb],
				  beta,
				  isign, 
				  target_cb);
      
      qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, D32_ap.getGeometry());


      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, target_cb);
      Phi res = zero;
      clov_qdp_ap.apply(res,psi,isign,target_cb);
      res[rb[target_cb]] -= beta*chi2;

      // Check the difference per number in chi vector
      Phi diff = res-chi;

      Double diff_norm = sqrt( norm2( diff , rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm > tolerance<FT>::small ) ) {
	for(int i=0; i < rb[target_cb].siteTable().size(); i++){ 
	  for(int s =0 ; s < Ns; s++) { 
	    for(int c=0; c < Nc; c++) { 
	      QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+i).elem(s).elem(c) << endl;
	    }
	  }
	}
      }
      assertion( toBool( diff_norm < tolerance<FT>::small ) );

    }
   }
    
#endif

  // Disabling testing the even odd operator until recoded with new vectorization
#if 1
  QDPIO::cout << "Testing Even Odd Operator" << endl;
  t_boundary=(double)(-1);
  EvenOddCloverOperator<FT,V,S,compress> M(u_packed,  
					   clov_packed[1], 
					   invclov_packed[0],  
					   &geom,
					   t_boundary,
					   aniso_fac_s,
					   aniso_fac_t);
  Phi ltmp=zero;
  Real betaFactor=Real(0.25);
   // Apply optimized
  for(int isign=1; isign >= -1; isign -=2) {
    
      chi=zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      M( chi_s[1],
	 psi_s[1],
	 isign);

      qdp_unpack_spinor<> (chi_s[0], chi_s[1], chi, geom);


      // Apply QDP Dslash
      chi2 = zero;
      
      dslash(chi2,u_test,psi, isign, 0);
      invclov_qdp_ap.apply(clov_chi2, chi2, isign, 0);
      dslash(ltmp,u_test,clov_chi2, isign, 1);
      
      clov_qdp_ap.apply(chi2, psi, isign, 1);

      chi2[rb[1]] -= betaFactor*ltmp;

      // Check the difference per number in chi vector

      Phi diff = zero;
      diff[rb[1]]=chi2-chi;

      Double diff_norm = sqrt( norm2( diff , rb[1] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t isign = " << isign << "  diff_norm = " << diff_norm << endl;      
     // Assert things are OK...
      if ( toBool( diff_norm > tolerance<FT>::small ) ) {
	for(int i=0; i < rb[1].siteTable().size(); i++){ 
	  for(int s =0 ; s < Ns; s++) { 
	    for(int c=0; c < Nc; c++) { 
	      REAL re=  diff.elem(rb[1].start()+i).elem(s).elem(c).real();
	      REAL im=  diff.elem(rb[1].start()+i).elem(s).elem(c).imag();
	      if( toBool ( fabs(re) > tolerance<FT>::small ) || toBool( fabs(im) > tolerance<FT>::small ) )  { 
		QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[1].start()+i).elem(s).elem(c) << endl;
	      }

	    }
	  }
	}
      }
      assertion( toBool( diff_norm < tolerance<FT>::small ) );

    
  }

#endif

#if 1
  {
    chi = zero;
    qdp_pack_cb_spinor<>(chi, chi_s[1], geom,1);

    double rsd_target=rsdTarget<FT>::value;
    int max_iters=500;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;
    
    InvCG<FT,V,S,compress> solver(M, max_iters);
    solver.tune();

    double start = omp_get_wtime();
    solver(chi_s[1], psi_s[1],rsd_target, niters, rsd_final, site_flops, mv_apps,1,verbose);
    double end = omp_get_wtime();
    
    
    
    qdp_unpack_cb_spinor<>(chi_s[1], chi, geom,1);
    
    // Multiply back 
    // chi2 = M chi
    dslash(chi2,u_test,chi, 1, 0);
    invclov_qdp_ap.apply(clov_chi2, chi2, 1, 0);
    dslash(ltmp,u_test,clov_chi2, 1, 1);
    
    clov_qdp_ap.apply(chi2, chi, 1, 1);
    chi2[rb[1]] -= betaFactor*ltmp;
    
    
    // chi3 = M^\dagger chi2
    Phi chi3 = zero;
    dslash(chi3,u_test,chi2, -1, 0);
    invclov_qdp_ap.apply(clov_chi2, chi3, -1, 0);
    dslash(ltmp,u_test,clov_chi2, -1, 1);
    
    clov_qdp_ap.apply(chi3, chi2, -1, 1);
    chi3[rb[1]] -= betaFactor*ltmp;
    
    //  dslash(chi3,u,chi2, (-1), 1);
    // dslash(ltmp,u,chi3, (-1), 0);
    // chi3[rb[0]] = massFactor*chi2 - betaFactor*ltmp;
    
    Phi diff = chi3 - psi;
    QDPIO::cout << "True norm is: " << sqrt(norm2(diff, rb[1])/norm2(psi,rb[1])) << endl;
    
    int Nxh = Nx/2;
    unsigned long num_cb_sites=Nxh * Ny * Nz * Nt;
    
    unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;

    masterPrintf("GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));
  }
#endif

#if 1
  {
    chi = zero;
    qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
    
    double rsd_target=rsdTarget<FT>::value;
    int max_iters=500;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;
    
    InvBiCGStab<FT,V,S,compress> solver(M, max_iters);
    solver.tune();
    const int isign=1;
    double start = omp_get_wtime();
    solver(chi_s[1], psi_s[1], rsd_target, niters, rsd_final, site_flops, mv_apps,isign,verbose);
    double end = omp_get_wtime();
    
    
    
    qdp_unpack_cb_spinor<>(chi_s[1], chi, geom,1);
    
    // Multiply back 
    // chi2 = M chi
    dslash(chi2,u_test,chi, 1, 0);
    invclov_qdp_ap.apply(clov_chi2, chi2, 1, 0);
    dslash(ltmp,u_test,clov_chi2, 1, 1);
    
    clov_qdp_ap.apply(chi2, chi, 1, 1);
    chi2[rb[1]] -= betaFactor*ltmp;
    
  
    Phi diff = chi2 - psi;
    QDPIO::cout << "True norm is: " << sqrt(norm2(diff, rb[1])/norm2(psi,rb[1])) << endl;
    
    int Nxh = Nx/2;
    unsigned long num_cb_sites=Nxh * Ny * Nz * Nt;
    
    unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
    masterPrintf("GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));
  }
#endif

  


  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
  geom.free(A_cb0);
  geom.free(A_cb1);
  geom.free(A_inv_cb0);
  geom.free(A_inv_cb1);

}


void
testClovDslashFull::run(void) 
{
  typedef LatticeColorMatrixF UF;
  typedef LatticeDiracFermionF PhiF;

  typedef LatticeColorMatrixD UD;
  typedef LatticeDiracFermionD PhiD;

#if defined(QPHIX_SCALAR_SOURCE)
  if( precision == FLOAT_PREC ) { 
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if( compress12 ) { 
      runTest<float,1,1,true,UF, PhiF>();
    }
    else { 
      runTest<float,1,1,false,UF, PhiF>();
    }
  }
  if( precision == DOUBLE_PREC ) { 
    QDPIO::cout << "DOUBLE PRECISION TESTING " << endl;
    if( compress12 ) { 
      runTest<double,1,1,true,UF, PhiF>();
    }
    else { 
      runTest<double,1,1,false,UF, PhiF>();
    }
  }
#else




   if( precision == FLOAT_PREC ) { 
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if( compress12 ) { 
      runTest<float,VECLEN_SP,4,true,UF, PhiF>();
    }
    else { 
      runTest<float,VECLEN_SP,4,false,UF, PhiF>();
    }

    if( compress12 ) { 
      runTest<float,VECLEN_SP,8,true,UF, PhiF>();
    }
    else { 
      runTest<float,VECLEN_SP,8,false,UF, PhiF>();
    }

#if defined(QPHIX_MIC_SOURCE)
    if( compress12 ) { 
      runTest<float,VECLEN_SP,16,true,UF, PhiF>();
    }
    else { 
      runTest<float,VECLEN_SP,16,false,UF, PhiF>();
    }
#endif
#endif  // QPHIX_AVX_SOURCE|| QPHIX_MIC_SOURCE

  }


  if( precision == HALF_PREC ) { 
#if defined(QPHIX_MIC_SOURCE)
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if( compress12 ) { 
      runTest<half,VECLEN_HP,4,true,UF, PhiF>();
    }
    else { 
      runTest<half,VECLEN_HP,4,false,UF, PhiF>();
    }

    if( compress12 ) { 
      runTest<half,VECLEN_HP,8,true,UF, PhiF>();
    }
    else { 
      runTest<half,VECLEN_HP,8,false,UF, PhiF>();
    }


    if( compress12 ) { 
      runTest<half,VECLEN_HP,16,true,UF, PhiF>();
    }
    else { 
      runTest<half,VECLEN_HP,16,false,UF, PhiF>();
    }
#else
    QDPIO::cout << "Half precision tests are not available in this build. Currently only in MIC builds" << endl;
#endif

  }

  if( precision == DOUBLE_PREC ) { 
    QDPIO::cout << "DOUBLE PRECISION TESTING" << endl;
    
#if defined(QPHIX_AVX_SOURCE)
    // Only AVX can do DP 2
    if( compress12 ) { 
      runTest<double,VECLEN_DP,2,true,UD, PhiD>();
    }
    else { 
      runTest<double,VECLEN_DP,2,false,UD, PhiD>();
    }
#endif
    
    if( compress12 ) { 
      runTest<double,VECLEN_DP,4,true,UD, PhiD>();
    }
    else { 
      runTest<double,VECLEN_DP,4,false,UD, PhiD>();
    }
    
#if defined(QPHIX_MIC_SOURCE)
    // Only MIC can do DP 8
    if( compress12 ) { 
      runTest<double,VECLEN_DP,8,true,UD, PhiD>();
    }
    else { 
      runTest<double,VECLEN_DP,8,false,UD, PhiD>();
    }
#endif
  }


#endif // SCALAR SOURCE
}
