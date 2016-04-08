#include <iostream>
using namespace std;

#include "unittest.h"
#include "testDslashFull.h"
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/blas_new_c.h"
// Disabling Full M and CG tests until vectorized dslash 
// works better
#include "qphix/wilson.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"
#if 1
#include "./invbicgstab_test.h"
#endif

#include <omp.h>
#if 0
#include "qphix/memmap.h"
#endif

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
#include <immintrin.h>

#elif defined(QPHIX_AVX_SOURCE) 

#define VECLEN_SP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SCALAR_SOURCE)
#warning SCALAR_SOURCE
#define VECLEN_SP 1
#define VECLEN_DP 1

#elif defined(QPHIX_QPX_SOURCE)
#warning QPX_SOURCE
#define VECLEN_SP 4
#define VECLEN_DP 4

#endif

  // What we consider to be small enough...
int Nx, Ny, Nz, Nt, Nxh;
bool verbose = true;

template<typename F> 
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
const double rsdTarget<half>::value = (double)(1.0e-4);

template<>
const double rsdTarget<float>::value = (double)(1.0e-7);


template<>
const double rsdTarget<double>::value = (double)(1.0e-12);


void
testDslashFull::run(void) 
{
  RNG::savern(rng_seed);

  typedef multi1d<LatticeColorMatrixF> UF;
  typedef multi1d<LatticeColorMatrixD> UD;
  typedef LatticeDiracFermionF PhiF;
  typedef LatticeDiracFermionD PhiD;

  // Diagnostic information:
  const multi1d<int>& lattSize = Layout::subgridLattSize();
  Nx = lattSize[0];
  Ny = lattSize[1];
  Nz = lattSize[2];
  Nt = lattSize[3];

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

 


  QDPIO::cout << "Inititalizing QDP++ gauge field"  << endl;
  // Make a random gauge field 
  multi1d<LatticeColorMatrix> u(4);
  LatticeColorMatrix g;
  LatticeColorMatrix uf;
  for(int mu=0; mu < 4; mu++) { 
    uf = 1;   // Unit gauge

    Real factor=Real(0.09);
    gaussian(g);
    u[mu] = uf + factor*g;
    reunit(u[mu]);
  }



#if 1 // Save build time
  if( precision == FLOAT_PREC ) {
    
    QDPIO::cout << "SINGLE PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    {

      if( soalen == 1 ) { 
#if defined(QPHIX_SCALAR_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=1 " << endl;
	testDslashWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
//	testDslashAChiMBDPsiWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
//	testMWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
//	testCGWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
//	testBiCGStabWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
#endif
      }

      if( soalen == 4 ) { 
#if defined (QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=4 " << endl;
	testDslashWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
	testMWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
	testCGWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
	testBiCGStabWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
#endif
      }


      if( soalen == 8 ) {
#if defined (QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=8"  << endl;
	testDslashWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
	testMWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
	testCGWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
	testBiCGStabWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
#endif
      }

      if ( soalen == 16 ) { 
#if defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=16 " << endl;
	testDslashWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
	testMWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
	testCGWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
	testBiCGStabWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
#else 
	masterPrintf("SOALEN=16 not available");
	return;
#endif

      }


    }
  }

#endif // If 0

#if 1
  if (precision == HALF_PREC ) { 
#if defined(QPHIX_MIC_SOURCE)
    QDPIO::cout << "HALF PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    {
      if( soalen == 4 ) { 
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=4 " << endl;
	testDslashWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
	testMWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
	testCGWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
	testBiCGStabWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
      }

      if (soalen == 8 ) {
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=8 " << endl;
	testDslashWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
	testMWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
	testCGWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
	testBiCGStabWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
      }

      if( soalen == 16 ) {
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=16 " << endl;
	testDslashWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
	testDslashAChiMBDPsiWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
	testMWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
	testCGWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
	testBiCGStabWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
      }
    }
#else
    QDPIO::cout << " Half Prec is only supported on MIC Targets just now " << endl;
#endif
  }
#endif // If 0

#if 1
  if( precision == DOUBLE_PREC ) { 
    QDPIO::cout << "DOUBLE PRECISION TESTING:" << endl;
    UD u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    
    {

      if( soalen == 1) {
#if defined (QPHIX_SCALAR_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=1 " << endl;
	testDslashWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
	testDslashAChiMBDPsiWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
	testMWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
	testCGWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
	testBiCGStabWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
	
#endif
      }


      if( soalen == 2) {
#if defined (QPHIX_AVX_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=2 " << endl;
	testDslashWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
	testDslashAChiMBDPsiWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
	testMWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
	testCGWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
	testBiCGStabWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
	
#endif
      }

      if( soalen == 4 ) { 
#if defined (QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=4 " << endl;
	testDslashWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testDslashAChiMBDPsiWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testMWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testCGWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testBiCGStabWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);

#elif defined (QPHIX_QPX_SOURCE) 
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=4 " << endl;
	testDslashWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testDslashAChiMBDPsiWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testMWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testCGWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testBiCGStabWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
#endif
      }

      if( soalen == 8 ) { 
#if defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=8 " << endl;
	testDslashWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
	testDslashAChiMBDPsiWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
	testMWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
	
	testCGWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
	testBiCGStabWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
#endif
      }
      
    }
  }
#endif // If 0

#if 0
  {
    multi1d<LatticeColorMatrixD3> u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

    int t_bc=-1;
    

#if defined(QPHIX_MIC_SOURCE)
    if ( Nx % 32 == 0 ){
      testBiCGStabWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
      testRichardson<double,VECLEN_DP,8,true,half,VECLEN_HP,16,UD,PhiD>(u_in, t_bc);
    }
    else { 
      if ( Nx % 16 == 0) {
	testBiCGStabWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
	testRichardson<double,VECLEN_DP,8,true,half,VECLEN_HP,8,UD,PhiD>(u_in, t_bc);
      }
      else {
	masterPrintf("I havent set up that mixed precision solver combination\n");
      }
    }
#else 
    // AVX: Double SOALEN = 4
    if ( Nx % 16 == 0 ){
      testBiCGStabWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
      testRichardson<double,VECLEN_DP,4,true,float,VECLEN_SP,8,UD,PhiD>(u_in, t_bc);
    }
    else { 
      if ( Nx % 8 == 0) {
	testBiCGStabWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
	testRichardson<double,VECLEN_DP,4,true,float,VECLEN_SP,4,UD,PhiD>(u_in, t_bc);
      }
      else {
	masterPrintf("I havent set up that mixed precision solver combination\n");
      }
    }
#endif


  }
#endif

}



template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
testDslashFull::testDslash(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock Spinor;
  QDPIO::cout << "In testDslash" << endl;
  double aniso_fac_s = ((double)0.35);
  double aniso_fac_t = ((double)1.4);
  double t_boundary = ((double)t_bc);
  QDPIO::cout << "Testing Dslash...  T BCs = " << t_boundary << endl;

  Phi chi,chi2;
  Phi psi;
 
  QDPIO::cout << "Filling psi with noise: " << endl;  
  gaussian(psi);


  

#if 0
  Dslash<T,V,S,compress> D32(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt, t_boundary, aniso_fac_s, aniso_fac_t);
  
  Geometry<T,V,S,compress>& geom= D32.getGeometry();
#endif
  
  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  Dslash<T,V,S,compress> D32(&geom, t_boundary, aniso_fac_s, aniso_fac_t);
  
  Gauge* packed_gauge_cb0 = ( Gauge *)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = ( Gauge *)geom.allocCBGauge();
  Spinor* psi_even=( Spinor*)geom.allocCBFourSpinor();
  Spinor* psi_odd=( Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_even=( Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_odd=( Spinor*)geom.allocCBFourSpinor();
  
  
  QDPIO::cout << "Fields allocated" << endl;
 
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);

 
  QDPIO::cout << "done" << endl;
 


  
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };
  
  QDPIO::cout << " Packing fermions..." ;	
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
    
  QDPIO::cout << "done" << endl; 


  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  QDPIO::cout << "U field prepared" <<endl;
  
  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
		     Real(t_boundary), Real(1));

  QDPIO::cout << "BCs applied" << endl;
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  int isign=1;
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;

      chi = zero;
      
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
      
      // Apply Optimized Dslash
      D32.dslash(chi_s[target_cb],	
		 psi_s[source_cb],
		 u_packed[target_cb],
		 isign, 
		 target_cb);

      //      qdp_unpack_spinor<T,V,S,compress, Phi >(chi_even,chi_odd, chi, geom);
      qdp_unpack_spinor<>(chi_even,chi_odd, chi, geom);
    
      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, target_cb);
 
      // Check the difference per number in chi vector
      Phi diff = chi2 -chi;
      
      Double diff_norm = sqrt( norm2( diff, rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
      
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm > tolerance<T>::small ) ) {
	
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
		    if( toBool( fabs(dr) > tolerance<T>::small ) || toBool ( fabs(di) > tolerance<T>::small ) ) {
		      QDPIO::cout <<"(x,y,z,t)=(" << x <<"," <<y<<","<<z<<","<<t<<") site=" << ind << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+ind).elem(s).elem(c) 
				  << "  chi = " << chi.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << " qdp++ =" << chi2.elem(rb[target_cb].start()+ind).elem(s).elem(c)  << endl;
		    }
		  }
		}
	      } // x 
	    } // y 
	  } // z 
	} // t
	assertion( toBool( diff_norm <= tolerance<T>::small ) );
      }
	
      
    } // cb
  } // isign



  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);

}




template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
testDslashFull::testDslashAChiMBDPsi(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock   Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock  Spinor;

  QDPIO::cout << "In testAChiMBDPsi " << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);

  Phi psi, chi, chi2;
  QDPIO::cout << "Filling psi with random noise" << endl;
  gaussian(psi);

#if 0
  Dslash<T,V,S,compress> D32(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt, t_boundary, aniso_fac_s, aniso_fac_t);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT 
  Geometry<T,V,S,compress>& geom= D32.getGeometry();
#endif

  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT 
  Dslash<T,V,S,compress> D32(&geom,t_boundary, aniso_fac_s, aniso_fac_t);
  
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  // Over allocate, so that an unsigned load doesnt cause segfault accesing either end...
 Spinor* psi_even=(Spinor*)geom.allocCBFourSpinor();
 Spinor* psi_odd=(Spinor*)geom.allocCBFourSpinor();
 Spinor* chi_even=(Spinor*)geom.allocCBFourSpinor();
 Spinor* chi_odd=(Spinor*)geom.allocCBFourSpinor();


  QDPIO::cout << "Fields allocated" << endl;
 
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };
  
  QDPIO::cout << " Packing fermions..." ;	
  //  qdp_pack_spinor< T,V,S,compress, Phi >(psi, psi_even, psi_odd, geom);
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
    
  QDPIO::cout << "done" << endl; 

  QDPIO::cout << "Testing DslashAChiMBDPsi \n" << endl;
  QDPIO::cout << "T BCs = " << t_boundary << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }

  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi=zero;
      //      qdp_pack_spinor< T,V,S,compress,Phi >(chi, chi_even, chi_odd, geom);
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      double alpha = (double)(4.01); // Nd + M, M=0.01
      double beta = (double)(0.5);   // Operator is (Nd+M) - (1/2)Dslash
      

      // Apply Optimized Dslash
      D32.dslashAChiMinusBDPsi(chi_s[target_cb],	
			       psi_s[source_cb],
			       psi_s[target_cb],
			       u_packed[target_cb],
			       alpha, 
			       beta,
			       isign, 
			       target_cb);
      
      // qdp_unpack_spinor< T,V,S,compress,Phi>(chi_s[0], chi_s[1], chi, geom);
      qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);


      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, target_cb);
      Phi res;
      res[rb[source_cb]]=zero;
      res[rb[target_cb]]= alpha*psi - beta*chi2;

      // Check the difference per number in chi vector
      Phi diff = res-chi;

      Double diff_norm = sqrt( norm2( diff , rb[target_cb] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      if ( toBool( diff_norm >= tolerance<T>::small ) ) {
	for(int i=0; i < rb[target_cb].siteTable().size(); i++){ 
	  for(int s =0 ; s < Ns; s++) { 
	    for(int c=0; c < Nc; c++) { 
	      	      QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[target_cb].start()+i).elem(s).elem(c) << endl;
	    }
	  }
	}
      }
      assertion( toBool( diff_norm < tolerance<T>::small ) );

    }
   }
    

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);

}

template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
testDslashFull::testM(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock Spinor;

  QDPIO::cout << "in testM:" << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);


  Phi psi,chi,chi2;
  QDPIO::cout << "Filling source spinor with gaussian noise" << endl;
  gaussian(psi);

  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  // Over allocate, so that an unsigned load doesnt cause segfault accesing either end...
  Spinor* psi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* psi_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_odd=(Spinor*)geom.allocCBFourSpinor();


  QDPIO::cout << "Fields allocated" << endl;

  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };

  QDPIO::cout << " Packing fermions..." ;	
  //  qdp_pack_spinor< T,V,S,compress, Phi >(psi, psi_even, psi_odd, geom);
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
    
  QDPIO::cout << "done" << endl; 

  QDPIO::cout << "Testing M \n" << endl;
  QDPIO::cout << "T BCs = " << t_boundary << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  double Mass=0.1;
  EvenOddWilsonOperator<T,V,S,compress> M(Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);



  Phi ltmp=zero;
  Real massFactor=Real(4) + Real(Mass);
  Real betaFactor=Real(0.25)/massFactor;
   // Apply optimized
  for(int isign=1; isign >= -1; isign -=2) {
    
      chi=zero;
      //      qdp_pack_spinor< T,V,S, compress, Phi >(chi, chi_even, chi_odd, geom);
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      M( chi_s[0],
	 psi_s[0],
	 isign);

      //      qdp_unpack_spinor< T,V,S, compress,  Phi> (chi_s[0], chi_s[1], chi, geom);
      qdp_unpack_spinor<> (chi_s[0], chi_s[1], chi, geom);


      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2,u_test,psi, isign, 1);
      dslash(ltmp,u_test,chi2, isign, 0);
      chi2[rb[0]] = massFactor*psi - betaFactor*ltmp;

      // Check the difference per number in chi vector
      Phi diff = zero;
      diff[rb[0]]=chi2-chi;

      Double diff_norm = sqrt( norm2( diff , rb[0] ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t isign = " << isign << "  diff_norm = " << diff_norm << endl;      
     // Assert things are OK...
      if ( toBool( diff_norm >= tolerance<T>::small ) ) {
	for(int i=0; i < rb[0].siteTable().size(); i++){ 
	  for(int s =0 ; s < Ns; s++) { 
	    for(int c=0; c < Nc; c++) { 
	      Double re=  Double(diff.elem(rb[0].start()+i).elem(s).elem(c).real());
	      Double im=  Double(diff.elem(rb[0].start()+i).elem(s).elem(c).imag());
	      if( toBool ( fabs(re) > tolerance<T>::small) || toBool( fabs(im) > tolerance<T>::small ) ) { 
		QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c << " Diff = " << diff.elem(rb[0].start()+i).elem(s).elem(c) << endl;
	      }

	    }
	  }
	}
      }
      assertion( toBool( diff_norm < tolerance<T>::small ) );

    
  }



  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);

}


template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
testDslashFull::testCG(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock   Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock  Spinor;

  QDPIO::cout << "In testCG:" << endl;

  Phi psi, chi, chi2,chi3;
  QDPIO::cout << "Filling psi with gaussian noise" << endl;

  gaussian(psi);
  QDPIO::cout << "Norm2 || psi || = " << norm2(psi,rb[0]) << endl;
  QDPIO::cout << "Done" << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);

  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT 
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  // Over allocate, so that an unsigned load doesnt cause segfault accesing either end...
  Spinor* psi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* psi_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_odd=(Spinor*)geom.allocCBFourSpinor();


  QDPIO::cout << "Fields allocated" << endl;
 
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };

  QDPIO::cout << " Packing fermions..." ;	

  // qdp_pack_spinor< T,V,S,compress, Phi >(psi, psi_even, psi_odd, geom);
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
    
  QDPIO::cout << "done" << endl; 

  QDPIO::cout << "T BCs = " << t_boundary << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  double Mass=0.1;
  EvenOddWilsonOperator<T,V,S,compress> M(Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
  Phi ltmp=zero;
  Real massFactor=Real(4) + Real(Mass);
  Real betaFactor=Real(0.25)/massFactor;

  {
    chi = zero;
    //    qdp_pack_spinor<T,V,S,compress, Phi >(chi, chi_even, chi_odd, geom);
    qdp_pack_cb_spinor<>(chi, chi_even, geom,0);
    
    double rsd_target=rsdTarget<T>::value;
    int max_iters=200;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;

    InvCG<T,V,S,compress> solver(M, max_iters);
    solver.tune();
    double r2;
    norm2Spinor<T,V,S,compress>(r2,psi_even,geom,1);
    masterPrintf("psi has norm2 = %16.8e\n", r2);
    int isign=1;
    double start = omp_get_wtime();
    solver(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps, isign, verbose);
    double end = omp_get_wtime();
    
    
    
    //    qdp_unpack_spinor<T, V, S, compress,Phi >(chi_s[0], chi_s[1], chi, geom);
    qdp_unpack_cb_spinor<>(chi_s[0], chi, geom, 0);
    
    // Multiply back 
    // chi2 = M chi
    dslash(chi2,u_test,chi, isign, 1);
    dslash(ltmp,u_test,chi2, isign, 0);
    chi2[rb[0]] = massFactor*chi - betaFactor*ltmp;
    
    // chi3 = M^\dagger chi2
    dslash(chi3,u_test,chi2, (-isign), 1);
    dslash(ltmp,u_test,chi3, (-isign), 0);
    chi3[rb[0]] = massFactor*chi2 - betaFactor*ltmp;
    
    Phi diff = chi3 - psi;
    Double true_norm = sqrt(norm2(diff, rb[0])/norm2(psi,rb[0]));
    QDPIO::cout << "True norm is: " << true_norm << endl;
    assertion( toBool( true_norm < (rsd_target + tolerance<T>::small) ) );

    unsigned long num_cb_sites=Layout::vol()/2;
    unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;

    masterPrintf("GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));

  }



  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);

}


template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
testDslashFull::testBiCGStab(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock   Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock  Spinor;

  QDPIO::cout << "In testBiCGStab" << endl;

  //  double aniso_fac_s = (double)(0.35);
  //double aniso_fac_t = (double)(1.4);
  double aniso_fac_s = (double)(1);
  double aniso_fac_t = (double)(1);

  double t_boundary = (double)(t_bc);

  Phi psi, chi, chi2;
  gaussian(psi);

  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT 

  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  // Over allocate, so that an unsigned load doesnt cause segfault accesing either end...
  Spinor* psi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* psi_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* chi_odd=(Spinor*)geom.allocCBFourSpinor();


  QDPIO::cout << "Fields allocated" << endl;
 
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };

  QDPIO::cout << " Packing fermions..." ;	
  //  qdp_pack_spinor< T,V,S,compress, Phi >(psi, psi_even, psi_odd, geom);
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
    
  QDPIO::cout << "done" << endl; 
  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  double Mass=0.1;
  EvenOddWilsonOperator<T,V,S,compress> M(Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
  Phi ltmp=zero;
  Real massFactor=Real(4) + Real(Mass);
  Real betaFactor=Real(0.25)/massFactor;

  {
    float rsd_target=rsdTarget<T>::value;
    int max_iters=100;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;
    
    InvBiCGStab<T,V,S,compress> solver(M, max_iters);
    solver.tune();

    for(int isign=1; isign >= -1; isign -=2) {
      chi=zero;
      qdp_pack_cb_spinor<>(chi, chi_even, geom,0);
      masterPrintf("Entering solver\n");
      
      double start = omp_get_wtime();
      solver(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps,isign,verbose);
      double end = omp_get_wtime();
      
      qdp_unpack_cb_spinor<>(chi_s[0], chi, geom,0);
      
      // Multiply back 
      // chi2 = M chi
      dslash(chi2,u_test,chi, isign, 1);
      dslash(ltmp,u_test,chi2, isign, 0);
      chi2[rb[0]] = massFactor*chi - betaFactor*ltmp;
      
    
      Phi diff = chi2 - psi;
      Double true_norm =  sqrt(norm2(diff, rb[0])/norm2(psi,rb[0]));
      QDPIO::cout << "BiCGStab Solve isign=" << isign << " True norm is: " << true_norm << endl;
      assertion( toBool( true_norm < (rsd_target + tolerance<T>::small) ) );
      int Nxh = Nx/2;
      unsigned long num_cb_sites=Layout::vol()/2;
      unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;
      masterPrintf("BICGSTAB Solve isign=%d, iters=%d MV apps=%lu Site flops=%lu\n",isign, niters, mv_apps, site_flops);
      masterPrintf("BICGSTAB Solve isign=%d, time=%e (sec) GFLOPS=%e\n", isign, (end-start), 1.0e-9*(double)(total_flops)/(end -start));
    }

  }

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);

}



template<typename T1, int VEC1, int SOA1, bool compress, typename T2, int VEC2, int SOA2, typename U, typename Phi>
void 
testDslashFull::testRichardson(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T1,VEC1,SOA1,compress>::SU3MatrixBlock   GaugeOuter;
  typedef typename Geometry<T1,VEC1,SOA1,compress>::FourSpinorBlock  SpinorOuter;

  typedef typename Geometry<T2,VEC2,SOA2,compress>::SU3MatrixBlock   GaugeInner;
  typedef typename Geometry<T2,VEC2,SOA2,compress>::FourSpinorBlock  SpinorInner;

  QDPIO::cout << "In testRichardson" << endl;

  double aniso_fac_s = (double)(1);
  double aniso_fac_t = (double)(1);
  double t_boundary = (double)(t_bc);

  Phi psi, chi, chi2;
  gaussian(psi);

  Geometry<T1,VEC1,SOA1,compress> geom_outer(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

  Geometry<T2,VEC2,SOA2,compress> geom_inner(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT 

  GaugeOuter* packed_gauge_cb0 = (GaugeOuter*)geom_outer.allocCBGauge();
  GaugeOuter* packed_gauge_cb1 = (GaugeOuter*)geom_outer.allocCBGauge();

  QDPIO::cout << "Packing gauge field..." ;
  //  qdp_pack_gauge< T,V,S,compress, U >(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom_outer);

  GaugeOuter * u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  QDPIO::cout << "done" << endl;


  GaugeInner* packed_gauge_cb0_inner = (GaugeInner*)geom_inner.allocCBGauge();
  GaugeInner* packed_gauge_cb1_inner = (GaugeInner*)geom_inner.allocCBGauge();

  QDPIO::cout << "Packing inner gauge field..." ;
  qdp_pack_gauge<>(u, packed_gauge_cb0_inner,packed_gauge_cb1_inner, geom_inner);
  GaugeInner * u_packed_inner[2];
  u_packed_inner[0] = packed_gauge_cb0_inner;
  u_packed_inner[1] = packed_gauge_cb1_inner;
  QDPIO::cout << "done" << endl;

  // Over allocate, so that an unsigned load doesnt cause segfault accesing either end...
  SpinorOuter* psi_even=(SpinorOuter*)geom_outer.allocCBFourSpinor();
  SpinorOuter* psi_odd=(SpinorOuter*)geom_outer.allocCBFourSpinor();
  SpinorOuter* chi_even=(SpinorOuter*)geom_outer.allocCBFourSpinor();
  SpinorOuter* chi_odd=(SpinorOuter*)geom_outer.allocCBFourSpinor();

  QDPIO::cout << "Fields allocated" << endl;
 
  SpinorOuter *psi_s[2] = { psi_even, psi_odd };
  SpinorOuter *chi_s[2] = { chi_even, chi_odd };

  QDPIO::cout << " Packing fermions..." ;	
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom_outer);
    
  QDPIO::cout << "done" << endl; 
  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  double Mass=0.001;
  EvenOddWilsonOperator<T1,VEC1,SOA1,compress> M_outer(Mass, u_packed, &geom_outer, t_boundary, aniso_fac_s, aniso_fac_t);

  EvenOddWilsonOperator<T2,VEC2,SOA2,compress> M_inner(Mass, u_packed_inner, &geom_inner, t_boundary, aniso_fac_s, aniso_fac_t);


  Phi ltmp=zero;
  Real massFactor=Real(4) + Real(Mass);
  Real betaFactor=Real(0.25)/massFactor;

  {
    float rsd_target=rsdTarget<T1>::value;
    int max_iters=200;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;

    InvBiCGStab<T2,VEC2,SOA2,compress> inner_solver(M_inner, max_iters,isign);    
    InvRichardsonMultiPrec<T1,VEC1,SOA1,compress,T2,VEC2,SOA2,compress> outer_solver(M_outer,inner_solver,0.01,isign,max_iters);
    
    for(int isign=1; isign >= -1; isign -=2) {
      // BiCGStab Inner SOlver

    
      chi = zero;
      //      qdp_pack_spinor<T,V,S, compress,Phi >(chi, chi_even, chi_odd, geom);
      qdp_pack_cb_spinor<>(chi, chi_even, geom_outer,0);
      masterPrintf("Entering solver\n");
      
      double start = omp_get_wtime();
      outer_solver(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps,isign, verbose);
      double end = omp_get_wtime();
      
      //      qdp_unpack_spinor<T,V,S, compress, Phi >(chi_s[0], chi_s[1], chi, geom);
      qdp_unpack_cb_spinor<>(chi_s[0], chi, geom_outer,0);
      
      // Multiply back 
      // chi2 = M chi
      dslash(chi2,u_test,chi, isign, 1);
      dslash(ltmp,u_test,chi2, isign, 0);
      chi2[rb[0]] = massFactor*chi - betaFactor*ltmp;
      
    
      Phi diff = chi2 - psi;
      Double true_norm =  sqrt(norm2(diff, rb[0])/norm2(psi,rb[0]));
      QDPIO::cout << "RICHARDSON Solve isign=" << isign << " True norm is: " << true_norm << endl;
      assertion( toBool( true_norm < (rsd_target + tolerance<T1>::small) ) );
      unsigned long num_cb_sites=Layout::vol()/2;
      unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;
      
      masterPrintf("RICHARDSON Solve isign=%d time=%e (sec) GFLOPS=%e\n", isign, (end - start), 1.0e-9*(double)(total_flops)/(end -start));

    }

  }

  geom_outer.free(packed_gauge_cb0);
  geom_outer.free(packed_gauge_cb1);

  geom_inner.free(packed_gauge_cb0_inner);
  geom_inner.free(packed_gauge_cb1_inner);

  geom_outer.free(psi_even);
  geom_outer.free(psi_odd);
  geom_outer.free(chi_even);
  geom_outer.free(chi_odd);

}

