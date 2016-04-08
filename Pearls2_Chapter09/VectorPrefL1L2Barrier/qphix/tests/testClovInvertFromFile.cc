#include "unittest.h"
#include "testClovInvertFromFile.h"

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
#include "qphix/inv_richardson_multiprec.h"
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
#warning SCALAR_SOURCE
#define VECLEN_DP 1
#define VECLEN_SP 1

#elif defined(QPHIX_QPX_SOURCE)
#define VECLEN_DP 4
#define VECLEN_SP 4

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

    
void MesPlq(const multi1d<LatticeColorMatrixD>& u, Double& w_plaq, Double& link)
{
  w_plaq=link=Double(0);
  
  // Compute the average plaquettes
  for(int mu=1; mu < Nd; ++mu) {
    
    for(int nu=0; nu < mu; ++nu) {

      /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
      LatticeColorMatrix tmp_0 = shift(u[nu],FORWARD,mu) * adj(shift(u[mu],FORWARD,nu));
      
      /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
      LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);
      
      /* tmp = sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu))) */
      Double tmp = sum(real(trace(u[mu]*tmp_1)));
      
      w_plaq += tmp;
    }
  }
  // Normalize
  w_plaq *= 2.0 / double(Layout::vol()*Nd*(Nd-1)*Nc);
  
  // Compute the average link
  for(int mu=0; mu < Nd; ++mu)
    link += sum(real(trace(u[mu])));
  
  link /= double(12*Layout::vol());
}


template<typename FT, int V, int S, bool compress, typename FT2, int V2, int SOA2, typename U, typename Phi>
void
testClovInvertFromFile::runTest(double mass, double clov_coeff, const std::string& filename) 
{
  
  typedef typename Geometry<FT,V,S,compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT,V,S,compress>::CloverBlock Clover;
  
  typedef typename Geometry<FT2,V2,SOA2,compress>::FourSpinorBlock SpinorInner;
  typedef typename Geometry<FT2,V2,SOA2,compress>::SU3MatrixBlock GaugeInner;
  typedef typename Geometry<FT2,V2,SOA2,compress>::CloverBlock CloverInner;
  bool verbose = false;
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


  QDPIO::cout << "Reading field from file" << filename << endl;
  XMLReader file_xml, record_xml;
  QDPFileReader from(file_xml, filename, QDPIO_PARALLEL);
  read(from, record_xml, u);
  close(from);

  Double w_plaq, link;
  MesPlq(u,w_plaq,link);
  QDPIO::cout << "Plaquette on reading: w="<<w_plaq << " link=" << link << endl;
  QDPIO::cout << "Reunitarizing links\n";
  for(int mu=0; mu < 4; mu++) { 
    reunit(u[mu]);
  }
  Double w_plaq2, link2;
  MesPlq(u,w_plaq2, link2);
  QDPIO::cout << "Plaquette after reunit: w="<<w_plaq2 << " link=" << link2 << endl;

  QDPIO::cout << "Delta w_plaq" << w_plaq2-w_plaq << endl;


  // Set anisotropy parameters -- pick some random numbers
  double xi_0_f = 1.0;
  double nu_f = 1.0;
  
  // This is what makeFermCoeffs does under the hood.
  // Get the direct spatial and temporal anisotropy factors
  double aniso_fac_s=(double)(1);
  double aniso_fac_t=(double)(1);

  QDPIO::cout << "Setting Clover Term Parameters" << endl;
  
  CloverFermActParams clparam;
  AnisoParam_t aniso;
  
  // Aniso prarams
  aniso.anisoP=false;
  aniso.xi_0 = 1;
  aniso.nu = 1;
  aniso.t_dir = 3;

  // Set up the Clover params
  clparam.anisoParam = aniso;

  // Some mass
  // Now use the real mass 

  clparam.Mass = Real(mass);
  
  // Some random clover coeffs
  clparam.clovCoeffR=Real(clov_coeff);
  clparam.clovCoeffT=Real(clov_coeff);
  
  double t_boundary=double(-1);
  // Create Dslash
  Geometry<FT,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  Geometry<FT2,V2,SOA2,compress> geom_inner(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);


  // Make a random source
  QDPIO::cout << "Initializing QDP++ input spinor" << endl;
  Phi chi, clov_chi;
  Phi psi, chi2, clov_chi2;
  QDPIO::cout << "Filling psi with gaussian noise" << endl;
  gaussian(psi);
  
  
  QDPIO::cout << "Allocating packged gauge fields" << endl;
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();

  GaugeInner* packed_gauge_cb0_i = (GaugeInner*)geom_inner.allocCBGauge();
  GaugeInner* packed_gauge_cb1_i = (GaugeInner*)geom_inner.allocCBGauge();

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

  CloverInner* A_cb0_i=(CloverInner*)geom_inner.allocCBClov();
  CloverInner* A_cb1_i=(CloverInner*)geom_inner.allocCBClov();
  CloverInner* A_inv_cb0_i=(CloverInner*)geom_inner.allocCBClov();
  CloverInner* A_inv_cb1_i=(CloverInner*)geom_inner.allocCBClov();
  CloverInner* invclov_packed_i[2] = { A_inv_cb0_i, A_inv_cb1_i };
  CloverInner* clov_packed_i[2] = { A_cb0_i, A_cb1_i };

  QDPIO::cout << "Fields allocated" << endl;
 
  
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  qdp_pack_gauge<>(u, packed_gauge_cb0_i,packed_gauge_cb1_i, geom_inner);
 

  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  GaugeInner* u_packed_i[2];
  u_packed_i[0] = packed_gauge_cb0_i;
  u_packed_i[1] = packed_gauge_cb1_i;

  QDPIO::cout << "done" << endl;
  
  QDPIO::cout << " Packing fermions..." ;	
  Spinor *psi_s[2] = { psi_even, psi_odd };
  Spinor *chi_s[2] = { chi_even, chi_odd };

  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
  QDPIO::cout << "done" << endl; 

  QDPIO::cout << "Creating the Clover Term " << endl;

  // Clover term deals with anisotropy internally -- so use original u field.
  QDPCloverTermT<Phi, U> clov_qdp;
  
  QDPIO::cout << "Adding on boundary field" << endl;
  // Modify u (antiperiodic BC's) 
  u[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
		     Real(t_boundary), Real(1));

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
    qdp_pack_clover<>(invclov_qdp.getTriBuffer(), invclov_packed_i[cb], geom_inner, cb);
  }

  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(clov_qdp.getTriBuffer(), clov_packed[cb], geom, cb);
  }
  QDPIO::cout << "Done" << endl;

  for(int cb=0; cb < 2; cb++) { 
    qdp_pack_clover<>(clov_qdp.getTriBuffer(), clov_packed_i[cb], geom_inner, cb);
  }
  QDPIO::cout << "Done" << endl;



  int max_iters=5000;

  EvenOddCloverOperator<FT,V,S,compress> M(u_packed,  
					   clov_packed[1], 
					   invclov_packed[0],  
					   &geom,
					   t_boundary,
					   aniso_fac_s,
					   aniso_fac_t);

  EvenOddCloverOperator<FT2,V2,SOA2,compress> M_inner(u_packed_i,  
						      clov_packed_i[1], 
						      invclov_packed_i[0],  
						      &geom_inner,
						      t_boundary,
						      aniso_fac_s,
						      aniso_fac_t);
  int isign=1;

  InvBiCGStab<FT,V,S,compress> solver(M, max_iters);
  solver.tune();

  InvBiCGStab<FT2,V2,SOA2,compress> solver_inner(M_inner,max_iters);
  solver_inner.tune();

  // Delta = 0.05 for half prec
  InvRichardsonMultiPrec<FT,V,S,compress,FT2,V2,SOA2,compress> mixed_solver(M,solver_inner, 0.1,5000);


  for(int run=0; run < 5; run++) { 
    
    
    Phi ltmp;
    double rsd_target=1.0e-10;
    Real betaFactor=Real(0.25);
    
    int max_iters=5000;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;

    chi = zero;
    qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
    
    double start = omp_get_wtime();
    solver(chi_s[1], psi_s[1], rsd_target, niters, rsd_final, site_flops, mv_apps,isign,verbose);
    double end = omp_get_wtime();
    qdp_unpack_cb_spinor<>(chi_s[1], chi, geom,1);
    
    // Multiply back 
    // chi2 = M chi
    dslash(chi2,u,chi, 1, 0);
    invclov_qdp.apply(clov_chi2, chi2, 1, 0);
    dslash(ltmp,u,clov_chi2, 1, 1);
    
    clov_qdp.apply(chi2, chi, 1, 1);
    chi2[rb[1]] -= betaFactor*ltmp;
    
    Phi diff = chi2 - psi;
    QDPIO::cout << "True norm is: " << sqrt(norm2(diff, rb[1])/norm2(psi,rb[1])) << endl;
    
    int Nxh = Nx/2;
    unsigned long num_cb_sites=Layout::vol()/2;
    unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
    masterPrintf("BICGSTAB: run=%d iters=%d\n", run, niters);
    masterPrintf("BICGSTAB: Time for solve=%16.8e sec\n", (end-start));
    masterPrintf("BICGSTAB: GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));
  }


  for(int run=0; run < 5; run++) {
    chi = zero;
    Phi ltmp;
    double rsd_target=1.0e-10;
    Double betaFactor=Real(0.25);
    
    int max_iters=5000;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;
    qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
    
    double start = omp_get_wtime();
    mixed_solver(chi_s[1], psi_s[1], rsd_target, niters, rsd_final, site_flops, mv_apps,isign,verbose);
    double end = omp_get_wtime();
    
    
    qdp_unpack_cb_spinor<>(chi_s[1], chi, geom,1);
    
    // Multiply back 
    // chi2 = M chi
    dslash(chi2,u,chi, 1, 0);
    invclov_qdp.apply(clov_chi2, chi2, 1, 0);
    dslash(ltmp,u,clov_chi2, 1, 1);
    
    clov_qdp.apply(chi2, chi, 1, 1);
    chi2[rb[1]] -= betaFactor*ltmp;
    
    
    Phi diff = chi2 - psi;
    QDPIO::cout << "True norm is: " << sqrt(norm2(diff, rb[1])/norm2(psi,rb[1])) << endl;
    
    int Nxh = Nx/2;
    unsigned long num_cb_sites=Layout::vol()/2;
    
    unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
    masterPrintf("RICHARDSON: run=%d iters=%d\n", run, niters);
    masterPrintf("RICHARDSON Time for solve=%16.8e sec\n", (end-start));
    masterPrintf("RICHARDSON GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));
  }

  
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

  geom_inner.free(packed_gauge_cb0_i);
  geom_inner.free(packed_gauge_cb1_i);

  geom_inner.free(A_cb0_i);
  geom_inner.free(A_cb1_i);
  geom_inner.free(A_inv_cb0_i);
  geom_inner.free(A_inv_cb1_i);

  
}


void
testClovInvertFromFile::run(void) 
{
  typedef LatticeColorMatrixF UF;
  typedef LatticeDiracFermionF PhiF;

  typedef LatticeColorMatrixD UD;
  typedef LatticeDiracFermionD PhiD;


  // Kappa 0.13632 => Mass is -0.3321595
  // Kappa 0.1364  => Mass is -0.3343108

  // 32^3x 64 run
#if 0

  std::string filename("./qcdsf.632.01111.lime");
  double mass=-0.3321595;
  double clov_coeff=1.9192;
#if defined(QPHIX_MIC_SOURCE)
  runTest<double,VECLEN_DP,8,true,half,VECLEN_HP,16,UD, PhiD>(mass,clov_coeff,filename);
#else
  runTest<double,VECLEN_DP,4,true,float,VECLEN_SP,8,UD, PhiD>(mass,clov_coeff,filename);
#endif
#endif

#if 0
  std::string filename("./qcdsf.743.00600.lime");
  double mass=-0.3343108;
  double clov_coeff=1.9192;
#if defined(QPHIX_MIC_SOURCE)
  runTest<double,VECLEN_DP,8,true,half,VECLEN_HP,8,UD, PhiD>(mass,clov_coeff,filename);
#else
  runTest<double,VECLEN_DP,4,true,float,VECLEN_SP,8,UD, PhiD>(mass,clov_coeff,filename);
#endif
#endif

  std::string filename("./cl3_64_128_b5p0_m0p3550_m0p3550_cfg_544.lime");
  double mass=-0.3550;
  double clov_coeff=1.90497469553511;
#if defined(QPHIX_MIC_SOURCE)
  runTest<double,VECLEN_DP,8,true,half,VECLEN_HP,16,UD, PhiD>(mass,clov_coeff,filename);
#else
  runTest<double,VECLEN_DP,4,true,float,VECLEN_SP,8,UD, PhiD>(mass,clov_coeff,filename);
#endif


}
