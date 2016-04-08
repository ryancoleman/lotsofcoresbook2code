#ifndef QPHIX_GEOMETRY_H
#define QPHIX_GEOMETRY_H

#include "qphix/dslash_utils.h"
#include "qphix/print_utils.h"

#include <cstdlib>
#include <iostream>


using namespace std;

namespace QPhiX {

  struct CorePhase {
    int Ct;
    int Cyz;
    int startBlock;
  }; 

  typedef unsigned short half;
  
#if defined(QPHIX_MIC_SOURCE)
  float cvtHalf2Float(half val) {
    float ret;
    _mm512_mask_packstorelo_ps(&ret, 0x1, _mm512_mask_extloadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val, _MM_UPCONV_PS_FLOAT16, _MM_HINT_NONE));
    return ret;
  }
  
  half cvtFloat2Half(float val) {
    half ret;
    _mm512_mask_extpackstorelo_ps(&ret, 0x1, _mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val), _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NONE);
    return ret;
  }
  

  // rep: cast 'in' of type T2 into a 'T1' and return it. 
  template <typename T1, typename T2>
  T1 rep(const T2& in) 
  {
    if(sizeof(T1) != sizeof(T2)) {
      
      if(sizeof(T1) == 2) // we are converting float/double to half
	return cvtFloat2Half((float)in);
      else if(sizeof(T2) == 2) // we are converting half to float/double
	return (T1)cvtHalf2Float(in);
      else
	return (T1)in; // we are converting between float and double so just cast is enough
    }
    else {
      return static_cast<T1>(in);  // both T1 and T2 are same
    }
  }
  
#else 

  // rep: cast 'in' of type T2 into a 'T1' and return it. 
  template <typename T1, typename T2>
  T1 rep(const T2& in) 
  {
    return (T1)(in); 
  }
 

#endif


  template<typename T, int V, int S, const bool compressP>
  class Geometry {
  public:
    // Later change this to depend on compressP
    typedef T  FourSpinorBlock[3][4][2][S];
    typedef T  TwoSpinorBlock[3][2][2][V];
    typedef T  SU3MatrixBlock[8][ ( compressP ? 2 : 3 ) ][3][2][V];


    struct CloverBlock {
      T diag1[6][V];         // Real Diagonal part of block 1
      T off_diag1[15][2][V];  // Complex, off diagonal part of block 1
      T diag2[6][V];          // Real Diagonal part of block 2
      T off_diag2[15][2][V];  // Complex, off diagonal part of block 2
    };

    Geometry(const int latt_size[],
	     int By_,
	     int Bz_,
	     int NCores_,
	     int Sy_,
	     int Sz_,
	     int PadXY_,
	     int PadXYZ_,
	     int MinCt_)
      : Nd(4),  By(By_), Bz(Bz_), num_cores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), nsimt(Sy_*Sz_),  num_threads(NCores_*Sy_*Sz_)
    {   
      Nx_ = latt_size[0];
      Ny_ = latt_size[1];
      Nz_ = latt_size[2];
      Nt_ = latt_size[3];
      Nxh_ = Nx_/2;
      
      nvecs_ = Nxh()/ S;
      if (Nxh()% S != 0) nvecs_++;
      
      if ( V % S != 0 ) { 
	cerr << "Error: Geometry constructor: SOALEN="<< S <<" does not divide V=" << V << endl;
	abort();
      }
      ngy_ = V/S;
      
      // Padding constants
      Pxy = (nvecs_*Ny_+ PadXY);
      Pxyz = (Pxy*Nz_+ PadXYZ);
      
       
      // Allos sizes 
      spinor_bytes = (Pxyz * Nt_ + 1)*sizeof(FourSpinorBlock);
      gauge_bytes = ((Pxyz*Nt_*S)/V)*sizeof(SU3MatrixBlock);
      clover_bytes = ((Pxyz*Nt_*S)/V)*sizeof(CloverBlock);

      // This works out the phase breakdown
      int ly = Ny_ / By;
      int lz = Nz_ / Bz;
      int rem = ly * lz;
      int stblk = 0;
      n_phases = 0;
      int n_cores_per_minct = num_cores / MinCt;
      while(rem > 0) {
	int ctd = n_cores_per_minct / rem;
	int ctu = (n_cores_per_minct + rem - 1) / rem;
	CorePhase& p = getCorePhase(n_phases);
	p.Ct = (ctu <= 4 ? ctu : ctd)*MinCt;
	p.Cyz = num_cores / p.Ct;
	if(p.Cyz > rem) p.Cyz = rem;
	p.startBlock = stblk;
	stblk += p.Cyz;
	rem -= p.Cyz;
	//	masterPrintf("Phase %d: Cyz = %d Ct = %d, start = %d\n", n_phases, p.Cyz, p.Ct, p.startBlock);
	n_phases++;
      }
    }
    
    
    ~Geometry() {}

    inline   int Nxh() const   { return Nxh_; } // Keep
    inline   int Nx()  const   { return Nx_; } // Keep
    inline   int Ny()  const   { return Ny_; } // Keep
    inline   int Nz()  const   { return Nz_; } // Keep
    inline   int Nt()  const   { return Nt_; }  //Keep
    inline int nVecs() const { return nvecs_; }
    inline int nGY() const { return ngy_; }


    /*! \brief Checkerboarded FourSpinor Allocator
     *
     * Allocates a single checkerboard of a Four Spinor.
     * An extra spinor is allocated beyond what is required.
     */
    FourSpinorBlock* allocCBFourSpinor()
    {
            
      FourSpinorBlock *ret_val = (FourSpinorBlock *)BUFFER_MALLOC(spinor_bytes, 128);
      if ( ret_val == (FourSpinorBlock *)0x0 ) { 
	masterPrintf("Failed to allocate FourSpinorBlock\n");
	abort();
      }

      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;
      
      // change from number of bytes to number of T type elements
      size_t num_ft = spinor_bytes / sizeof(T);
      
      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
#if defined (__INTEL_COMPILER)      
#pragma simd 
#pragma vector nontemporal(ret_val_ft)
#endif
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] =rep<T,double>(0.0);
      }
      
      return ret_val+1;
    }


    
    void free(FourSpinorBlock* p) 
    {
      FourSpinorBlock* freeme=p-1;
      BUFFER_FREE(freeme,spinor_bytes);
    }
  
  /*! \brief Checkerboard Gauge Field Allocation
   *
   * This function allocates memory for a single checkerboard of 
   * a gauge field
   */
    SU3MatrixBlock* allocCBGauge()
    {
      SU3MatrixBlock *ret_val = (SU3MatrixBlock *)BUFFER_MALLOC(gauge_bytes, 128);
      if ( ret_val == (SU3MatrixBlock *)0x0 ) { 
	masterPrintf("Failed to allocate SU3MatrixBlock\n");
	abort();
      }

      // For AVX we should loop and zero it here....
      // later on.
      
      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;
      
      // change from number of bytes to number of T type elements
      size_t num_ft = gauge_bytes / sizeof(T);
      
      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector nontemporal(ret_val_ft)
#endif
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] = rep<T,double>(0.0);
      }
      
      return ret_val;
    }

    void free(SU3MatrixBlock *p)
    {
      BUFFER_FREE(p, gauge_bytes);
    }

    CloverBlock* allocCBClov()
    {
      CloverBlock *ret_val = (CloverBlock *)BUFFER_MALLOC(clover_bytes, 128);
      if ( ret_val == (CloverBlock *)0x0 ) { 
	masterPrintf("Failed to allocate CloverBlock\n");
	abort();
      }

      // For AVX we should loop and zero it here....
      // later on.

      // Zero the field.
      // Cast the pointer.
      T *ret_val_ft = (T *)ret_val;

      // change from number of bytes to number of T type elements
      size_t num_ft = clover_bytes / sizeof(T);

      // Zero it all (including) (especially) the pad regions.
      // FIXME: this is not NUMA friendly necessarily
#if defined (__INTEL_COMPILER)
#pragma simd
#pragma vector nontemporal(ret_val_ft)
#endif
#pragma omp parallel for
      for(int i=0; i < num_ft; i++) {
	ret_val_ft[i] = rep<T,double>(0.0);
      }

      return ret_val;
    }

    void free(CloverBlock* p) 
    {
      BUFFER_FREE(p,clover_bytes);
    }


    int getBy() const { return By; }
    int getBz() const { return Bz; }
    int getSy() const { return Sy; }
    int getSz() const { return Sz; }
    int getPadXY() const { return PadXY; }
    int getPadXYZ() const { return PadXYZ; }
    int getPxy() const { return Pxy; }
    int getPxyz() const { return Pxyz; }
    int getNSIMT() const { return nsimt; }
    int getNumThreads() const { return num_threads; }
    int getNumCores() const { return num_cores; }
    int getMinCt() const { return MinCt; }
    int getVolCB() const { return Nxh_*Ny_*Nz_*Nt_ ; }
    CorePhase& getCorePhase(int i) { return phase[i]; }
    const CorePhase& getCorePhase(int i) const { return phase[i]; }
    int getNumPhases() const { return n_phases; }
  private:
    

    int Nxh_;
    int Nx_;
    int Ny_;
    int Nz_;
    int Nt_;

    const int Nd;
    const int By;
    const int Bz;
    const int Sy;
    const int Sz;
    const int PadXY;
    const int PadXYZ;
    int Pxy;
    int Pxyz;
    int MinCt; // Minimum no of cores in T dir
    //  MinCt = 1 for single socket/Xeon Phi
    //  MinCt = 2 for dual socket
    //  MinCt = 4 for quad socket

    const int nsimt;
    const int num_threads;
    const int num_cores;




    int nvecs_;
    int ngy_;

    // Dhiraj's new core mapping
    static const int MAX_PHASES=128;
    CorePhase phase[MAX_PHASES];
    int n_phases;
    int minCt; // Minimum no of cores in T dir

    size_t gauge_bytes;
    size_t spinor_bytes;
    size_t clover_bytes;

    //  minCt = 1 for single socket/Xeon Phi
    //  minCt = 2 for dual socket
    //  minCt = 4 for quad socket
  };
} // Namespace

#endif
