#ifndef QPHIX_DSLASH_DEF_H
#define QPHIX_DSLASH_DEF_H

#include "qphix/dslash_utils.h"
#include "qphix/geometry.h"
#include "qphix/comm.h"

namespace QPhiX
{

  /* Specialize - Dslash of float */
  template<typename FT, int veclen, int soalen, bool compress12>  /* The teplate is the floating point type */
  class Dslash {
  public:
    // Pack gauges... 
    typedef typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  SU3MatrixBlock;
    typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock FourSpinorBlock;
    typedef typename Geometry<FT,veclen,soalen,compress12>::TwoSpinorBlock TwoSpinorBlock;
    

    Dslash(Geometry<FT,veclen, soalen,compress12>* geom_,
	   double t_boundary_,
	   double aniso_coeff_S_,
	   double aniso_coeff_T_);

    /* Destructor */
    ~Dslash();
    
    /* Apply the operator */
    void dslash(FourSpinorBlock* res, 
		const FourSpinorBlock* psi, 
		const SU3MatrixBlock* u,
		int isign,
		int cb);
    
    void dslashAChiMinusBDPsi (FourSpinorBlock* res, 
			       const FourSpinorBlock* psi, 
			       const FourSpinorBlock* chi,
			       const SU3MatrixBlock* u,
			       double alpha, 
			       double beta,
			       int isign,
			       int cb);



    Geometry<FT,veclen,soalen,compress12>& getGeometry(void) { return (*s); }
  private:
    Geometry<FT,veclen,soalen,compress12>* s;
    Comms<FT,veclen,soalen,compress12>* comms;

    int log2veclen;
    int log2soalen;
    void init();


    const  int n_threads_per_core;

    const int By;
    const int Bz;
    const int NCores;
    const int Sy;
    const int Sz;
    const int PadXY;
    const int PadXYZ;
    const int MinCt;

    // Antiperiodic Boundary
    const double t_boundary;
    bool amIPtMin;
    bool amIPtMax;

    // Anisotropy
    double aniso_coeff_S;
    double aniso_coeff_T;


    Barrier* gBar;
    Barrier*** barriers;

    BlockPhase* block_info;
    int *tmpspc_all; // Space for YZ offsets etc
    
    // Masks for comms in X & Y
    unsigned int xbmask_x0_xodd[2];
    unsigned int xfmask_xn_xodd[2];
    unsigned int ybmask_y0;
    unsigned int yfmask_yn;
    
    // Hide Free Constructor 
    Dslash();
    void DyzPlus(int tid,
		 const FourSpinorBlock* psi, 
		 FourSpinorBlock* res,
		 const SU3MatrixBlock* u, 
		 int cb);

    void DyzMinus(int tid,
		  const FourSpinorBlock* psi, 
		  FourSpinorBlock* res,
		  const SU3MatrixBlock* u, 
		  int cb);


    void DyzPlusAChiMinusBDPsi(int tid,
			       const FourSpinorBlock* psi, 
			       const FourSpinorBlock* chi,
			       FourSpinorBlock* res,
			       const SU3MatrixBlock* u, 
			       double alpha, 
			       double beta, 
			       int cb);


    void DyzMinusAChiMinusBDPsi(int tid,
				const FourSpinorBlock* psi, 
				const FourSpinorBlock* chi,
				FourSpinorBlock* res,
				const SU3MatrixBlock* u, 
				double alpha, 
				double  beta,
				int cb);




    void DPsiPlus(const SU3MatrixBlock *u, const FourSpinorBlock *psi_in, FourSpinorBlock *res_out, int cb);


    void DPsiMinus(const SU3MatrixBlock *u, const FourSpinorBlock *psi_in, FourSpinorBlock *res_out, int cb);



    void DPsiPlusAChiMinusBDPsi(const SU3MatrixBlock *u, 
				const FourSpinorBlock *psi_in, 
				const FourSpinorBlock *chi, 
				FourSpinorBlock *res_out, 
				double alpha, double beta, 
				int cb);


    void DPsiMinusAChiMinusBDPsi(const SU3MatrixBlock *u, 
				 const FourSpinorBlock *psi_in, 
				 const FourSpinorBlock *chi, 
				 FourSpinorBlock *re_out, 
				 double alpha, double beta, 
				 int cb);



#ifdef QPHIX_DO_COMMS
#warning defining face packers
     // PACK FACE FOR SENDING
    void packFaceDir(int tid,
		     const FourSpinorBlock *psi,
		     FT *res,
		     int cb, int dir, int fb, int isPlus);

    //  RECEIVE AND COMPLETE FACE 
    void completeFaceDir(int tid,
			 const FT* psi,
			 FourSpinorBlock* res,
			 const SU3MatrixBlock* u,
			 const double beta, 
			 int cb, int dir, int fb, int isPlus);
#endif



  }; // Class



} // Namespace
#include "qphix/dslash_body.h"

#ifdef QPHIX_DO_COMMS
#include "qphix/face.h"
#endif
#endif
