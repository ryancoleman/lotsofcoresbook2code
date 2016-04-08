#ifndef QPHIX_WILSON_H
#define QPHIX_WILSON_H

#include "qphix/linearOp.h"
#include "qphix/dslash_def.h"

namespace QPhiX
{ 

  template<typename FT, int veclen, int soalen, bool compress12>
    class EvenOddWilsonOperator : public EvenOddLinearOperator<FT,veclen,soalen,compress12> {
  public: 
    typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock FourSpinorBlock;
    typedef typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock SU3MatrixBlock;

    EvenOddWilsonOperator(const double Mass_,
			  SU3MatrixBlock* u_[2],
			  Geometry<FT,veclen,soalen,compress12>* geom_,
			  double t_boundary_, 
			  double aniso_fac_s_,
			  double aniso_fac_t_): Mass(Mass_), D(new Dslash<FT, veclen,soalen,compress12>(geom_, t_boundary_, aniso_fac_s_, aniso_fac_t_))
      {

	Geometry<FT,veclen, soalen, compress12>& geom = D->getGeometry();
	u[0] = u_[0];
	u[1] = u_[1];
	tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
	
	mass_factor_alpha = (double)4 + Mass;
	mass_factor_beta = (double)(0.25)/ mass_factor_alpha; 
    }


    ~EvenOddWilsonOperator() {
      Geometry<FT,veclen, soalen, compress12>& geom = D->getGeometry();
      geom.free(tmp);
      delete D;
    }


    void operator()(FourSpinorBlock *res, const FourSpinorBlock* in, int isign) {
    
      D->dslash(tmp, in, u[1], isign, 1);
      D->dslashAChiMinusBDPsi(res, tmp, in, u[0], mass_factor_alpha, mass_factor_beta, isign, 0);
    }

    Geometry<FT,veclen, soalen, compress12>& getGeometry() { return D->getGeometry(); }
  private:
    double Mass;
    Dslash<FT, veclen,soalen, compress12>* D;
    const SU3MatrixBlock *u[2];
    FourSpinorBlock *tmp;
    double mass_factor_alpha; 
    double mass_factor_beta;
  }; // Class
}; // Namespace

    

#endif
