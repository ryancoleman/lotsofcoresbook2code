#ifndef QPHIX_LINEAR_OPERATOR_H
#define QPHIX_LINEAR_OPERATOR_H


#include <qphix/geometry.h>

namespace QPhiX { 

  template<typename FT, int veclen, int soalen, bool compress> 
  class EvenOddLinearOperator {
  public:
    typedef typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock FourSpinorBlock;
    virtual void operator()(FourSpinorBlock *res, const FourSpinorBlock* in, int isign) = 0;

    virtual Geometry<FT,veclen,soalen,compress>& getGeometry(void)=0;

  };






}



#endif
