#ifndef QPHIX_ABS_SOLVER_H
#define QPHIX_ABS_SOLVER_H

#include "qphix/geometry.h"

namespace QPhiX {

  template<typename FT, int V, int S, bool compress12>
  class AbstractSolver {
  public:
    typedef typename Geometry<FT,V,S,compress12>::FourSpinorBlock Spinor;
    virtual void operator()(Spinor* x,
			    const Spinor *rhs,
			    const double RsdTarget, 
			    int& niters, 
			    double& rsd_sq_final,
			    unsigned long& site_flops,
			    unsigned long& mv_apps,
			    int isign,
			    bool verboseP) = 0;

    virtual void tune() = 0;
    virtual  Geometry<FT,V,S,compress12>& getGeometry() =  0;
  };



}




#endif 
