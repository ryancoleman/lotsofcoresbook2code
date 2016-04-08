// $Id: taproj.cc,v 1.3 2004-08-10 14:38:20 edwards Exp $
// TAPROJ

/*! \file
 *  \brief Take the traceless antihermitian projection of a color matrix
 */

#include "qdp.h"
#include "examples.h"

using namespace QDP;

//! Take the traceless antihermitian projection of a color matrix
/*!
 * \ingroup gauge
 *
 *  a = (1/2)[a - a_dag] - Tr[(1/2)*(a - a_dag)]/Nc
 *
 * that is the anti-hermitian traceless part of a 
 *
 * Arguments:
 *
 *  \param a        LatticeColorMatrix          (Modify)
 */

void taproj(LatticeColorMatrix& a)
{
  START_CODE();

  // a = a - a_dagger  --- a -> antihermitian matrix
  LatticeColorMatrix aux_1 = a;
  a -= adj(aux_1);
 
  if (Nc > 1) {
    // tmp = Im Tr[ a ]
    LatticeReal tmp = imag(trace(a));
    
    // a = a - (1/Nc) * Im Tr[ a] = a - (1/Nc)*tmp
    tmp *= (Real(1)/Real(Nc));

    // Is this a fill or a UnitMatrix*I?
    LatticeColorMatrix aux = cmplx(0, tmp);
    a -= aux;
  }

  // Normalisation to make taproj idempotent
  a *= (Real(1)/Real(2));
  
  END_CODE();
}
