// -*- C++ -*-

/*! \file
 * \brief Yet another random number generator
 *
 * This is an RNG independent from the vanilla QDP RNG.
 *
 * Seeds are four 12 bit integers
 */

#ifndef QDP_RANNYU_H
#define QDP_RANNYU_H

#include "qdp.h"

namespace QDP
{
  namespace RANNYU
  {
    //! The RNG. Has side effects.
    double random();

    //! Seed has been set by default - this allows one to override it
    void setrn(const multi1d<int>& iseed);

    //! Recover the seed
    multi1d<int> savern();

    //! Hold state and RNG.
    struct RNGState_t
    {
      multi1d<int> seed;   /*!< The state of the RNG. */
      double       ran;    /*!< Output of the RNG. */
    };

    //! The RNG, but packaging all the state with it. No side effects.
    void random(RNGState_t& ran_state);

  } // namespace RANNYU

} // namespace QDP

#endif
