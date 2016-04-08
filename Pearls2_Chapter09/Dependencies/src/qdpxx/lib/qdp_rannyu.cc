//! Yet another random number generator
/*!
 *      It is a linear congruential with modulus m = 2**48, increment c = 1,
 *      and multiplier a = (2**36)*m1 + (2**24)*m2 + (2**12)*m3 + m4. 
 *      The multiplier is stored in common (see subroutine setrn)
 *      and is set to a = 31167285 (recommended by Knuth, vol. 2,
 *      2nd ed., p. 102).
 *
 *     Multiplier is 31167285 = (2**24) + 3513*(2**12) + 821.
 *        Recommended by Knuth, vol. 2, 2nd ed., p. 102.
 *     (Generator is linear congruential with odd increment
 *        and maximal period, so seed is unrestricted: it can be
 *        either even or odd.)
 */

#include "qdp_rannyu.h"
#include "qdp.h"

namespace QDP
{
  namespace RANNYU
  {
    namespace
    {
      int m[4] = {0, 1, 3513, 821};

      multi1d<int> ran_seed;
      int __default_seed[4] = {13, 15, 1, 17};
      bool inited = false;
      
      double twom12 = 1/4096.0;
    }

    // Hide the actual RNG
    namespace
    {
      void __rand(double& ran, multi1d<int>& ll)
      {
	int ii[4];

	ii[0] = ll[0]*m[3] + ll[1]*m[2] + ll[2]*m[1] + ll[3]*m[0];
	ii[1] = ll[1]*m[3] + ll[2]*m[2] + ll[3]*m[1];
	ii[2] = ll[2]*m[3] + ll[3]*m[2];
	ii[3] = ll[3]*m[3] + 1;
	ll[3] = ii[3] & 4095;
	ii[2] = ii[2] + (ii[3] >> 12);
	ll[2] = ii[2] & 4095;
	ii[1] = ii[1] + (ii[2] >> 12);
	ll[1] = ii[1] & 4095;
	ll[0] = (ii[0] + (ii[1] >> 12)) >> 12;
	ran = twom12*((double)ll[0] + twom12*((double)ll[1] + twom12*((double)ll[2] + twom12*((double)ll[3]))));
      }
    }


    //! The RNG. Has side effects.
    double random()
    {
      double ran;

      // Initialize the global seed if needed
      if (! inited)
      {
	ran_seed.resize(4);
	for(int i=0; i < 4; ++i)
	  ran_seed[i] = __default_seed[i];

	inited = true;
      }

      __rand(ran, ran_seed);
      return ran;
    }


    //! Seed has been set by default - this allows one to override it
    void setrn(const multi1d<int>& iseed)
    {
      if (iseed.size() != 4)
      {
	QDPIO::cerr << __func__ << ": rannyu seed is not length 4\n";
	QDP_abort(1);
      }
	
      ran_seed = iseed;
    }


    //! Recover the seed
    multi1d<int> savern()
    {
      return ran_seed;
    }


    // The RNG, but packaging all the state with it. No side effects.
    void random(RNGState_t& ran_state)
    {
      if (ran_state.seed.size() != 4)
      {
	QDPIO::cerr << __func__ << ": rannyu seed is not length 4\n";
	QDP_abort(1);
      }
	
      __rand(ran_state.ran, ran_state.seed);
    }

  } // namespace RANNYU

} // namespace QDP

