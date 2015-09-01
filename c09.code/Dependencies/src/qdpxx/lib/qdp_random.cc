//
// Random number generator support


#include "qdp.h"

namespace QDP {

// Random number generator namespace
/* 
 * A collection of routines and data for supporting random numbers
 * 
 * It is a linear congruential with modulus m = 2**47, increment c = 0,
 * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
 */

namespace RNG
{
  //! Global (current) seed
  Seed ran_seed;
  //! RNG multiplier
  Seed ran_mult;
  //! RNG multiplier raised to the volume+1
  Seed ran_mult_n;
  //! The lattice of skewed RNG multipliers
  LatticeSeed *lattice_ran_mult;

    //! Find the number of bits required to represent x.
  int numbits(int x)
  {
    int num = 1;
    int iceiling = 2;
    while (iceiling <= x)
    {
      num++;
      iceiling *= 2;
    }

    return num;
  }


  //! Initialize the random number generator with a default seed
  void initDefaultRNG()
  {
    RNG::initRNG();

    Seed seed = 11;
    RNG::setrn(seed);
  }


  //! Initialize the internals of the random number generator
  void initRNG()
  {
    int old_profile_level = setProfileLevel(0);

    /* Multiplier used. Use big integer arithmetic */
    Seed seed_tmp3;
    Seed seed_tmp2;
    Seed seed_tmp1;
    Seed seed_tmp0;

    seed_tmp3 = 1222;
    seed_tmp2 = (seed_tmp3 << 12) | 1498;
    seed_tmp1 = (seed_tmp2 << 12) | 712;
    seed_tmp0 = (seed_tmp1 << 12) | 1645;

    ran_mult = seed_tmp0;

    // Find the number of bits it takes to represent the total lattice volume.
    // NOTE: there are no lattice size restrictions here.
    int nbits = numbits(Layout::vol());

    /* Get the lattice coordinate of each site (note the origin is 0) and
     *   build up a lexicographic ordering for the lattice. The definition
     *   here is totally arbitrary and only this routine needs to worry
     *   about it.  The lexicographic value of site K is
     *
     *     lexoc(k) = sum_{i = 1, ndim} x(k,i)*L^i     +   1
     */
    LatticeInteger lexoc;
    lexoc = Layout::latticeCoordinate(Nd-1);

    for(int m=Nd-2; m>=0; --m)
    {
      lexoc *= Layout::lattSize()[m];
      lexoc += Layout::latticeCoordinate(m);
    }

    lexoc += 1;

    /*
     * Setup single multiplier ( a^1 ) on each site 
     */
    LatticeSeed laa;
    laa = ran_mult;

    /*
     * Calculate the multiplier  a^n  where n = lexicographic numbering of the site.
     *   Put one into each an_f, then multiply them by  a^(2^i) under the context
     *   flag where i is the bit number of the lexicographic numbering.
     *   In other words, the very first site is multiplied by a.
     */
    LatticeSeed lattice_ran_mult_tmp;
    lattice_ran_mult_tmp = 1;

    LatticeSeed laamult;
    LatticeBoolean lbit;

    for(int i=0; i<nbits; ++i)
    {
      lbit = (lexoc & 1) > 0;

      laamult = lattice_ran_mult_tmp * laa;
      copymask(lattice_ran_mult_tmp,lbit,laamult);

      lexoc >>= 1;
      laamult = laa * laa;
      laa = laamult;
    }

    // Calculate separately the multiplier for the highest lexicographically ordered site.
    // NOTE: I'm changing the meaning here slightly, but in an important way.
    // Technically, ran_mult_n = ran_mult^{vol} . Instead, I'm going to throw
    // away an rng call after every lattice call. So, I will DEFINE
    //
    //   ran_mult_n = ran_mult^{vol + 1}
    //
    bool bit;
    Seed aa;
    Seed aamult;

    int ibit = Layout::vol();
    aa = ran_mult;
//    ran_mult_n = 1;    // produces def    ran_mult_n = ran_mult^{vol}
    ran_mult_n = ran_mult;   // produces def  ran_mult_n = ran_mult^{vol+1}

    for(int i=0; i<nbits; ++i)
    {
      bit = (ibit & 1) > 0;

      aamult = ran_mult_n * aa;
      if (bit)
	ran_mult_n = aamult;

      ibit >>= 1;
      aamult = aa * aa;
      aa = aamult;
    }

    lattice_ran_mult = new LatticeSeed;
    if( lattice_ran_mult == 0x0 ) { 
      QDP_error_exit("Unable to allocate ran_mult\n");
    }

    *lattice_ran_mult = lattice_ran_mult_tmp;
    QDPIO::cout << "Finished init of RNG" << std::endl; 

    setProfileLevel(old_profile_level);
  }


  //! Initialize the random number generator seed
  void setrn(const Seed& seed)
  {
    ran_seed = seed;
  }


  //! Return a copy of the random number seed
  void savern(Seed& seed)
  {
    seed = ran_seed;
  }


  //! Scalar random number generator. Done on the front end. */
  /*! 
   * It is linear congruential with modulus m = 2**47, increment c = 0,
   * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
   */
  float sranf(Seed& seed, Seed& skewed_seed, const Seed& seed_mult)
  {
    /* Calculate the random number and update the seed according to the
     * following algorithm
     *
     * FILL(twom11,TWOM11);
     * FILL(twom12,TWOM12);
     * i3 = ran_seed(3)*ran_mult(0) + ran_seed(2)*ran_mult(1)
     *    + ran_seed(1)*ran_mult(2) + ran_seed(0)*ran_mult(3);
     * i2 = ran_seed(2)*ran_mult(0) + ran_seed(1)*ran_mult(1)
     *    + ran_seed(0)*ran_mult(2);
     * i1 = ran_seed(1)*ran_mult(0) + ran_seed(0)*ran_mult(1);
     * i0 = ran_seed(0)*ran_mult(0);
     *
     * ran_seed(0) = mod(i0, 4096);
     * i1          = i1 + i0/4096;
     * ran_seed(1) = mod(i1, 4096);
     * i2          = i2 + i1/4096;
     * ran_seed(2) = mod(i2, 4096);
     * ran_seed(3) = mod(i3 + i2/4096, 2048);
     *
     * sranf = twom11*(TO_REAL32(VALUE(ran_seed(3)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(2)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(1)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(0)))))));
     */
    Real _sranf;
    float _ssranf;
    Seed ran_tmp;

    _sranf = seedToFloat(skewed_seed);
    cast_rep(_ssranf, _sranf);

    ran_tmp = seed * seed_mult;
    seed = ran_tmp;

    ran_tmp = skewed_seed * seed_mult;
    skewed_seed = ran_tmp;

    return _ssranf;
  }


  //! Scalar random number generator. Done on the front end. */
  /*! 
   * It is linear congruential with modulus m = 2**47, increment c = 0,
   * and multiplier a = (2**36)*m3 + (2**24)*m2 + (2**12)*m1 + m0.  
   */
  void sranf(float* d, int N, Seed& seed, ILatticeSeed& skewed_seed, const Seed& seed_mult)
  {
    /* Calculate the random number and update the seed according to the
     * following algorithm
     *
     * FILL(twom11,TWOM11);
     * FILL(twom12,TWOM12);
     * i3 = ran_seed(3)*ran_mult(0) + ran_seed(2)*ran_mult(1)
     *    + ran_seed(1)*ran_mult(2) + ran_seed(0)*ran_mult(3);
     * i2 = ran_seed(2)*ran_mult(0) + ran_seed(1)*ran_mult(1)
     *    + ran_seed(0)*ran_mult(2);
     * i1 = ran_seed(1)*ran_mult(0) + ran_seed(0)*ran_mult(1);
     * i0 = ran_seed(0)*ran_mult(0);
     *
     * ran_seed(0) = mod(i0, 4096);
     * i1          = i1 + i0/4096;
     * ran_seed(1) = mod(i1, 4096);
     * i2          = i2 + i1/4096;
     * ran_seed(2) = mod(i2, 4096);
     * ran_seed(3) = mod(i3 + i2/4096, 2048);
     *
     * sranf = twom11*(TO_REAL32(VALUE(ran_seed(3)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(2)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(1)))
     *       + twom12*(TO_REAL32(VALUE(ran_seed(0)))))));
     */
    ILatticeReal _sranf;
    Seed ran_tmp1;
    ILatticeSeed ran_tmp2;

    _sranf = seedToFloat(skewed_seed);
    for(int i=0; i < N; ++i)
    {
      cast_rep(d[i], getSite(_sranf,i));
    }

    ran_tmp1 = seed * seed_mult;
    seed = ran_tmp1;

    ran_tmp2 = skewed_seed * seed_mult;
    skewed_seed = ran_tmp2;
  }

};

} // namespace QDP;
