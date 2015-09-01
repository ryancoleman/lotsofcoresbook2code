// -*- C++ -*-

/*! @file
 * @brief Forward declarations for QDP
 */

namespace QDP
{

  // IO
  class TextReader;
  class TextWriter;
  class BinaryReader;
  class BinaryWriter;


  // Forward declarations
  //! dest  = random
  template<class T1, class T2>
  inline void
  fill_random(float& d, T1& seed, T2& skewed_seed, const T1& seed_mult);

  //! dest  = random
  template<class T1, class T2>
  inline void
  fill_random(double& d, T1& seed, T2& skewed_seed, const T1& seed_mult);

  //! dest  = random
  template<class T1, class T2, int N>
  inline void
  fill_random(float* d, T1& seed, T2& skewed_seed, const T1& seed_mult);

  //! dest  = random
  template<class T1, class T2, int N>
  inline void
  fill_random(double* d, T1& seed, T2& skewed_seed, const T1& seed_mult);


  namespace RNG 
  {
//  float sranf(Seed&, Seed&, const Seed&);
  }

  
  // Inner
  template<class T> class IScalar;
  template<class T, int N> class ILattice;

  // Reality
  template<class T> class RScalar;
  template<class T> class RComplex;

  // Primitives
  template<class T> class PScalar;
  template <class T, int N, template<class,int> class C> class PMatrix;
  template <class T, int N, template<class,int> class C> class PVector;
  template <class T, int N> class PColorVector;
  template <class T, int N> class PSpinVector;
  template <class T, int N> class PColorMatrix;
  template <class T, int N> class PSpinMatrix;
  template <class T> class PSeed;

  template<int N> class GammaType;
  template<int N, int m> class GammaConst;

  template<int N> class GammaTypeDP;
  template<int N, int m> class GammaConstDP;

  // Outer
  template<class T> class OScalar;
  template<class T> class OLattice;

  // Outer types narrowed to a subset
  template<class T> class OSubScalar;
  template<class T> class OSubLattice;

  // Main type
  template<class T, class C> class QDPType;

  // Expression class for QDP
  template<class T, class C> class QDPExpr;

  // Main type narrowed to a subset
  template<class T, class C> class QDPSubType;

  // Simple scalar trait class
  template<class T> struct SimpleScalar;
  template<class T> struct InternalScalar;
  template<class T> struct LatticeScalar;
  template<class T> struct PrimitiveScalar;
  template<class T> struct RealScalar;
  template<class T> struct WordType;
  template<class T> struct SinglePrecType;
  template<class T> struct DoublePrecType;

  // Empty leaf functor tag
  struct ElemLeaf
  {
    inline ElemLeaf() { }
  };

  // Empty print tag
  struct PrintTag;

  // Used for nearest neighbor shift (a map)
  class ArrayBiDirectionalMap;


} // namespace QDP

  
