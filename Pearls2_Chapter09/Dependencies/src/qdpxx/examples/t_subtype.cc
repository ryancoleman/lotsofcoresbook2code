// -*- C++ -*-
//
// $Id: t_foo.cc,v 1.44 2007-03-15 03:14:55 edwards Exp $
//
/*! \file
 *  \brief Silly little internal test code
 */


#include "qdp.h"

using namespace QDP;

typedef OSubLattice< PScalar< PColorVector< RComplex<REAL>, Nc> > > SubLatticeColorVector;


class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc() {}
  int operator() (const multi1d<int>& coordinate) const ;
  int numSubsets() const ;
};


int TimeSliceFunc::operator() (const multi1d<int>& coordinate) const
{
  return coordinate[Nd-1];
}

int TimeSliceFunc::numSubsets() const
{
  return Layout::lattSize()[Nd-1];
}



void test_tslice()
{
  Set tslice;
  tslice.make(TimeSliceFunc());

  LatticeColorVector a;
  random(a);

  const int Nt = Layout::lattSize()[Nd-1];

  // Default create a set of Nt SLCV. The default constructor does not allocate memory
  multi1d<SubLatticeColorVector> vec(Nt);

  // Now attach a subset to each of the SLCV.
  // Then copy the according timeslice of 'a' into each SLCV
  for (int i=0 ; i < Nt ; ++i ) {
    vec[i].setSubset( tslice[i] );
    zero_rep(vec[i]);   // does something only on nodes where the subset is non-empty 
    vec[i] = a;         // dito
  }

  // Now, do some contractions.
  multi1d<Complex> contr( Nt );
  for (int i=0 ; i < Nt ; ++i ) {
    contr[i] = innerProduct( a , vec[i] );
  }

  // Another way to construct SLCV:
  SubLatticeColorVector tmp( tslice[Nt-1] , a );  // This attaches the subset to 'tmp' and copies the according TS from 'a'.
  zero_rep(a);
  a = tmp; // This writes into the subset that was attached to tmp

  // The OSubLattice is *not* fleshed into PETE. Thus only those operations that are explicitly implemented work, i.e. like innerProduct
  // This does not work:
  // vec[0] = vec[1] + vec[2];

  if (0)
  {
    XMLFileWriter foo("out.xml");
    push(foo,"foo");
    write(foo,"contr",contr);
    write(foo,"tmp",a);
    foo.close();
  }



  
}




void stuff()
{
  LatticeColorVector a;
  LatticeColorVector b;
  LatticeColorMatrix lcm;

#if 1
  zero_rep(a);
  random(b);

  a[rb[0]] = b;

  ColorVector cs,cs2;
  cs[rb[0]] = cs2;


  LatticeComplex z;
  pokeColor(lcm[rb[0]],z,2,1);
#endif


  SubLatticeColorVector s( rb[1] , a );
  s = b;

  LatticeColorVector wr;
  zero_rep(wr);
  wr[rb[1]] = b;

  LatticeColorVector b2;
  zero_rep(b2);
  b2[rb[1]] = b;

  SubLatticeColorVector s2(s);
  s2 = s;
  a[rb[0]] = s;
  s = a[rb[1]];


#if 1
  multi1d<SubLatticeColorVector> vec(10);

  for(int i=0;i<10;i++)
    vec[i].setSubset(rb[i&1]);

  for(int i=0;i<10;i++)
    vec[i] = a;

  vec[0].setSubset( rb[1] );
#endif

  Complex aa = innerProduct(s,a);
  aa = innerProduct(a,s);
  aa = innerProduct(s,s);
}




int main(int argc, char *argv[])
{
  QDP_initialize(&argc, &argv);

  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements

  Layout::setLattSize(nrow);
  Layout::create();



  test_tslice();




  QDP_finalize();
  //exit(0);
}
