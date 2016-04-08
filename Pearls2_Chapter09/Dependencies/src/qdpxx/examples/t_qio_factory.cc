// $Id: t_qio_factory.cc,v 1.4 2004-03-29 21:29:31 edwards Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"

using namespace QDP;

typedef PSpinMatrix< PColorMatrix < RComplex<REAL>, Nc>, Ns > MyProp;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int dims[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = dims;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  LatticePropagator quark_propagator = zero;

  int site_size;
  site_size = sizeof(MyProp);

  QDPIO::cout << "Site Data Size = " << site_size << endl << flush ;

  // Get a buffer for a site-s worth of data 
  char *buf = new char [ site_size ];
  if( buf == 0x0 ) { 
    QDPIO::cerr << "Unable to allocate buf" << endl;
    QDP_abort(1);
  }

  // Fill the first 128 members of the lattice propagator 
  int i;
  MyProp* site;
  for(i = 0; i < 128; i++) {
    site = &quark_propagator.elem(i);
    memset((void *)site, i, site_size);
  }

  // Now try to get the last 128 sites with the Factory Get function 
  int coords[4] = {0, 0, 0, 0};
  int j;

  int x, y, z, t;
  for(x = 0; x < dims[0]/2; x++) { 
    coords[0] = 2*x;
    
    QDPFactoryGet<LatticePropagator>(buf, coords, &quark_propagator);
    for(i = 0, j = 0; i < site_size; i++) { 
      cout << (int)buf[i] << " ";
      j++;
      if ( j > 30) { 
	cout << endl << flush;
	j = 0;
      }
    }
    cout << endl;
  }
  
  // Time to bolt
  QDP_finalize();

  exit(0);
}
