// $Id: t_entry.cc,v 1.5 2004-11-22 19:31:30 edwards Exp $

/*! \file
 *  \brief Test entry/exit routines
 *
 */

#include <iostream>
#include <cstdio>

#include "qdp.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  LatticeFermion la;
  random(la);

  // Suspend communications
  // NOTE: this is not necessary for an extraction
  QDP_suspend();

  // Now can do stuff outside of QDP, e.g. call communications directly
  // but QDP is still ``turned on''
  int x=1;
  x++;

  // Resume communications
  // NOTE: this is not necessary for an extraction
  QDP_resume();

  // Must make sure all communications are allowed
  XMLFileWriter toxml("t_entry.xml");
  push(toxml,"t_entry");

  // Pull out the data
//  Subset& even = rb[0];

  // WARNING: the variable "sa" will be node dependent, e.g. on each node
  // there will be different values. So, the write statement below of 
  // sa will only print that from the **primary** node. However, that
  // is the point of extract/insert - it allows the user to manipulate their
  // data into a QDP field on each node, from their own private data format.
  multi1d<Fermion> sa(Layout::sitesOnNode());
  sa = zero;
  QDP_extract(sa, la, even);


  // Fiddle with the innards of this new field sa . This is for demonstration.
  // Pull out some goodies
  int color_index = 0, spin_index = 0;
  ColorVector colorVec = peekSpin(sa[0],spin_index);   // use a temporary - not needed
  Complex sitecomp = peekColor(colorVec,color_index);
  Real re = real(sitecomp) + 1.0;
  Real im = imag(sitecomp) + 2.0;
  sitecomp = cmplx(re,im);

  ColorVector sitecolor = zero;

  // Shove the modified goodies back into the local temp variable
  pokeSpin(sa[0],
	   pokeColor(sitecolor,sitecomp,color_index),
	   spin_index);


  // Stuff sa back into the QDP field. Zero out QDP field for testing
  la = zero;
  QDP_insert(la, sa, even);

  push(toxml, "Site_field");
  write(toxml,"sa", sa);
  pop(toxml);

  push(toxml, "Lattice_field");
  write(toxml,"la", la);
  pop(toxml);

  pop(toxml);
  toxml.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
