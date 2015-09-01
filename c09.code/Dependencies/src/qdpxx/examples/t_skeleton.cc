// $Id: t_skeleton.cc,v 1.5 2003-06-07 19:09:32 edwards Exp $
/*! \file
 *  \brief Skeleton of a QDP main program
 */

#include "qdp.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  for(int i=0; i < Nd; ++i)
    nrow[i] = 2;         // Set the lattice size to 2^4

  // Insert the lattice size into the Layout
  // There can be additional calls here like setting the number of
  // processors to use in an SMP implementation
  Layout::setLattSize(nrow);

  // Create the lattice layout
  // Afterwards, QDP is useable
  Layout::create();

  // Do some wonderful and amazing things - impress your friends


  // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
