// $Id: t_cugauge.cc,v 1.1 2005-02-21 15:27:05 bjoo Exp $
/*! \file
 *  \brief Skeleton of a QDP main program
 */

#include "qdp.h"
#include "qdp_iogauge.h"
#include "examples.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  nrow[0]=16;
  nrow[1]=16;
  nrow[2]=16;
  nrow[3]=32;

  // Insert the lattice size into the Layout
  // There can be additional calls here like setting the number of
  // processors to use in an SMP implementation
  Layout::setLattSize(nrow);

  // Create the lattice layout
  // Afterwards, QDP is useable
  Layout::create();
  XMLFileWriter xml("t_cugauge.xml");

  // Do some wonderful and amazing things - impress your friends
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;
  QDPIO::cout << "Trying to read config" << std::endl;
  XMLReader gauge_xml;
  readArchiv(gauge_xml, u, "test_lattice");

  QDPIO::cout << "Dump the gauge xml" << std::endl;
  xml << gauge_xml;
                                                                                
  // Try out the plaquette routine
  QDPIO::cout << "Start mesplq\n";
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
  QDPIO::cout << "link = " << link << std::endl;

  // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
