// $Id: test_ildglat.cc,v 1.1 2006-01-13 16:41:10 bjoo Exp $
/*! \file
 *  \brief Skeleton of a QDP main program
 */

#include "qdp.h"
#include "qdp_iogauge.h"
#include <iostream>
#include <string>
#include "examples.h"


using namespace QDP;

typedef struct { 
  std::string ILDG_file_name;
} UserInput;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  UserInput p;

  try { 
    XMLReader param("./test_ildglat.xml");
    

    XMLReader paramtop(param, "/ildglat");
    read(paramtop, "ILDG_file", p.ILDG_file_name);
    read(paramtop, "nrow", nrow);

  } catch(const std::string& e) { 
    QDPIO::cout << "Caught exception while reading XML: " << e << endl;
    QDP_abort(1);
  }

  // Insert the lattice size into the Layout
  // There can be additional calls here like setting the number of
  // processors to use in an SMP implementation
  Layout::setLattSize(nrow);
  
  // Create the lattice layout
  // Afterwards, QDP is useable
  Layout::create();

  // Try to read the NERSC Archive file
  multi1d<LatticeColorMatrix> u(Nd);

  XMLReader record_in_xml;
  XMLReader file_in_xml;


  // Set QIO Verbosity to DEBUG
  QIO_verbose(QIO_VERB_DEBUG);



  QDPFileReader ildg_reader(file_in_xml, p.ILDG_file_name, QDPIO_SERIAL);
 
  ildg_reader.read(record_in_xml, u);
  Double w_plaq, s_plaq, t_plaq, link;

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq " << w_plaq << endl;
  QDPIO::cout << "link " << link << endl;

    // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
