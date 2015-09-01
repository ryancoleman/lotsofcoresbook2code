// $Id: t_nersc.cc,v 1.8 2005-08-22 21:17:27 edwards Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "qdp_iogauge.h"
#include "examples.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  std::string filename = "t_nersc.cfg";
//  string filename = "tmp.cfg";
//  string filename = "u_TEST_4223.101";

  XMLFileWriter xml("t_nersc.xml");
  push(xml, "t_nersc");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

#if 1
  {
    multi1d<LatticeColorMatrix> u(Nd);
    Double w_plaq, s_plaq, t_plaq, link;

    QDPIO::cout << "Start gaussian\n";
    for(int m=0; m < u.size(); ++m)
      gaussian(u[m]);

    // Reunitarize the gauge field
    QDPIO::cout << "Start reunit\n";
    for(int m=0; m < u.size(); ++m)
      reunit(u[m]);

    QDPIO::cout << "Start mesplq\n";
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
    QDPIO::cout << "link = " << link << std::endl;

    // Write out the results
    push(xml,"Initial_observables");
    write(xml,"w_plaq", w_plaq);
    write(xml,"link", link);
    pop(xml);

    // Now write the gauge field in NERSC format
    QDPIO::cout << "Trying to write NERSC Archive  " << filename << std::endl;
    writeArchiv(u, filename);
  }
#endif

  {
    multi1d<LatticeColorMatrix> u(Nd);
    Double w_plaq, s_plaq, t_plaq, link;

    QDPIO::cout << "Trying to read back config" << std::endl;

    XMLReader gauge_xml;
    readArchiv(gauge_xml, u, filename);
 
    QDPIO::cout << "Dump the gauge xml" << std::endl;
    push(xml, "Gauge_xml");
    xml << gauge_xml;
    pop(xml);

    // Try out the plaquette routine
    QDPIO::cout << "Start mesplq\n";
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
    QDPIO::cout << "link = " << link << std::endl;

    // Write out the results
    push(xml,"Final_observables");
    write(xml,"w_plaq", w_plaq);
    write(xml,"link", link);
    pop(xml);
  }

  pop(xml);
  xml.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
