// $Id: nersc2ildg.cc,v 1.3 2006-01-12 02:17:39 bjoo Exp $
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
  std::string NERSC_file_name;
  std::string ILDG_file_name;
  std::string dataLFN;
  int output_size;
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
    XMLReader param("./nersc2ildg.xml");
    

    XMLReader paramtop(param, "/nersc2ildg");
    read(paramtop, "NERSC_file", p.NERSC_file_name);
    read(paramtop, "ILDG_file", p.ILDG_file_name);
    read(paramtop, "dataLFN", p.dataLFN);
    read(paramtop, "nrow", nrow);
    read(paramtop, "output_size", p.output_size);

  } catch(const std::string& e) { 
    QDPIO::cout << "Caught exception while reading XML: " << e << endl;
    QDP_abort(1);
  }

  if(! ( p.output_size == 32 || p.output_size == 64 )  ) {
    QDPIO::cerr << "Output Size must be either 32 or 64 (bits)" << endl;
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
  readArchiv(u, p.NERSC_file_name);

  // Do checks here.
  XMLBufferWriter file_metadata;
  push(file_metadata, "file_metadata");
  write(file_metadata, "annotation", "NERSC Config Converted by QDP++ NERS2ILDG");
  pop(file_metadata);

  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  QDPFileWriter ildg_out(file_metadata,  
			 p.ILDG_file_name,
			 QDPIO_SINGLEFILE,
			 QDPIO_SERIAL,
  			 p.dataLFN);

  XMLBufferWriter record_metadata;
  push(record_metadata, "record_metadata");
  write(record_metadata, "annotation", "NERSC Config Record Converted by QDP++ NERSC2ILDG");
  pop(record_metadata);

  switch( p.output_size ) { 
  case 32:
    {
      multi1d<LatticeColorMatrixF> u_single_out(Nd);
      for(int mu=0; mu < Nd; mu++) {
	u_single_out[mu] = u[mu];
      }
      ildg_out.write(record_metadata, u_single_out);
    }
    break;
  case 64:
    {
      multi1d<LatticeColorMatrixD> u_double_out(Nd);
      for(int mu=0; mu < Nd; mu++) {
	u_double_out[mu] = u[mu];
      }
      ildg_out.write(record_metadata, u_double_out);
    }
    break;
  default: 
    {
      QDPIO::cerr << "Output precision must be either 32 or 64. You entered " << p.output_size << endl;
      QDP_abort(1);
      break;
    }
  }
  ildg_out.close();

  // Reread the ILDG File
  XMLReader record_in_xml;
  XMLReader file_in_xml;
  QDPFileReader ildg_back_in(file_in_xml, p.ILDG_file_name, QDPIO_SERIAL);

  multi1d<LatticeColorMatrix> u_back_in(Nd);
  ildg_back_in.read(record_in_xml, u_back_in);

  record_in_xml.print(cout);
  cout.flush();

  MesPlq(u_back_in, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "Read Back Plaquette " << w_plaq << endl;

    // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
