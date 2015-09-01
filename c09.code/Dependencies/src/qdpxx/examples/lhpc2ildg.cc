// $Id: lhpc2ildg.cc,v 1.1 2007-12-05 19:44:20 bjoo Exp $
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
  std::string dataLFN;
  int output_size;
} UserInput;

int main(int argc, char *argv[])
{
  // Put the machine into a known state

  std::string input_xml="";
  std::string output_xml="";

  for(int i=0; i < argc; i++) { 
    if ( strncmp(argv[i], "-i", 2) == 0 ) {
      input_xml = argv[++i];
    }
    if ( strncmp(argv[i], "-o", 2) == 0 ) {
      output_xml = argv[++i];
    }
  }

  if( input_xml == "" ) { 
    input_xml = "DATA";
  }
  if( output_xml == "" ) { 
    output_xml = "XMLDAT";
  }

  QDP_initialize(&argc, &argv);
  QDPIO::cout << "Input file =" << input_xml << endl;
  QDPIO::cout << "Output file=" << output_xml << endl;


  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  UserInput p;

  try { 
    XMLReader param(input_xml);
    

    XMLReader paramtop(param, "/lhpc2ildg");
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
  
  XMLReader file_in_xml;
  XMLReader record_in_xml;

  QDPFileReader lhpc_in(file_in_xml,	
			p.ILDG_file_name,
			QDPIO_SERIAL);

  lhpc_in.read(record_in_xml, u);
  lhpc_in.close();

  XMLFileWriter xmlout(output_xml);
  push(xmlout, "lhpc2ildg");

  // Compute plaquette
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "plaquette = " << w_plaq << endl;
  write(xmlout, "plaquette", w_plaq);

  // Get code_version from file metadata if it exists
  std::string code_ver;
  if ( file_in_xml.count("/HMC/ProgramInfo/code_version") > 0 ) { 
    std::string chroma_version;
    std::string qdp_version;
    read(file_in_xml, "/HMC/ProgramInfo/code_version/chroma", chroma_version);
    read(file_in_xml, "/HMC/ProgramInfo/code_version/qdp", qdp_version);
    
    std::ostringstream version_ostr;
    version_ostr << "chroma_" << chroma_version << "_qdp_" << qdp_version;
    code_ver = version_ostr.str();
  }
  else { 
    code_ver = "unrecorded";
  }
  
  QDPIO::cout << "code_version = " << code_ver << endl;
  write(xmlout, "code_version", code_ver);
  
  // Get run date if it exists
  std::string run_date;
  if( file_in_xml.count("/HMC/ProgramInfo/run_date") > 0 ) { 
    read(file_in_xml, "/HMC/ProgramInfo/run_date", run_date);
  }
  else { 
    run_date = "unrecorded";
  }
  QDPIO::cout << "run_date = " << run_date <<  endl;
  write(xmlout, "run_date", run_date);
  
  // Dump the LFN
  write(xmlout, "lfn", p.dataLFN);

  // Dump local filename
  write(xmlout, "local", p.ILDG_file_name);

  // Dump output filename
  std::string outfilename=p.ILDG_file_name+".ildg";
  write(xmlout, "outfile", outfilename);

  // Dump config
  XMLBufferWriter file_out;
  XMLBufferWriter record_out;

  // Keep file and record XMLs
  file_out << file_in_xml;
  record_out << record_in_xml;
  
  // Open output file
  QDPFileWriter ildg_out(file_out,  
			 outfilename,
			 QDPIO_SINGLEFILE,
			 QDPIO_SERIAL,
  			 p.dataLFN);


  switch( p.output_size ) { 
  case 32:
    {
      multi1d<LatticeColorMatrixF> u_single_out(Nd);
      for(int mu=0; mu < Nd; mu++) {
	u_single_out[mu] = u[mu];
      }
      ildg_out.write(record_out, u_single_out);
    }
    break;
  case 64:
    {
      multi1d<LatticeColorMatrixD> u_double_out(Nd);
      for(int mu=0; mu < Nd; mu++) {
	u_double_out[mu] = u[mu];
      }
      ildg_out.write(record_out, u_double_out);
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
  XMLReader record_back_in_xml;
  XMLReader file_back_in_xml;
  QDPFileReader ildg_back_in(file_back_in_xml, outfilename, QDPIO_SERIAL);

  multi1d<LatticeColorMatrix> u_back_in(Nd);
  ildg_back_in.read(record_back_in_xml, u_back_in);

  record_in_xml.print(cout);
  cout.flush();

  MesPlq(u_back_in, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "Read Back Plaquette " << w_plaq << endl;

  pop(xmlout);
  xmlout.close();
    // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
