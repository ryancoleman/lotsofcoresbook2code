#include <fstream>
#include <ostream>
#include <istream>
#include <string>

#include <xml_array_type.h>
#include <xml_struct_writer.h>
#include <xml_array_writer.h>

using namespace std;

namespace XMLStructWriterAPI {

XMLFileArrayWriter* XMLFileStructWriter::arrayChild(const string& tagname, 
						    const string& elem_name, 
						    ArrayType t)
{
  XMLFileArrayWriter* child = new(nothrow) XMLFileArrayWriter(output_stream, 
		  				    tagname, elem_name, 
						    t, indent_level+1, false);
  if( child == 0x0 ) { 
    std::cerr << "Failed to allocate child: ArrayChild " << endl <<flush;
    exit(-1);
  }

  return child;
}

XMLFileStructWriter* XMLFileArrayWriter::elementStruct(void) 
{
  XMLFileStructWriter* child=new(nothrow) XMLFileStructWriter(output_stream,
		 				    elem_qname, indent_level+1, false); 

  if( child == 0x0 ) { 
    std::cerr << "Failed to allocate child: elementStruct() " << endl << flush;
    exit(-1);
  }

  elems_written++; 
  return child;
}

};
