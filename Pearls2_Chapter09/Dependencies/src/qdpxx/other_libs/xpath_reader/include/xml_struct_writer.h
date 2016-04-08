#ifndef XML_STRUCT_WRITER_H
#define XML_STRUCT_WRITER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <xml_array_type.h>
#include <cstdlib>

#define INDENT_SPACES 2


namespace XMLStructWriterAPI { 

  class XMLFileArrayWriter;

  // Base class for struct writers
  class XMLStructWriterBase {
  public:
    XMLStructWriterBase(const std::string& tagname, int _indent_level=0) 
    { 
      indent_level = _indent_level;
      qname=tagname;
    }

    virtual ~XMLStructWriterBase() {};

    void writeSimple(const std::string& tagname, const int& value) {
      writeSimpleTag(tagname, value);
    }

    void writeSimple(const std::string& tagname, const float& value) {
      writeSimpleTag(tagname, value);
    }
    
    void writeSimple(const std::string& tagname, const std::string& value) {
      writeSimpleTag(tagname, value);
    }

    void writeSimple(const std::string& tagname, const bool& value) {
      writeSimpleTag(tagname, value);
    }  

    template <typename T>
      void writeSimpleTag(const std::string& tagname, T& value)
      {
	std::ostream &os=getOstream();
	std::string indent((indent_level+1)*INDENT_SPACES, ' ');
	os << std::endl << indent << "<" << tagname << ">";
	os << std::boolalpha <<  value << "</" << tagname << ">";
      }  

  protected:
    
    virtual std::ostream& getOstream() = 0;

    void writePrologue(std::ostream& os) const {
      os << "<?xml version=\"1.0\"?>" << std::endl;
      os << "<!-- Written by XMLSimpleWriter class by Balint Joo -->";
      os << std::endl;
      os.flush();
    }

    bool simple_array;
    int indent_level;
    std::string qname;
    int elems_written;
  };

  class XMLFileStructWriter : public XMLStructWriterBase {
  public:
    XMLFileStructWriter(std::ofstream& _os, 
			const std::string& tagname, 
			int indent_level=0, 
			bool write_prologue=false) : 
      XMLStructWriterBase(tagname, indent_level), output_stream(_os) {
      
      if ( write_prologue == true ) { 
	writePrologue(_os);
      }

      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "<" << qname << ">";

    }

    ~XMLFileStructWriter(void) {
      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "</" << qname << ">";
    }

    XMLFileStructWriter* structChild(const std::string& tagname) {
      XMLFileStructWriter* child=new(std::nothrow) XMLFileStructWriter(output_stream, tagname, indent_level+1, false);

      if( child == 0x0 ) { 
	std::cerr << "Failed to allocate child: structChild" << std::endl << std::flush;
	exit(-1);
      }

      return child;
    }

    
    XMLFileArrayWriter* arrayChild(const std::string& tagname, const std::string& elem_name, ArrayType t);
    
  private:
    std::ofstream& output_stream;
    std::ostream& getOstream(void) { 
      return output_stream;
    }
    
  };

  class XMLBufferStructWriter : public XMLStructWriterBase {
  public:
    XMLBufferStructWriter(std::ostringstream& _os, 
			  const std::string& tagname, 
			  int indent_level=0, 
			  bool write_prologue=false) : 
      XMLStructWriterBase(tagname, indent_level), output_stream(_os) {
      
      if ( write_prologue == true ) { 
	writePrologue(_os);
      }

      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "<" << qname << ">";

    }

    ~XMLBufferStructWriter(void) {
      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "</" << qname << ">";
    }

    XMLBufferStructWriter* structChild(const std::string& tagname) {
      XMLBufferStructWriter* child=new(std::nothrow) XMLBufferStructWriter(output_stream, tagname, indent_level+1, false);
      if( child == 0x0 ) {
	std::cerr << "Failed to allocate child: structChild() " << std::endl << std::flush ;
	exit(-1);
      }

      return child;
    }
  
  private:
    std::ostringstream& output_stream;
    std::ostream& getOstream(void) { 
      return output_stream;
    }
    
  };

};



#endif
