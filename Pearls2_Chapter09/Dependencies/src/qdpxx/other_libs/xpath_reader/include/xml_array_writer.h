#ifndef XML_ARRAY_WRITER_H
#define XML_ARRAY_WRITER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <xml_array_type.h>



#define INDENT_SPACES 2

namespace XMLStructWriterAPI { 

  class XMLFileArrayWriter;
  class XMLFileStructWriter;
 
  // Base class for struct writers
  class XMLArrayWriterBase {
  public:
    XMLArrayWriterBase(const std::string& tagname, 
		       const std::string& elem_name, 
		       ArrayType atype, 
		       int _indent_level=0) 
    { 
      indent_level = _indent_level;
      qname=tagname;
      elem_qname=elem_name;
      elems_written=0;
      type=atype;
    }

    virtual ~XMLArrayWriterBase() { };

    void elementSimple(const int& value) {
      if (type == SIMPLE_INT) { 
	writeSimpleTag(value);
      }
      else {
	throw "Array member type mismatch error: wanted SIMPLE_INT";
      }
    }

    void elementSimple(const float& value) {
      if (type == SIMPLE_FLOAT) { 
	writeSimpleTag(value);
      }
      else {
	throw "Array member type mismatch error: wanted SIMPLE_FLOAT";
      }
    }

    void elementSimple(const std::string& value) {
      if (type == SIMPLE_STRING) { 
	writeSimpleTag(value);
      }
      else {
	throw "Array member type mismatch error: wanted SIMPLE_STRING";
      }
    }

    void elementSimple(const bool& value) {
      if (type == SIMPLE_BOOL) { 
	writeSimpleTag(value);
      }
      else {
	throw "Array member type mismatch error: wanted SIMPLE_BOOL";

      }
    }

    template <typename T>
      void writeSimpleTag(T& value)
      {
	std::ostream &os=getOstream();
	std::string indent((indent_level+1)*INDENT_SPACES, ' ');
	os << std::endl << indent << "<" << elem_qname << ">";
	os << std::boolalpha <<  value << "</" << elem_qname << ">";
	elems_written++;
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
    std::string elem_qname;
    ArrayType type;
    int elems_written;
  };

  class XMLFileArrayWriter : public XMLArrayWriterBase {
  public:
    XMLFileArrayWriter(std::ofstream& _os, 
		       const std::string& tagname, 
		       const std::string& elem_name,
		       ArrayType atype,
		       int indent_level=0, 
		       bool write_prologue=false) : 
      XMLArrayWriterBase(tagname,elem_name, atype, indent_level), output_stream(_os) {
      
      if ( write_prologue == true ) { 
	writePrologue(_os);
      }
      
      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "<" << qname << ">";
      
    }
    
    ~XMLFileArrayWriter(void) {
      std::string indent(indent_level*INDENT_SPACES, ' ');
      output_stream  << std::endl << indent << "</" << qname << ">";
    }
    
    
    XMLFileArrayWriter* elementArray(ArrayType t) {
      XMLFileArrayWriter* child=new(std::nothrow) XMLFileArrayWriter(
						       output_stream,
						       elem_qname,
						       "element",
						       t,
						       indent_level+1,
						       false);
      if( child == 0x0 ) { 
	std::cerr << " Failed to allocate child " << std::endl << std::flush;
	exit(-1);
      }

      elems_written++;
      return child;
    }
    XMLFileArrayWriter* elementArray(const std::string& elem_name, ArrayType t) {
      XMLFileArrayWriter* child=new(std::nothrow) XMLFileArrayWriter(
						       output_stream,
						       elem_qname,
						       elem_name,
						       t,
						       indent_level+1,
						       false);
      if( child == 0x0 ) { 
	std::cerr << " Failed to allocate child " << std::endl << std::flush;
	exit(-1);
      }

      elems_written++;
      return child;
    }

    XMLFileStructWriter* elementStruct(void); 
    
  private:
    std::ofstream& output_stream;
    std::ostream& getOstream(void) { 
      return output_stream;
    }
    
  };

  /*
  class XMLBufferStructWriter : public XMLStructWriterBase {
  public:
    XMLBufferStructWriter(ostringstream& _os, 
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
        std::cerr << "Failed to new child " << std::endl << std::flush;
	exit(-1);
      }
      return child;
    }
  
  private:
    ostringstream& output_stream;
    std::ostream& getOstream(void) { 
      return output_stream;
    }
    
  };
  */
};



#endif
