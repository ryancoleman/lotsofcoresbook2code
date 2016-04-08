#ifndef XML_WRITER_H
#define XML_WRITER_H

#include <xml_attribute.h>
#include <xml_tcomplex.h>
#include <xml_array.h>
#include <xml_simplewriter.h>
#include <sstream>
#include <fstream>

namespace XMLWriterAPI {

  // Namespace composition. 
  using XMLArray::Array;
  using XMLTComplex::TComplex;


  // Base class for the XMLWriter classes
  class XMLWriterBase : public XMLSimpleWriter {
  public:
    XMLWriterBase(void) {}
    ~XMLWriterBase(void) {}

    // Overloaded Writer Functions
    /* These wrappers are needed since there is overloading of the name "write"
     * below for complex and array. The openTag's do not need to be re-overloaded */
    void write(const std::string& output) {XMLSimpleWriter::write(output);}
    void write(const int& output) {XMLSimpleWriter::write(output);}
    void write(const unsigned int& output) {XMLSimpleWriter::write(output);}
    void write(const short int& output) {XMLSimpleWriter::write(output);}
    void write(const unsigned short int& output) {XMLSimpleWriter::write(output);}
    void write(const long int& output) {XMLSimpleWriter::write(output);}
    void write(const unsigned long int& output) {XMLSimpleWriter::write(output);}
    void write(const float& output) {XMLSimpleWriter::write(output);}
    void write(const double& output) {XMLSimpleWriter::write(output);}
    void write(const bool& output) {XMLSimpleWriter::write(output);}
   
    // Additional overloaded Writer Functions
    template <typename T>
      void write(TComplex<T>& output) {
      openTag((std::string)"cmpx");
      openTag((std::string)"re");
      write(output.real());
      closeTag();
      openTag("im");
      write(output.imag());
      closeTag();
      closeTag();
    }

    template <typename T> 
      void write(const std::string& sizeName, 
		 const std::string& elemName, 
		 const std::string& indexName,
		 const unsigned int& indexStart,
		 Array<T>& a) {
      
      AttributeList alist;
      alist.push_back(Attribute((std::string)"sizeName",  sizeName));
      alist.push_back(Attribute((std::string)"elemName",  elemName));
      alist.push_back(Attribute((std::string)"indexName", indexName));
      alist.push_back(Attribute((std::string)"indexStart", indexStart));
      
      // Write the array - tag
      openTag((std::string)"array", alist);

      openTag(sizeName);
      write(a.size());
      closeTag();

      unsigned int index;
      for(index=0; index < a.size(); index++) {
	alist.clear();
	alist.push_back(Attribute(indexName, index + indexStart));
	openTag(elemName, alist);
	write(a[index]);
	closeTag();
      }

      closeTag(); // Array
    }

    template < typename T > 
    void write(Array<T>& a) {
      write((std::string)"size",
	    (std::string)"element",
	    (std::string)"index",
	    0, a);
    }
  };

  // A Class to write XML to a std::string
  // Supplies an ostd::stringstream as the outstream for the BaseClass
  class XMLSimpleStringWriter : public XMLWriterBase {
  public:

    // Constructor
    XMLSimpleStringWriter(bool write_prologue=true) { 
      if( write_prologue ) { 
	writePrologue(output_stream);
      }
      indent_level=0;
    }  

    ~XMLSimpleStringWriter(void) {
      bool done=false;
      
      while( ! done ) { 
	try { 
	  closeTag();
	}
	catch(std::string& error) { 
	  done = true;
	}
      }
    }

    // Get the string -- not the stream...
    std::string str(void) const { 
      return output_stream.str();
    }
        
  private:

    // The output stream...
    std::ostringstream output_stream;

    // The function that supplies the stream to the parent...
    std::ostream& getOstream(void) { 
      return output_stream;
    }
  };

 
  class XMLSimpleFileWriter : public XMLWriterBase {
  public:
    XMLSimpleFileWriter(const std::string& filename, bool write_prologue=true) {
    
      output_stream.open(filename.c_str(), std::ofstream::out);
      if( write_prologue ) {
	writePrologue(output_stream);
      }
      indent_level=0;
    }
        
    void close() 
      {
	output_stream.close();
      }

    ~XMLSimpleFileWriter(void) {
      bool done=false;
      
      while( ! done ) { 
	try { 
	  closeTag();
	}
	catch(std::string& error) { 
	  done = true;
	}
      }

      output_stream.close();
    }

  private:
    std::ofstream output_stream;
    std::ostream& getOstream(void) { 
      return output_stream;
    }
  };
  
};

#endif



