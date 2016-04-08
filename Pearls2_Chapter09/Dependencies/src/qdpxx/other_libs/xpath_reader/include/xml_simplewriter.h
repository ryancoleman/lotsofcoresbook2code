#ifndef XML_SIMPLEWRITER_H
#define XML_SIMPLEWRITER_H

#include <xml_attribute.h>
#include <ostream>
#include <string>
#include <stack>


namespace XMLWriterAPI {

  // Base class for the XMLSimpleWriter classes
  class XMLSimpleWriter {
  public:
    XMLSimpleWriter() {
      doctag_written = false;
      primitive_last = false;
    }

    virtual ~XMLSimpleWriter() {};

    void openTag(const std::string& tagname);
    void openTag(const std::string& nsprefix, const std::string& tagname);
    void openTag(const std::string& tagname, AttributeList& al);

    void openTag(const std::string& nsprefix,
		 const std::string& tagname, 
		 AttributeList& al);

    void closeTag();

    void emptyTag(const std::string& tagname);
    void emptyTag(const std::string& nsprefix, const std::string& tagname);
    void emptyTag(const std::string& tagname, AttributeList& al);

    void emptyTag(const std::string& nsprefix,
		  const std::string& tagname, 
		  AttributeList& al);
    

   

    // Overloaded Writer Functions
    void write(const std::string& output);
    void write(const int& output);
    void write(const unsigned int& output);
    void write(const short int& output);
    void write(const unsigned short int& output);
    void write(const long int& output);
    void write(const unsigned long int& output);
    void write(const float& output);
    void write(const double& output);
    void write(const bool& output);

    // Write XML string
    void writeXML(const std::string& output);


  protected:
    // Protected functions -- for overloading access by derived classes
    
    // Get the internal ostream
    virtual std::ostream& getOstream() = 0;

    // Write the prologue
    void writePrologue(std::ostream& os) const;
    unsigned int indent_level;

  private:
    // Check we have written a document tag...
    // apparently we need this before writing anything else
    bool doctag_written;
    bool primitive_last;

    // A stack to hold names.
    std::stack<std::string> namestack;

    // The universal tag-write function. Copes with empty tags too
    // all the openTag, emptyTag functions call this basically
    void dumpTag(const std::string& nsprefix,
		 const std::string& tagname, 
		 AttributeList& al,
		 bool is_empty=false);
     
    // The universal data-write. All the write functions call this
    template< typename T>
      void
      writePrimitive(const T& output);


    void
      writePrimitive(const float& output);
    

    void
      writePrimitive(const double& output);
  };

};

#endif



