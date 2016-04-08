// -*- C++ -*-

/*! @file
 * @brief XML IO support
 */

#ifndef QDP_XMLIO_H
#define QDP_XMLIO_H

#include <string>
#include <sstream>
#include <stack>
#include <list>
#include <vector>
#include <map>

#include "xml_simplewriter.h"
#include "basic_xpath_reader.h"

namespace QDP 
{

  // Forward declarations
  class XMLReader;
  class XMLWriter;
  class XMLBufferWriter;
  class XMLFileWriter;
  class XMLArrayWriter;


  /*! @ingroup io
   * @{
   */

  //--------------------------------------------------------------------------------
  //! XML reader class
  /*!
    This is used to read data from an XML file using Xpath.

    Note that only the primary node opens and reads XML files. Results from
    Xpath queries are broadcast to all nodes.
  */
  class XMLReader : protected XMLXPathReader::BasicXPathReader
  {
  public:
    //! Empty constructor
    XMLReader();

    //! Construct from contents of file
    /*!
      Opens and reads an XML file.
      \param filename The name of the file.
    */
    XMLReader(const std::string& filename);

    //! Construct from contents of stream
    XMLReader(std::istream& is);

    //! Construct from contents of a XMLBufferWriter
    XMLReader(const XMLBufferWriter& mw);

    //! Clone a reader but with a possibly different path
    XMLReader(XMLReader& old, const std::string& xpath);
    ~XMLReader();

    /* The meaning of these should be clear to you */

    //! Opens and reads an XML file.
    /*!
      \param filename The name of the file
      \post Any previously opened file is closed.
    */
    void open(const std::string& filename);

    //! Opens and reads an XML file.
    /*!
      \param id The input stream of the file
      \post Any previously opened file is closed      
    */
    void open(std::istream& is);

    //! Reads content of a  XMLBufferWriter
    void open(const XMLBufferWriter& mw);

    //! Queries whether the binary file is open
    /*!
      \return true if the binary file is open; false otherwise.
    */
    bool is_open();

    //! Queries whether the XML data has been obtained from another XMLReader
    /*!
      A private method allows this XMLReader to be copy the contents of
      another.
    */
    bool is_derived() const;

    //! Closes the last file opened
    void close();
    
    /* So should these, there is just a lot of overloading */
    //! Xpath query
    void get(const std::string& xpath, std::string& result);
    //! Xpath query
    void get(const std::string& xpath, int& result);
    //! Xpath query
    void get(const std::string& xpath, unsigned int& result);
    //! Xpath query
    void get(const std::string& xpath, short int& result);
    //! Xpath query
    void get(const std::string& xpath, unsigned short int& result);
    //! Xpath query
    void get(const std::string& xpath, long int& result);
    //! Xpath query
    void get(const std::string& xpath, unsigned long int& result);
    //! Xpath query
    void get(const std::string& xpath, float& result);
    //! Xpath query
    void get(const std::string& xpath, double& result);
    //! Xpath query
    void get(const std::string& xpath, bool& result);

    //! read integer attributes so that I can read Matrices in XML
    void getAttribute(const std::string& xpath,
		      const std::string& attrib_name, 
		      int& result);

    //! Set a replacement of a primitive
    template<typename T>
    void set(const std::string& xpath, const T& to_set) 
      {
	if (Layout::primaryNode())
	{  
	  BasicXPathReader::set<T>(xpath, to_set);
	}
      }


    //! Return the entire contents of the Reader as a stream
    void print(std::ostream& is);
        
    //! Print the current context
    void printCurrentContext(std::ostream& is);
        
    //! Count the number of occurances from the Xpath query
    int count(const std::string& xpath);

    void registerNamespace(const std::string& prefix, const std::string& uri);

  private:
    //! Hide the = operator
    void operator=(const XMLReader&) {}
  
    //! Hide the copy constructor
    XMLReader(const XMLReader&) {}
  
    void open(XMLReader& old, const std::string& xpath);
  protected:
    // The universal data-reader. All the read functions call this
    template<typename T>
    void
    readPrimitive(const std::string& xpath, T& output);
  
    // The universal attribute reader
    template<typename T>
    void
    readAttrPrimitive(const std::string& xpath, 
		      const std::string& attrib_name,
		      T& result);
  private:
    bool  iop;  //file open or closed?
    bool  derived; // is this reader derived from another reader?
  };


  // Time to build a telephone book of basic primitives
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, std::string& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, char& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, unsigned int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, short int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, unsigned short int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, long int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, unsigned long int& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, float& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, double& input);
  //! Xpath query
  void read(XMLReader& xml, const std::string& s, bool& input);


  //! Read a XML multi1d element
  template<class T>
  inline
  void read(XMLReader& xml, const std::string& s, multi1d<T>& input)
  {
    XMLReader arraytop(xml, s);

    std::ostringstream error_message;
  
    // Count the number of elements
    std::string elem_base_query = "elem";
	
    int array_size;
    try {
      array_size = arraytop.count(elem_base_query);
    }
    catch( const std::string& e) { 
      error_message << "Exception occurred while counting " << elem_base_query 
		    << " during array read " << std::endl;
    }
      
    // Now resize the array to hold the no of elements.
    input.resize(array_size);

    // Get the elements one by one
    for(int i=0; i < input.size(); i++) 
    {
      std::ostringstream element_xpath;

      // Create the query for the element 
      element_xpath << elem_base_query << "[" << (i+1) << "]";

      // recursively try and read the element.
      try 
      {
	read(arraytop, element_xpath.str(), input[i]);
      } 
      catch (const std::string& e) 
      {
	error_message << "Failed to match element " << i
		      << " of array with query " << element_xpath.str()
		      << std::endl
		      << "Query returned error: " << e;
	throw error_message.str();
      }
    }
  }


  // Specialized versions for basic types
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<unsigned int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<short int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<unsigned short int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<long int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<unsigned long int>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<float>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<double>& input);
  template<>
  void read(XMLReader& xml, const std::string& s, multi1d<bool>& input);


  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //! Read a XML Array element
  template<class T>
  inline
  void read(XMLReader& xml, const std::string& s, std::vector<T>& input)
  {
    XMLReader arraytop(xml, s);

    std::ostringstream error_message;
    std::string elemName = "elem";
  
    // Count the number of elements
    std::string elem_base_query = elemName;
	
    int array_size;
    try {
      array_size = arraytop.count(elem_base_query);
    }
    catch( const std::string& e) { 
      error_message << "Exception occurred while counting " << elem_base_query 
		    << " during array read " << s << std::endl;
      arraytop.close();
      throw error_message.str();
    }
      
    // Now resize the array to hold the no of elements.
    input.resize(array_size);

    // Get the elements one by one
    for(int i=0; i < input.size(); i++) 
    {
      std::ostringstream element_xpath;

      // Create the query for the element 
      element_xpath << elem_base_query << "[" << (i+1) << "]";

      // recursively try and read the element.
      try {
	read(arraytop, element_xpath.str(), input[i]);
      } 
      catch (const std::string& e) 
      {
	error_message << "Failed to match element " << i
		      << " of array  " << s << "  with query " << element_xpath.str()
		      << std::endl
		      << "Query returned error: " << e;
	arraytop.close();
	throw error_message.str();
      }
    }

    // Arraytop should self destruct but just to be sure.
    arraytop.close();
  }


  // Specialized versions for basic types
  void read(XMLReader& xml, const std::string& s, std::vector<int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<unsigned int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<short int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<unsigned short int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<long int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<unsigned long int>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<float>& input);
  void read(XMLReader& xml, const std::string& s, std::vector<double>& input);
//  void read(XMLReader& xml, const std::string& s, std::vector<bool>& input);  // does not seem to exist



  //---------------------------------------------------------------
  //! Read a XML Array element
  template<class T>
  inline
  void read(XMLReader& xml, const std::string& s, std::list<T>& input)
  {
    XMLReader arraytop(xml, s);

    std::ostringstream error_message;
    std::string elemName = "elem";
  
    // Count the number of elements
    std::string elem_base_query = elemName;
	
    int array_size;
    try {
      array_size = arraytop.count(elem_base_query);
    }
    catch( const std::string& e) { 
      error_message << "Exception occurred while counting " << elem_base_query 
		    << " during array read " << s << std::endl;
      arraytop.close();
      throw error_message.str();
    }
      
    // Clear out the original list
    input.clear();

    // Get the elements one by one
    for(int i=0; i < input.size(); i++) 
    {
      std::ostringstream element_xpath;

      // Create the query for the element 
      element_xpath << elem_base_query << "[" << (i+1) << "]";

      // recursively try and read the element.
      try {
	T thingy;
	read(arraytop, element_xpath.str(), thingy);

	input.push_back(thingy);
      } 
      catch (const std::string& e) 
      {
	error_message << "Failed to match element " << i
		      << " of array  " << s << "  with query " << element_xpath.str()
		      << std::endl
		      << "Query returned error: " << e;
	arraytop.close();
	throw error_message.str();
      }
    }

    // Arraytop should self destruct but just to be sure.
    arraytop.close();
  }


  // Specialized versions for basic types
  void read(XMLReader& xml, const std::string& s, std::list<int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<unsigned int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<short int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<unsigned short int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<long int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<unsigned long int>& input);
  void read(XMLReader& xml, const std::string& s, std::list<float>& input);
  void read(XMLReader& xml, const std::string& s, std::list<double>& input);



  //---------------------------------------------------------------
  //! Read a XML Array1dO element
  template<class T>
  inline
  void read(XMLReader& xml, const std::string& s, Array1dO<T>& input)
  {
    read(xml, s, input.ref());
  }


  //---------------------------------------------------------------
  //! Read a XML Array element
  template<typename K, typename V>
  inline
  void read(XMLReader& xml, const std::string& s, std::map<K,V>& input)
  {
    XMLReader arraytop(xml, s);

    std::ostringstream error_message;
    std::string elemName = "elem";
  
    int array_size;
    try {
      array_size = arraytop.count(elemName);
    }
    catch( const std::string& e) { 
      error_message << "Exception occurred while counting " << elemName
		    << " during array read " << s << std::endl;
      arraytop.close();
      throw error_message.str();
    }
      
    // Get the elements one by one
    for(int i=0; i < array_size; i++) 
    {
      std::ostringstream element_xpath;

      // Create the query for the element 
      element_xpath << elemName << "[" << (i+1) << "]";

      // recursively try and read the element.
      try {
	XMLReader xml_elem(arraytop, element_xpath.str());

	K key;
	V val;

	read(xml_elem, "Key", key);
	read(xml_elem, "Val", val);

	input.insert(std::make_pair(key, val));
      } 
      catch (const std::string& e) 
      {
	error_message << "Failed to match element " << i
		      << " of array  " << s << "  with query " << element_xpath.str()
		      << std::endl
		      << "Query returned error: " << e;
	arraytop.close();
	throw error_message.str();
      }
    }

    // Arraytop should self destruct but just to be sure.
    arraytop.close();
  }



  //---------------------------------------------------------------
  //! Read a XML Array element
  template<typename T1, typename T2>
  inline
  void read(XMLReader& xml, const std::string& s, std::pair<T1,T2>& input)
  {
    XMLReader paramtop(xml, s);

    read(paramtop, "First", input.first);
    read(paramtop, "Second", input.second);
  }



  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  //! Metadata output class
  /*!
    Use this to write XML.When closing tags, you do not have to specify which
    tag to close since this class will remember the order in which you opened
    the tags and close them in reverse order to ensure well-formed XML.

    Note that only the primary node writes XML.
  */
  class XMLWriter : protected XMLWriterAPI::XMLSimpleWriter
  {
  public:
    XMLWriter();

    // Virtual destructor
    virtual ~XMLWriter();

    //! Writes an opening XML tag
    /*!
      \param tagname The name of the tag
    */
    virtual void openSimple(const std::string& tagname);
    virtual void closeSimple();

    //! Writes an opening XML tag    
    /*!
      \param tagname The name of the tag
    */
    virtual void openStruct(const std::string& tagname);
    virtual void closeStruct();

    //! Writes an opening XML tag    
    /*!
      \param tagname The name of the tag
    */
    void openTag(const std::string& tagname);

    //! Writes an opening XML tag    
    /*!
      \param nsprefix A namespace prefix for the tag 
      \param tagname The name of the tag
    */
    void openTag(const std::string& nsprefix, const std::string& tagname);

    //! Writes an opening XML tag    
    /*!
      \param tagname The name of the tag
      \param al A list of attributes for this tag
    */
    void openTag(const std::string& tagname, XMLWriterAPI::AttributeList& al);

    //! Writes an opening XML tag    
    /*!
      \param nsprefix A namespace prefix for the tag 
      \param tagname The name of the tag
      \param al A list of attributes for this tag      
    */
    void openTag(const std::string& nsprefix,
		 const std::string& tagname, 
		 XMLWriterAPI::AttributeList& al);

    //! Closes a tag
    void closeTag();

    //! Writes an empty tag
    /*!
      \param tagname The name of the tag
    */
    void emptyTag(const std::string& tagname);

    //! Writes an empty tag
    /*!
      \param nsprefix A namespace prefix for the tag 
      \param tagname The name of the tag
    */
    void emptyTag(const std::string& nsprefix, const std::string& tagname);

    //! Writes an empty tag
    /*!
      \param tagname The name of the tag
      \param al A list of attributes for this tag            
    */
    void emptyTag(const std::string& tagname, XMLWriterAPI::AttributeList& al);

    //! Writes an empty tag
    /*!
      \param nsprefix A namespace prefix for the tag 
      \param tagname The name of the tag
      \param al A list of attributes for this tag            
    */
    void emptyTag(const std::string& nsprefix,
		  const std::string& tagname, 
		  XMLWriterAPI::AttributeList& al);
    

    // Overloaded Writer Functions

    //! Write tag contents
    void write(const std::string& output);
    //! Write tag contents
    void write(const int& output);
    //! Write tag contents
    void write(const unsigned int& output);
    //! Write tag contents
    void write(const short int& output);
    //! Write tag contents
    void write(const unsigned short int& output);
    //! Write tag contents
    void write(const long int& output);
    //! Write tag contents
    void write(const unsigned long int& output);
    //! Write tag contents
    void write(const float& output);
    //! Write tag contents
    void write(const double& output);
    //! Write tag contents
    void write(const bool& output);

    // Write all the XML to std::string
    void writeXML(const std::string& output);

    friend class XMLArrayWriter;
  };


  //! Push a group name
  /*! Write an opening tag
    \param xml The writer
    \param s the name of the tag
  */
  void push(XMLWriter& xml, const std::string& s);

  //! Pop a group name
  /*! Write an closing tag */

  void pop(XMLWriter& xml);

  //! Write something from a reader
  void write(XMLWriter& xml, const std::string& s, const XMLReader& d);
  XMLWriter& operator<<(XMLWriter& xml, const XMLReader& d);

  //! Write something from a XMLBufferWriter
  void write(XMLWriter& xml, const std::string& s, const XMLBufferWriter& d);
  XMLWriter& operator<<(XMLWriter& xml, const XMLBufferWriter& d);

  // Time to build a telephone book of basic primitives
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const std::string& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const char* output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const char& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const unsigned int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const short int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const unsigned short int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const long int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const unsigned long int& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const float& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const double& output);
  //! Write a opening tag, contents and a closing tag
  /*!
    \param xml The writer
    \param s the tag name
    \param output The  contents
  */
  void write(XMLWriter& xml, const std::string& s, const bool& output);

  // Versions that do not print a name
  XMLWriter& operator<<(XMLWriter& xml, const std::string& output);
  XMLWriter& operator<<(XMLWriter& xml, const char* output);
  XMLWriter& operator<<(XMLWriter& xml, const char& output);
  XMLWriter& operator<<(XMLWriter& xml, const int& output);
  XMLWriter& operator<<(XMLWriter& xml, const unsigned int& output);
  XMLWriter& operator<<(XMLWriter& xml, const short int& output);
  XMLWriter& operator<<(XMLWriter& xml, const unsigned short int& output);
  XMLWriter& operator<<(XMLWriter& xml, const long int& output);
  XMLWriter& operator<<(XMLWriter& xml, const unsigned long int& output);
  XMLWriter& operator<<(XMLWriter& xml, const float& output);
  XMLWriter& operator<<(XMLWriter& xml, const double& output);
  XMLWriter& operator<<(XMLWriter& xml, const bool& output);

  //! Write a opening tag, array contents and a closing tag
  /*!
    Each element of the array is written in a "elem" tag.
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const multi1d<T>& s1)
  {
    // Write the array name
    xml.openTag(s);

#if 0
    // This stuff is for schemas
    XMLWriterAPI::AttributeList alist;
    alist.push_back(XMLWriterAPI::Attribute("minOccurs", s1.size()));
    alist.push_back(XMLWriterAPI::Attribute("maxOccurs", s1.size()));
      
    xml.openTag("complexType");
    xml.openTag("sequence", alist);
#endif

    for(int index=0; index < s1.size(); index++) 
    {
      write(xml, "elem", s1[index]);  // Possibly grab user defines here
    }

#if 0
    // This stuff is for schemas
    xml.closeTag(); // sequence
    xml.closeTag(); // complexType
#endif
    xml.closeTag(); // Array name
  }


  // Writers for arrays of basic types
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<unsigned int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<short int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<unsigned short int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<long int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<unsigned long int>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<float>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<double>& output);
  //! Writes an array of data
  /*!
    All the data are written inside a single tag pair
    \param xml The writer
    \param s the tag name
    \param s1 The array of contents
  */
  template<>
  void write(XMLWriter& xml, const std::string& s, const multi1d<bool>& output);

  //! Write an expression
  template<class RHS, class C>
  void write(XMLWriter& xml, const std::string& s, const QDPExpr<RHS,C>& d)
  {
    return write(xml, s, C(d));
  }

  //! XML OScalar output
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const OScalar<T>& d)
  {
    xml.openTag(s);
    xml << d;
    xml.closeTag();
  }

  //! XML OLattice output
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const OLattice<T>& d)
  {
    xml.openTag(s);
    xml << d;
    xml.closeTag();
  }


#if 0
  // NEED TO FIX THIS
  //! Write a XML multi2d element
  template<class T> 
  inline
  void write(XMLWriter& xml, const std::string& s, const multi2d<T>& s1)
  {
    for(int j=0; j < s1.size1(); ++j)
      for(int i=0; i < s1.size2(); ++i)
      {
	std::ostringstream ost;
	if (Layout::primaryNode()) 
	  ost << s << "[ " << i << " ][ " << j << " ]";
	write(xml, ost.str(), s1[i][j]);
      }
  }

#endif



  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //! Write a XML vector element
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const std::vector<T>& s1)
  {
    // Write the array name
    xml.openTag(s);

    for(unsigned index=0; index < s1.size(); index++) 
    {
      write(xml, "elem", s1[index]);  // Possibly grab user defines here
    }

    xml.closeTag(); // Array name
  }


  // Writers for arrays of basic types
  void write(XMLWriter& xml, const std::string& s, const std::vector<int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<unsigned int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<short int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<unsigned short int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<long int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<unsigned long int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<float>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<double>& output);
  void write(XMLWriter& xml, const std::string& s, const std::vector<bool>& output);


  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //! Write a XML list element
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const std::list<T>& s1)
  {
    // Write the array name
    xml.openTag(s);

    for(typename std::list<T>::const_iterator t=s1.begin(); t != s1.end(); ++t)
    {
      write(xml, "elem", *t);  // Possibly grab user defines here
    }

    xml.closeTag(); // Array name
  }


  // Writers for arrays of basic types
  void write(XMLWriter& xml, const std::string& s, const std::list<int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<unsigned int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<short int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<unsigned short int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<long int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<unsigned long int>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<float>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<double>& output);
  void write(XMLWriter& xml, const std::string& s, const std::list<bool>& output);


  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //! Write a XML Array1dO element
  template<class T>
  inline
  void write(XMLWriter& xml, const std::string& s, const Array1dO<T>& s1)
  {
    write(xml, s, s1.ref());
  }


  //---------------------------------------------------------------
  //---------------------------------------------------------------
  //! Write a map
  template<typename K, typename V>
  inline
  void write(XMLWriter& xml, const std::string& path, const std::map<K,V>& param)
  {
    push(xml, path);

    for(typename std::map<K,V>::const_iterator pp = param.begin(); 
	pp != param.end(); 
	++pp)
    {
      push(xml, "elem");
      write(xml, "Key", pp->first);
      write(xml, "Val", pp->second);
      pop(xml);
    }

    pop(xml);
  }


  //! Write a pair
  template<typename T1, typename T2>
  inline
  void write(XMLWriter& xml, const std::string& path, const std::pair<T1,T2>& param)
  {
    push(xml, path);

    write(xml, "First", param.first);
    write(xml, "Second", param.second);

    pop(xml);
  }


  //--------------------------------------------------------------------------------
  //! Writes XML metadata to a buffer
  class XMLBufferWriter : public XMLWriter
  {
  public:

    /*! No prologue written */
    XMLBufferWriter();
  
    //! Construct from a string
    explicit XMLBufferWriter(const std::string& s);

    //! Destroy
    ~XMLBufferWriter();

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    // Return root element as a string
    std::string printCurrentContext() const;
        
    //! Flush the buffer
    void flush() {}

    //! Return true if some failure occurred in previous IO operation
    bool fail() const {return false;}

  private:
    // The output stream...
    std::ostringstream output_stream;

    // The function that supplies the stream to the parent...
    std::ostream& getOstream(void) {return output_stream;}
  };


  //--------------------------------------------------------------------------------
  //! Writes XML metadata to a file
  /*!
    \ingroup io
  */

  class XMLFileWriter : public XMLWriter
  {
  public:

    XMLFileWriter();

    //! Constructor from a filename
    /*!
      \param filename The name of the file
      \param write_prologue Whether to write the standard opening line of
      XML files. Defaults to true.
    */
    explicit XMLFileWriter(const std::string& filename, bool write_prologue=true)
      {
	open(filename, write_prologue);
      }

    ~XMLFileWriter();

    //! Queries whether the binary file is open
    /*!
      \return true if the binary file is open; false otherwise.
    */
    bool is_open();

    //!Opens a file
    /*!
      \param filename The name of the file
      \param write_prologue Whether to write the standard opening line of
      XML files. Defaults to true.
    */
    void open(const std::string& filename, bool write_prologue=true);

    //! Flush the buffer
    void flush();

    //! Return true if some failure occurred in previous IO operation
    bool fail() const;

    //! Closes the last  file  opened.
    void close();
        
  private:
    std::ofstream output_stream;
    std::ostream& getOstream(void) {return output_stream;}
  };



  //--------------------------------------------------------------------------------
  //! Writes metadata to an array which serves as a handle for another XML object
  class XMLArrayWriter : public XMLWriter
  {
  public:

    /*! No prologue written
     * @param xml_out  previous XMLWriter object - used for handle source
     * @param size     size of array - default unbounded
     */
    explicit XMLArrayWriter(XMLWriter& xml_out, int size = -1) : 
      output_xml(xml_out), array_size(size)
      {
	initP = arrayTag = false;
	elements_written = 0;
	indent_level = xml_out.indent_level;
	simpleElements = false; // do not know this yet
      }
  

    ~XMLArrayWriter();

    // Flush the buffer
    //  void flush();

    //! Closes the array writer
    /*! It is an error to close before all data is written, unless unbounded */
    void close();
       
    //! Returns the size of the array
    int size() const {return array_size;}

    void openArray(const std::string& tagname);
    void closeArray();

    //  virtual void openSimple(const std::string& tagname);
    //  virtual void closeSimple();

    void openStruct(const std::string& tagname);
    void closeStruct();

  private:
    std::string qname;
    std::string elem_qname;

    bool arrayTag;         // set once array tag is written
    bool initP;            // set once we know how the array is composed
    bool simpleElements;   // true if elements will all be simple types

    // output stream is actually the original stream
    XMLWriter& output_xml; 

    int array_size;        // total array element size
    int elements_written;  // elements written so far

    // A stack to hold context.
    enum ElementType {SIMPLE, STRUCT};
    std::stack<ElementType> contextStack;

    std::ostream& getOstream(void) {return output_xml.getOstream();}
  };

  //! Push a group name
  void push(XMLArrayWriter& xml);

  //! Pop a group name
  void pop(XMLArrayWriter& xml);


  /*! @} */   // end of group io

} // namespace QDP

#endif
