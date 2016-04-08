/* ID: $Id: basic_xpath_reader.cc,v 1.20 2009-08-25 15:41:16 colin Exp $
 *
 * File: basic_xpath_reader.cc
 * 
 * Purpose: All the nasty libxml stuff should be hidden in here 
 * as well as a lot of the overloading/wrappint in BasicXPathReader
 */

#include <basic_xpath_reader.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <ios>
#include <iostream>
#include <sstream>

/* need this for iostream and string */
using namespace std;

namespace XMLXPathReader {
  bool xpath_is_initialised = false;

  void initXPath(void) {
    xmlXPathInit();
    xpath_is_initialised = true;
  }
};
/* This is our namespace */
using namespace XMLXPathReader;


BasicXPathReader::BasicXPathReader(void) 
{
  if( xpath_is_initialised != true ) { 
    initXPath();
  }

  docref = 0x0;
  xpath_context = NULL;
  query_result = NULL;
}

BasicXPathReader::BasicXPathReader(std::istream& is) 
{ 
  if( xpath_is_initialised != true ) 
  { 
    initXPath();
  }

  docref = 0x0;
  xpath_context = NULL;
  query_result = NULL; 
  try {
    open(is);
  }
  catch( const std::string& e) {
    throw;
  }
}
    
BasicXPathReader::BasicXPathReader(const std::string& filename)  
{ 
  if( xpath_is_initialised != true ) { 
    initXPath();
  }

  docref = 0x0;
  xpath_context = NULL;
  query_result = NULL;
  try { 
    open(filename);
  }
  catch( const std::string& e) {
    throw;
  }
}
    
BasicXPathReader::~BasicXPathReader(void) 
{
  close();
}
    
/* This should let you clone an XPath reader */
BasicXPathReader::BasicXPathReader(BasicXPathReader& old, const std::string& xpath) 
{

  docref=0x0;
  xpath_context = NULL;
  query_result = NULL;

  open(old, xpath);
}
    

void BasicXPathReader::open(const std::string& filename) 
{
  // There are two possibilities to opening
  //
  // i) Reader is already open
  // ii) The reader is not yet open.
  //
  //
  // case i) Reader is already open.
  close();
      
  // Ok, open the document
  try 
  {
    docref = new(std::nothrow) XMLDocument(filename);
	
    if( docref == 0x0 ) {
      std::ostringstream error_stream;
      error_stream << "Unable to create new XML Document using open(filename) filename=" << filename << ")";
      throw error_stream.str();
    }
	
    // Make sure the document feels referred to.
    docref->increaseRefcount();
    query_result = (xmlXPathObjectPtr) NULL;
    xmlDocPtr doc = docref->getDoc();
    xpath_context = xmlXPathNewContext(doc);
    xpath_context->node = xmlDocGetRootElement(doc);
    // snarfNamespaces(xmlDocGetRootElement(doc), xpath_context);
  }
  catch ( const std::string& e ) 
  {
    throw;
  }
}
    
void BasicXPathReader::open(std::istream &is) 
{
  // There are two possibilities to opening
  //
  // i) Reader is already open
  // ii) The reader is not yet open.
  //
  //
  // case i) Reader is already open.
  close();

  // Ok, open the document
  try 
  {
    docref = new(std::nothrow) XMLDocument(is);
	
    if( docref == 0x0 ) {
      std::ostringstream error_stream;
      error_stream << "Unable to create new XML Document using open(istream& is)";
      throw error_stream.str();
    }
	
    // Make sure the document feels referred to.
    docref->increaseRefcount();
	
    query_result = (xmlXPathObjectPtr) NULL;
    xmlDocPtr doc = docref->getDoc();
    xpath_context = xmlXPathNewContext(doc);
    xpath_context->node = xmlDocGetRootElement(doc);
    // snarfNamespaces(xmlDocGetRootElement(doc), xpath_context);
  }
  catch ( const std::string& e ) 
  {
    throw;
  }
}

/* This should let you clone an XPath reader */
void BasicXPathReader::open(BasicXPathReader& old, const std::string& xpath) 
{
  close();

  if( old.docref != 0x0) 
  {
    docref = old.docref;
	
    query_result = (xmlXPathObjectPtr) NULL;
    xmlDocPtr doc = docref->getDoc();
    xpath_context = xmlXPathNewContext(doc);

    // Now execute the xpath query and set the context node
    std::ostringstream error_message;
	
    try { 
      old.evaluateXPath(xpath);
    }
    catch ( const std::string& e ) { 
      throw e;
    }
	
    // Check that the query returned non-empty result
    try { 
      old.checkQuery(xpath);
    }
    catch( const std::string& e) { 
      throw;
    }
	
    //Check that the node set contains only 1 element
    if( old.query_result->nodesetval->nodeNr != 1 ) 
    {
      error_message << "XPath Query: " << xpath << " did not return unique node."
		    << " nodes returned = " << old.query_result->nodesetval->nodeNr 
		    << std::endl;
      throw error_message.str();
    }
	
    // Check that the node returned is an element
    xpath_context->node = old.query_result->nodesetval->nodeTab[0];

    std::list<NameSpace>::iterator iter;
    for( iter = old.nslist.begin(); iter != old.nslist.end(); iter++) 
    {
      registerNamespace(iter->getPrefix(), iter->getURI());
    }
	
    docref->increaseRefcount();
  }
  else { 
    throw "Attempting to clone a closed Reader";
  }
}
      
void BasicXPathReader::close(void) 
{ 
  if( docref != 0x0 ) 
  {
    // Reader is already open. 
    // We detach from the current document:
	
    // decrement its refcount
    docref->decreaseRefcount();
	
    // if decrementing the refcount means there are no more references
    // then delete the object -- this ought to call the destructor
    if ( docref->getRefcount() == 0 ) 
    {
#ifdef DEBUG_XML_REFCOUNT
      std::cout << "Reader: docrefs refcount reached 0. Deleting" << std::endl;
#endif
      delete docref;
    }
	
    // We are now officially not open 
    docref = 0x0;
	
    // Clean up any left-over query result
    if ( query_result != NULL ) { 
      xmlXPathFreeObject(query_result);
    }
	
    // Clean up the XPath content
    if( xpath_context != NULL ) { 
      xmlXPathFreeContext(xpath_context);
    }

    nslist.clear();
  }
  else {
    // Else we are not open so closing will do sod all
  }      
}


/* So should these, there is just a lot of overloading */
void BasicXPathReader::get(const std::string& xpath, std::string& result) 
{
      
  // Try and get the std::string for the xpath
  try {
    getPrimitiveString(xpath, result);
  }
  catch(const std::string& e) { 
    // Pass up exception if it occurs
    throw e;
  }
}

void BasicXPathReader::get(const std::string& xpath, int& result) 
{ 
  getPrimitive<int>(xpath, result, "int");
}

void BasicXPathReader::get(const std::string& xpath, unsigned int& result) 
{
  getPrimitive<unsigned int>(xpath, result, "unsigned int");
}
    
void BasicXPathReader::get(const std::string& xpath, short int& result) 
{
  getPrimitive<short int>(xpath, result,"short int");
}

void BasicXPathReader::get(const std::string& xpath, unsigned short int& result) 
{
  getPrimitive<unsigned short int>(xpath, result, "unsigned short int");
}

void BasicXPathReader::get(const std::string& xpath, long int& result) 
{
  getPrimitive<long int>(xpath, result, "long int");
}
    
void BasicXPathReader::get(const std::string& xpath, unsigned long int& result) 
{
  getPrimitive<unsigned long int>(xpath, result, "unsigned long int");
}
    
void BasicXPathReader::get(const std::string& xpath, float& result) 
{
  getPrimitive<float>(xpath, result, "float");
}

void BasicXPathReader::get(const std::string& xpath, double& result) 
{
  getPrimitive<double>(xpath, result, "double");
}

void BasicXPathReader::get(const std::string& xpath, bool& result) 
{
  //  getPrimitive<bool>(xpath, result, "bool");
  std::string result_string;
  getPrimitiveString(xpath, result_string);
  if ( result_string == "true" ) { 
    result = true;
  }
  else if( result_string == "false") { 
    result = false;
  }
  else { 
    std::ostringstream error_message;
      
    error_message << "Failed to extract boolean from result string "<< result_string <<" on xpath quary " << xpath << endl;
    throw error_message.str();
  }
          
}

/* So should these, especially if you read the introductory comments */
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    std::string& result) 
{
  getAttributeString(xpath_to_node, attrib_name, result);
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    int& result) 
{
  getPrimitiveAttribute<int>(xpath_to_node, attrib_name, result, "int");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    unsigned int& result) 
{
  getPrimitiveAttribute<unsigned int>(xpath_to_node, attrib_name, result, 
				      "unsigned int");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    short int& result) 
{
  getPrimitiveAttribute<short int>(xpath_to_node, attrib_name, result, 
				   "short int");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    unsigned short int& result) 
{
  getPrimitiveAttribute<unsigned short int>(xpath_to_node, 
					    attrib_name, result, 
					    "unsigned short int");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    long int& result) 
{
  getPrimitiveAttribute<long int>(xpath_to_node, attrib_name, result, 
				  "long int");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    unsigned long int& result) 
{
  getPrimitiveAttribute<unsigned long int>(xpath_to_node, attrib_name, result, 
					   "unsigned long int");
}

void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    float& result) 
{
  getPrimitiveAttribute<float>(xpath_to_node, attrib_name, result, 
			       "float");
}
    
void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    double& result) 
{
  getPrimitiveAttribute<double>(xpath_to_node, attrib_name, result, 
				"double");
}

void BasicXPathReader::getAttribute(const std::string& xpath_to_node,
				    const std::string& attrib_name, 
				    bool& result) 
{
  getPrimitiveAttribute<bool>(xpath_to_node, attrib_name, result, 
			      "bool");
}

int BasicXPathReader::count(const std::string& xpath) 
{
  std::ostringstream error_message;
  int ret_val;
      
  // Construct query 
  std::string query="count("+xpath+")";
      
  // Run query 
  // This is the bit that checks that what we are doing is legit
      
  try { 
    evaluateXPath(query);
  }
  catch(const std::string& e) {
    throw e;
  }
      
      
  // Check that the result is a number otherwise something is very wrong! 
  if( query_result->type != XPATH_NUMBER ) 
  {
    error_message << "XPath Query: " << query
		  << " did not succeed" << std::endl;
    throw error_message.str();
  }
  else 
  {
    // It is  a number -- goodie. Cast to int. This should be OK, 
    // as you cant count a non-integral no of things.
    ret_val = (int)(query_result->floatval);
  }
      
  // Post query cleanup
  if ( query_result != NULL ) 
  {
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }  
      
  // return the number 
  return ret_val;      
}

void BasicXPathReader::evaluateXPath(const std::string& xpath) 
{
  std::ostringstream error_message;
      
  // Guard against queries of empty readers
  if ( docref == 0 ) 
  {
    error_message << "Attempting to execute a query on an empty Reader";
    throw error_message.str();
  }
      
  // Make sure previos queries didn't leave anything kicking about 
  if( query_result != NULL ) 
  {
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }
      
  // Call libxml2 to run the xpath query 
  query_result = xmlXPathEval((const xmlChar *)xpath.c_str(),
			      xpath_context);
      
  // if the query is null throw up our hands in despair 
  if( query_result == NULL ) 
  {
    error_message << "XPath Query " << xpath << " Failed! " << std::endl;
    throw error_message.str();
  }
}


void 
BasicXPathReader::printDoc(ostream &os)
{  
  if( docref != 0x0 ) 
  {
    xmlChar *buffer=(xmlChar *)NULL;
    int buflen;
    ostringstream error_stream;
    xmlDocPtr doc = docref->getDoc();
    
    if( doc != NULL ) 
    {
      xmlDocDumpFormatMemory(doc, &buffer, &buflen, 1);
      if( buffer ==(xmlChar *)NULL ) 
      {
        error_stream << "xmlDocDumpMemory produced NULL XML-char";
        throw error_stream.str();
      }
      buffer[buflen]='\0';
      os << buffer << endl;
      xmlFree(buffer);
    }
  }
}

void
BasicXPathReader::print(ostream& os)
{
  if( docref != 0x0 ) {
    printXPathNode(os, ".");
  }
}      


void
BasicXPathReader::printChildren(ostream& os)
{
  if( docref != 0x0 ) { 
    printXPathNode(os,"./*");
  }
}

void 
BasicXPathReader::printRoot(ostream& os)
{
  if( docref != 0x0 ) { 
    xmlDocPtr doc = docref->getDoc();
    printNode(os,xmlDocGetRootElement(doc) );
  }
}

void 
BasicXPathReader::printXPathNode(ostream& os, const string& xpath_to_node)
{
  ostringstream error_message;

  // Evaluate the Xpath to the node
  try { 
    evaluateXPath(xpath_to_node);
  }
  catch (const string& e) {
    throw e;
  }

  // check for general non nullness
  try { 
    checkQuery(xpath_to_node);
  }
  catch (const string& e) {
    throw e;
  }

  for(int i = 0; i < query_result->nodesetval->nodeNr; i++) 
  {
    /* Check that the node returned is an element */
    xmlNodePtr res_node = query_result->nodesetval->nodeTab[i];
  
    if ( res_node != NULL) 
    {
      // print it 
      printNode(os, res_node);
      os << endl;
    }
  }
  // clean up
  if ( query_result != NULL ) 
  {
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }    
}

void
BasicXPathReader::printNode(ostream& os, xmlNodePtr node)
{
  if( docref != 0x0 ) 
  {
    xmlDocPtr doc = docref->getDoc();
    xmlBufferPtr xmlBuf;
    ostringstream error_message;
    
    xmlBuf = xmlBufferCreate();
    if ( xmlBuf == (xmlBufferPtr)NULL ) 
    {
      error_message << "Failed to create buffer in printNode" << endl;
      throw error_message.str();
    }
    
#if 0
    // RGE: This was the original: updated below
    int size;
    size = xmlNodeDump(xmlBuf, doc, node, 2, 1);
    
    if ( size == -1 ) 
    {
      error_message << "xmlNodeDump failed most heinously" << endl;
      throw error_message.str();
    }
    
    os.write((char *)xmlBufferContent(xmlBuf), size);
#else
    xmlNodeDump(xmlBuf, doc, node, 2, 1);
    
    os << string((char *)xmlBufferContent(xmlBuf));
#endif
    
    xmlBufferFree(xmlBuf);
  }
}



/*! pretty print the query_result object onto an os */
void
BasicXPathReader::printQueryResult(ostream& os)
{
  if( docref != 0x0 ) 
  {
    int nset;
    int i;
    ostringstream error_message;
    
    if( query_result == NULL ) 
    {
      os << "Query Result is NULL" << endl;
    }
    else 
    {
      xmlDocPtr doc = docref->getDoc();
      
      /* Print different stuff depending on the query result type */
      switch(query_result->type) 
      {
      case XPATH_UNDEFINED:
	os << "======== QUERY RESULT IS UNDEFINED =========" << endl;
	break;
      case XPATH_NODESET:
	os << "======== QUERY RESULT IS A NODE_SET ========" << endl;
	os << "====== NODE DUMP START ========" << endl;
	os.flush();
	
	/* Dump a node set -- if it is not null */
	if( query_result->nodesetval == NULL  ) 
	{
	  os << "====== NODE SET IS NULL =======" << endl;
	}
	else 
	{
	  /* Get the number of nodes in the node set */
	  nset = query_result->nodesetval->nodeNr;
	  os << "Nset is " << nset << endl;
	  if( nset > 0 ) 
	  {
	    /* If there is more than one, than go through the nodes
	     * and dump each one */
	    for(i=0; i < nset; i++) 
	    {
	      os.flush();
	      switch(query_result->nodesetval->nodeTab[i]->type) { 
	      case XML_DOCUMENT_NODE:
		// if query was / we get back the whole document
		os << "NodeSet contains Document node" <<endl;
		
		xmlDocDump(stdout,doc );
		break;
	      case XML_ELEMENT_NODE:
		// otherwise use different function to dump
		xmlElemDump(stdout, doc, query_result->nodesetval->nodeTab[i]);
		break;
	      default:
		// We may get back other nodes, but I don't care about them.
		os << "Much Weirdness" << endl;
	      }
	    }
	  }
	}
	fflush(stdout);
	
	os << endl;
	os << "======= NODE DUMP END =========" << endl;
	os.flush();     
	break;
      case XPATH_BOOLEAN:
	os << "======== QUERY RESULT IS A BOOLEAN ========" << endl;
	os << "===== Value is: " << boolalpha << (bool)((query_result->boolval > 0) ? true : false)  << endl;
	break;
      case XPATH_NUMBER:
	os << "======== QUERY RESULT IS A NUMBER =========" << endl;
	os << "===== Value is: " << (query_result->floatval) << endl;
	break;
      case XPATH_STRING:
	os << "======== QUERY RESULT IS A STRING =========" << endl;
	os << "======== Value is: " << (query_result->stringval) <<endl;
	break;
      case XPATH_POINT:
	os << "======== QUERY RESULT IS XPATH_POINT ========" <<endl;
	break;
      case XPATH_RANGE:
	os << "======== QUERY RESULT IS XPATH_RANGE ========" << endl;
	break;
      case XPATH_LOCATIONSET:
	os << "======== QUERY RESULT IS XPATH_LOCATIONSET =====" << endl;
	break;
      case XPATH_USERS:
	os << "======== QUERY RESULT IS XPATH_USERS ==========" << endl;
	break;
      case XPATH_XSLT_TREE:
	os << "======== QUERY RESULT IS XPATH_XSLT_TREE ======" << endl;
	break;
      default:
	os << "======== DONT KNOW WHAT TO DO WITH QUERY_RESULT =======" <<endl;
	break;
      }
    }
  }
  else 
  {
    // No Document -- Nothing to print.
  }
}

//------------------------------------------------
// private

#if 0
/*! This contentious piece of code goes through the whole "tree"
 *  and registers every single namespace it finds into the XPath context.
 *  This may or may not be standard compliant as is */
void 
BasicXPathReader::snarfNamespaces(xmlNodePtr current_node,
				  xmlXPathContextPtr xpath_context) 
{
      
  std::ostringstream error_message;
  xmlNsPtr *ns_ptr_array;
  xmlNsPtr *shadow;
  xmlDocPtr doc;
  int i;

  doc = docref->getDoc();

  // Get the namespaces applying to the current node
  ns_ptr_array = xmlGetNsList(doc,current_node);
    
  // Walk list and register xpath
  if(  ns_ptr_array != NULL ) {
    shadow = ns_ptr_array;
    while( *shadow != NULL ) {
      if( (*shadow)->prefix != NULL ) { 
	xmlXPathRegisterNs(xpath_context, (*shadow)->prefix, (*shadow)->href);
      }
      shadow++; // Next
    }
   
    /* xmlFreeNsList(*ns_ptr_array); */
  }
      
}
#endif



/* Templated getPrimitive function. Called by all the overloaded
   get wrappers. This level of indirection also allows that the 
   wrappers submit a string "ptype" containing their type, so that
   you can check what type fails to be matched. I may eliminated 
   this and simply have a templated get --- Dunno for now.
   -- For any given T, it does the following:
   i) calls getPrimitiveString to get the element content
   ii) dumps the content string into an istringstream
   iii) tries to read back from the istringstream into T,
   converting it if possible. Failure throws exception.
   iv)  checks that the istringstream contains no more non
   whitespaces after T, ie T is unique in the tag.
*/

template< typename T >
void
BasicXPathReader::getPrimitive(const std::string& xpath, T& result, const char* ptype) 
{
  std::string pstring;

  // Get the std::string value of the query 
  std::ostringstream error_message;
  try {
    getPrimitiveString(xpath, pstring);
  }
  catch(const std::string& e) { 
    throw e;
  }
	
  // Make an istringstream from the std::string value 
  std::istringstream is(pstring);
  std::ios::iostate state;
	
  // Save the exception mask of the istringstream 
  state = is.exceptions();
	
  // Force the istringstream to throw an exception if it fails 
  // (to convert the string to our type)
  is.exceptions(std::ios::failbit);
  is.setf(ios_base::boolalpha);

  // Try to read the type from the istringstream. 
  //   bool-s should be "true" or "false" stirngs 
  try { 
    is  >> result;
  }
  catch(std::istringstream::failure e) { 
    error_message << "XPath Query: " << xpath 
		  << ": Error: Failed to convert query result string: " 
		  << pstring << " to " << ptype;
    throw error_message.str();
  }
	
  // Turn off exceptions on failure
  is.exceptions(state);
	
  // Look for non whitespaces in whats left in the stream 
  // (ie in what's left from the string
  if( findNonWhitespace(is) ) {
    error_message << "XPath Query: " << xpath << ": Error: query result "  << pstring << " does not contain a single " << ptype;
    throw error_message.str();
  }
}
    

/* Templated getAttribute function, much like getPrimitive<T>. 
 * apart from the case of string, all the overloaded getAttribute
 * functions wrap this. It does the following:
 * 
 *    i) call getAttributeString to match the attribute
 *   ii) converts to desired type in the same way as getPrimitive
 */
template< typename T >
void
BasicXPathReader::getPrimitiveAttribute(const std::string& xpath_to_node,
					const std::string& attrib_name,
					T& result,
					const char* ptype)
{
  std::string stringval;
  std::ostringstream error_message;
	
  // Get the string value of the attribute
  try {
    getAttributeString(xpath_to_node, attrib_name, stringval);
  }
  catch (const std::string& e) {
    throw e;
  }
	
  // Put it into an istringstream
  std::istringstream is(stringval);
	
  // Save the default exception mask (throw no exceptions)
  std::ios::iostate state;
  state = is.exceptions();
	
  // make the istringstream throw a wobbly on conversion failures 
  is.exceptions(std::ios::failbit);
	
  // Do the conversion, booleans must be "true" or "false"
  try { 
    is >> std::boolalpha >> result;
  }
  catch(std::istringstream::failure e) { 
    error_message << "Failed to convert attribute string: " << stringval
		  << " to type " << ptype << " for attribute " << attrib_name
		  << " on node identified by XPath: " << xpath_to_node;
    throw error_message.str();
  }
	
  // Restore excpetion mask (no exceptions)
  is.exceptions(state);
	
  // look for non-whitespaces following.
  if( findNonWhitespace(is) ) {
    error_message << "Attribute string: " << stringval
		  << " for attribute " << attrib_name
		  << " on node identified by XPath: "
		  << " doesn't contain a single " << ptype;
    throw error_message.str();
  }
}


void BasicXPathReader::registerNamespace(const std::string& prefix, const std::string& uri) 
{
  if( docref != 0x0 ) 
  {
    if( xpath_context != 0x0 ) 
    {
      xmlXPathRegisterNs(xpath_context, (xmlChar *)prefix.c_str(),
			 (xmlChar *)uri.c_str());
      nslist.push_back(NameSpace(prefix, uri));
    }
  }
}


/* Ensure that the query returned something non-null */
void BasicXPathReader::checkQuery(const std::string& xpath) 
{
  std::ostringstream error_message;
  // Check that the result of the XPath is a node set
  // Note: don't have to worry about query_result being NULL,
  // that was already checked in evaluate XPath
  if( query_result != NULL ) { 
	  
    if( query_result->type != XPATH_NODESET ) {
      error_message << "XPath Query: " << xpath << " didn't return a node_set" 
		    << std::endl;
      throw error_message.str();
    }
	  
    // Check that the nodesetval is not NULL
    if( query_result->nodesetval == 0x0 ) { 
      error_message << "XPath Query: " << xpath << " returned a NULL node set"
		    << std::endl;
      throw error_message.str();
    }
  }
  else { 
    throw "query_result is NULL in checkQuery";
  }
      
}
 

/* Ensure that the query has returned something non-null, 
   and that the query result consists of a unique node with
   simple content */
void BasicXPathReader::checkQueryPrimitive(const std::string& xpath) 
{
  std::ostringstream error_message;
      
  // first check that the query result is a nonnempty node set
  try { 
    checkQuery(xpath);
  } catch(const std::string& e) {
    throw e;
  }
      
  // Check that the node set contains only 1 element
  if( query_result->nodesetval->nodeNr != 1 ) {
    error_message << "XPath Query: " << xpath << " did not return unique node."
		  << " nodes returned = " << query_result->nodesetval->nodeNr 
		  << std::endl;
    throw error_message.str();
  }
      
  // Check that the node returned is an element 
  xmlNodePtr res_node = query_result->nodesetval->nodeTab[0];
      
  if ( res_node->type != XML_ELEMENT_NODE ) {
    error_message << "XPath Query: " << xpath << " returned a non-element node"
		  << std::endl;
    throw error_message.str();
  }
      
  // Check that the element returned is simple (ie contains no more elements)
  if( res_node->children != NULL) { 
    for( res_node=res_node->children; res_node != NULL; 
	 res_node=res_node->next ) {
      if( res_node->type == XML_ELEMENT_NODE ) { 
	error_message << "XPath Query: " << xpath << " returned non-simple node" << std::endl;
	throw error_message.str();
      }
    }
  }
}
    

/* Get a std::string from a path expression that matches checkQueryPrimitive.
 * this std::string is copied into result, after which the query result
 * is freed (hopefully properly).
 */
void BasicXPathReader::getPrimitiveString(const std::string& xpath, std::string& result)
{
  std::ostringstream error_message;

  // Run the query 
  try { 
    evaluateXPath(xpath);
  }
  catch(const std::string& e) {
    std::cerr << e;
    throw e;
  }
      
  try { 
    checkQueryPrimitive(xpath);
  }
  catch(const std::string& e) {
    throw e;
  }
      
      
  // get the "simple" content as a std::string.
  xmlChar *element_content = xmlNodeGetContent(query_result->nodesetval->nodeTab[0]);
      
  // check it for non-nullness */
  if(element_content != NULL ) {
    // COPY it into result 
    result = (std::string)((const char *)element_content);
    // Free memory malloced by libxml2
    xmlFree(element_content); 
  }
  else {
    // otherwise the query ran, but the node returned is empty
    result = "";
  }
      
  // Clean up post query
  if ( query_result != NULL ) { 
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }
}


/* Get a std::string from a path expression that matches checkQueryPrimitive.
 * this std::string is copied into result, after which the query result
 * is freed (hopefully properly).
 */
void BasicXPathReader::setPrimitiveString(const std::string& xpath, const std::string& to_write){
  std::ostringstream error_message;

  // Run the query 
  try { 
    evaluateXPath(xpath);
  }
  catch(const std::string& e) {
    std::cerr << e;
    throw e;
  }
      
  try { 
    checkQueryPrimitive(xpath);
  }
  catch(const std::string& e) {
    throw e;
  }
      
  xmlChar *new_content = xmlCharStrndup(to_write.c_str(), to_write.length());

  // get the "simple" content as a std::string.
  xmlNodeSetContent(query_result->nodesetval->nodeTab[0],new_content);
  
    
  xmlFree(new_content);

  // Clean up post query
  if ( query_result != NULL ) { 
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }
}

    
/* Get the attribute std::string from a query, that satisfies 
   checkQuery, and is a unique element node ( but doesn't have to have 
   simple content). It copies the attribute into result after which 
   the query result is freed (hopefully properly).
*/
void BasicXPathReader::getAttributeString(const std::string& xpath_to_node,
					  const std::string& attrib_name,
					  std::string& result) 
{
  std::ostringstream error_message;

  // First check that the xpath to node identifies a unique node
      
  // run the query 
  // This is the bit that checks that we have a document to query 
  try {
    evaluateXPath(xpath_to_node);
  }
  catch (const std::string& e) {
    throw e;
  }
      
  // check for general non nullness
  try { 
    checkQuery(xpath_to_node);
  }
  catch (const std::string& e) {
    throw e;
  }
      
  // check node for uniqueness
  if( query_result->nodesetval->nodeNr != 1 ) { 
    error_message << "XPath Query: " << xpath_to_node
		  << " does not identify a unique node" ;
    throw error_message.str();
  }
      
  // OK, so at this point a unique node is identified by xpath_to_node
  //
  // construct query to get the attribute
  std::string query = "string(" + xpath_to_node + "/@" + attrib_name + ")";
      
  // run the query 
  try { 
    evaluateXPath(query);
  }
  catch(const std::string& e) {
    error_message << "Failed to find attribute: " << attrib_name 
		  << " on node selected by XPath Query: " << xpath_to_node;
    throw error_message.str();
  }
      
  // check that we have a result
  if( query_result == NULL ) {
    error_message << "Query for attribute: " << attrib_name 
		  << " on node selected by XPath Query: " << xpath_to_node
		  << " returned NULL query result";
    throw error_message.str();
  }
      
  // check that the result is a std::string
  if( query_result->type != XPATH_STRING ) {
    error_message << "Query for attribute: " << attrib_name
		  << " on node selected by XPath Query: " << xpath_to_node
		  << " didn't return a string value";
    throw error_message.str();
  }
      
  // if the string is null then it presumably is an empty string(?)
  if( query_result->stringval == NULL ) {
    result = "";
  }
  else {
    // otherwise COPY the result string into result.
    result = std::string((const char *)(query_result->stringval));
  }
      
  // Post query cleanup
  if ( query_result != NULL ) { 
    xmlXPathFreeObject(query_result);
    query_result = NULL;
  }
}


/* Once you snarfed your data out of s, call this function to check
   there are no further whitespaces (ie another string) following 
   -- inelegant, but I don't know an elegant way without casting 
   -- doubles etc. */
bool BasicXPathReader::findNonWhitespace(std::istringstream& s) 
{
  char search;
	
  // While the stream is not empty 
  while(!s.eof()) 
  {
    //getchar 
    s.get(search);
	
    // if its not a space return true -- (ie found non-whitespace) 
    if( !isspace(search) ) { 
      return true;
    }
  }
	
  // Ok we reached eof without finding a non-whitespace.
  return false;
}
    
