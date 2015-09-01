/* ID: $Id: basic_xpath_reader_test.cc,v 1.6 2003-09-10 12:07:29 bjoo Exp $ 
 * File: basic_xpath_reader_test.cc 
 * 
 * This file was a test to see how well I understood libxml xpath, 
 * and allowed me to test evaluateXPath etc. 
 *
 * essentially it is an "XPath Shell" -- you can type in XPath
 * and see what comes back 
 */

#include <basic_xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
using namespace XMLXPathReader;

int main(int argc, char *argv[])
{

  BasicXPathReader reader;

  // Check command line arguments
  if(argc != 2 ) { 
    cerr << "Usage: " << argv[0] << " <xml_file>" << endl;
    return EXIT_FAILURE;
  }

  // Read the XML
  string filename(argv[1]);
  cout << "Reading via " << filename << " via the file interface... " <<endl;

  
  try { 
    reader.open(filename);
  } 
  catch (string &error_mesg) {
    cerr << error_mesg << endl;
    throw;
  }
  
  cout << "Reader Open Complete" << endl;

  // Loop until "exit" or "quit"
  // is entered
  bool finished = false;
  while(!finished) { 
    string xpath_expression;
    cout << "Enter XPATH expression (exit or quit to finish)" << endl;
    cin >> xpath_expression;
    
    if ( (xpath_expression != "exit") && (xpath_expression != "quit") ) {
      if ( xpath_expression == "#setns" )  {
        cout << "Enter prefix" <<endl;
        string prefix;
        cin >> prefix;
        cout << "Enter URI " << endl;
        string uri;
        cin >> uri;
	
        reader.registerNamespace(prefix, uri);
      }
      else {
	try { 
	  reader.evaluateXPath(xpath_expression);
	}
	catch(string &error_mesg) { 
	  // Catch errors but don't bomb. 
	  // let the user try again...
	  cerr << error_mesg;
	}      
      } 
    }
    else {
      // exit or quit was entered 
      finished=true;
    }
    
    // If no error occurred, and also we didn't quit then print the query
    // results.
    reader.printQueryResult(cout);
  }
  
  // bye
  cout << "Closing reader." << endl;
  reader.close();
  
  return EXIT_SUCCESS;
}

