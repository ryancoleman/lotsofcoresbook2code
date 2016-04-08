/* ID: $Id: basic_xpath_reader_test3.cc,v 1.1 2004-04-27 11:22:59 bjoo Exp $
 *
 * file: basic_xpath_reader2.cc
 *
 * a second example of using BasicXPathReader -- now to bind 
 * primitive types 
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

  // Get the reader
  BasicXPathReader reader;

  // Open the file
  try { 
    reader.open("edwards.xml");
  } 
  catch (string &error_mesg) {
    cerr << error_mesg << endl;
    throw;
  }
  cout << "Document Open Complete" << endl;

  cout << "Creating reader from path /foo/bar" << endl;
  BasicXPathReader reader2(reader, "/foo/bar");

  cout << "Printing Document " << endl;
  reader2.printDoc(cout);
  cout << endl;

  cout << "Printing Root " << endl;
  reader2.printRoot(cout);
  cout << endl;
  cout << endl;

  cout << "Printing self" << endl;
  reader2.print(cout);
  cout << endl;

  cout << "Printing Children " << endl;
  reader2.printChildren(cout);
  cout << endl;

  return EXIT_SUCCESS;
}

