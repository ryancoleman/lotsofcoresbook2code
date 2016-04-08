#include <xml_writer.h>
#include <xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using namespace XMLWriterAPI;
using namespace XMLXPathReader;
using namespace std;

int main(int argc, char *argv[]) 
{

  XMLSimpleStringWriter swriter(true);

  // Write some text 
  swriter.openTag("jim");
  swriter.write((string)"Foobar");
  swriter.closeTag();
  
  try { 
    swriter.openTag("jim");
    swriter.write((string)"Foobar");
    swriter.closeTag(); 
  }
  catch( const string& error) { 
    cout << "Error: " << error << endl;
  }
  
  cout << "Content of document: " << endl;
  cout << swriter.str();
  cout << "Finished with swriter2" << endl;

  /* XMLSimpleFileWriter swriter2("jimbob.xml"); */
  XMLSimpleStringWriter swriter2;
  
  swriter2.openTag("foo");
  
  AttributeList alist;
  alist.push_back(Attribute("xmlns:bjns", (string)"http://www.ph.ed.ac.uk/~bj"));
  alist.push_back(Attribute("fred", (string)"jim"));
  swriter2.openTag((string)"bar", alist);

  alist.clear();
  alist.push_back(Attribute("kappa", (double)0.1355));
  
  swriter2.openTag("bjns", "funktastic");
  swriter2.emptyTag("fred", alist);
  
  swriter2.closeTag(); // funktastic
  swriter2.closeTag(); // bar
  
  Array<double> x(3);

  x[0] = 4;
  x[1] = 5;
  x[2] = 6;
  swriter2.openTag("x");
  swriter2.write(x); 
  swriter2.closeTag(); // x 
  
  
 
  Array< TComplex<string> > sa(2);

  sa[0] = TComplex<string>("re0", "im0");
  sa[1] = TComplex<string>("re1", "im1");

  
  swriter2.openTag("stringArray");
  swriter2.write("length", "fred", "idx", 5, sa);
  swriter2.closeTag();
  

  swriter2.openTag("unnormalizedString");
  swriter2.write((string)"        Spaces before me spaces behind me      ");
  swriter2.closeTag();

  swriter2.closeTag();

  cout << "Swriter2 Contains:" <<endl;
  cout << swriter2.str() << endl;

  istringstream foo(swriter2.str());

  XPathReader reader;
  try {
    reader.open(foo);
    
  }
  catch( const string& error) { 
    cout << "Reader Error: " << error << endl;
  }


  Array<double> y;
  reader.getXPath("/foo/x", y);

  cout << "Read back array into y" << endl;
  cout << "y.size = " << y.size() << endl;
  int i;
  int size=y.size();
  for(i=0; i < size; i++) {
    cout << " y [ " << i << " ] = " << y[i] << endl;
  }
  cout << endl;

  Array< TComplex <string> > back;
  reader.getXPath("/foo/stringArray", back);

  cout << "Read back array into back" << endl;
  cout << "back.size = " << back.size() << endl;
  
  size=back.size();
  for(i=0; i < size; i++) {
    cout << " back [ " << i << " ] = " << "Re: " <<  back[i].real() << " " << "Im : " << back[i].imag()  << endl;
  }
  cout << endl;

  string unnormalized_str;

  reader.getXPath<string>("//unnormalizedString", unnormalized_str);

  cout << "Unnormalized string is: \"" << unnormalized_str << "\"" << endl;

  reader.close();

  return(EXIT_SUCCESS);
}

