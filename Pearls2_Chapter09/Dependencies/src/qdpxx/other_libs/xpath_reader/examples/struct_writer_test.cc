#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <xml_struct_writer.h>
#include <xml_array_writer.h>

using namespace std;
using namespace XMLStructWriterAPI;

int main(int argc, char *argv[]) 
{

  ofstream f_out("struct_writer_test.xml");
  XMLFileStructWriter* root=new XMLFileStructWriter(f_out, "root",0,true);
  root->writeSimple("fred", (string)"I am simple");

  XMLFileStructWriter* child=root->structChild("pooh");
  
  child->writeSimple((string)"jim", (int)5);
  child->writeSimple((string)"bob", (float)5.4);
  
  delete child;

  XMLFileArrayWriter* rootarr=root->arrayChild("Ginger", "Biscuits", SIMPLE_STRING);
  rootarr->elementSimple((string)"Biscuit The First");
  rootarr->elementSimple((string)"Biscuit the Second");
  rootarr->elementSimple((string)"Ham");
  rootarr->elementSimple((string)"Grits");

  delete rootarr;

  child = root->structChild("Freak");

  XMLFileStructWriter* childchild=child->structChild("jimmi");
  childchild->writeSimple((string)"Scooby", (int)10);
  childchild->writeSimple((string)"Shaggy", (string)"Zoiks");

  delete childchild;
  delete child;
  delete root;

  ostringstream buffer;
  XMLBufferStructWriter* broot=new XMLBufferStructWriter(buffer, "root", 0, false);
  broot->writeSimple("jim", (string)"I am a simple string      ");
  XMLBufferStructWriter* bchild=broot->structChild("fred");
  bchild->writeSimple("kappa",(float)0.1400);
  delete bchild;
  delete broot;

  cout << buffer.str();
  
  ofstream f_out2("array_writer_test.xml");
  XMLFileArrayWriter* arroot=new XMLFileArrayWriter(f_out2, "arr", "el", SIMPLE_INT, 0, false);
  
  arroot->elementSimple(4);
  arroot->elementSimple(5);
  arroot->elementSimple(6);

  try {
    arroot->elementSimple((float)4.6);
  }
  catch(const char *e) {
    cout << "Caught Deliberate error: " << e;
    cout << endl;
  }

  delete arroot;

  arroot = new XMLFileArrayWriter(f_out2, "arrarr", "jim", STRUCT, 0, false);
  XMLFileArrayWriter* jim=arroot->elementArray(SIMPLE_FLOAT);
  jim->elementSimple((float)4.0);
  jim->elementSimple((float)5.1);
  jim->elementSimple((float)6.2);
  delete jim;
  
  jim=arroot->elementArray(SIMPLE_FLOAT);
  jim->elementSimple((float)3.2);
  jim->elementSimple((float)4.3);
  jim->elementSimple((float)5.4);
  delete jim;
  delete arroot;

  arroot = new XMLFileArrayWriter(f_out2, "arrstruct", "guff", STRUCT, 0, false);
  XMLFileStructWriter *guff=arroot->elementStruct();
  guff->writeSimple("fred", (string)"Flintstone");
  guff->writeSimple("barney", (string)"Rubble");
  delete guff;

  guff=arroot->elementStruct();
  guff->writeSimple("fred", (string)"Flintstone2");
  guff->writeSimple("barney", (string)"Rubble2");
  delete guff;

  delete arroot;
  

}


