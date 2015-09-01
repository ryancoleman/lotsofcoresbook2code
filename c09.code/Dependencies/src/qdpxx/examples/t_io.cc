// t_io.cc

#include "qdp.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  LatticeReal a;
  Double d = 17;
  std::string astring = "hello world";
  random(a);

  // Test writing
  {
    BinaryFileWriter tobinary("t_io.bin");
    write(tobinary, a);
    write(tobinary, d);
    write(tobinary, astring);
    QDPIO::cout <<  "WriteBinary: t_io.bin:checksum = " << tobinary.getChecksum() << std::endl;
    tobinary.close();
  }

  LatticeReal aa;
  Double dd = 0.0;
  random(aa);

  // Test reading
  {
    BinaryFileReader frombinary("t_io.bin");
    read(frombinary, aa);
    read(frombinary, dd);
    read(frombinary, astring, 100);
    QDPIO::cout <<  "ReadBinary: t_io.bin:checksum = " << frombinary.getChecksum() << std::endl;
    frombinary.close();
  }

  // Test seeks
  {
    BinaryFileReader frombinary("t_io.bin");
    BinaryFileReader::pos_type pos = frombinary.currentPosition();
    QDPIO::cout <<  "ReadBinary: t_io.bin: position = " << pos << std::endl;
    read(frombinary, aa); 
    QDPIO::cout <<  "ReadBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl; 
    pos = frombinary.currentPosition();
    read(frombinary, dd);
    QDPIO::cout <<  "ReadBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl;
    read(frombinary, astring, 100);
    QDPIO::cout <<  "ReadBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl;
    QDPIO::cout <<  "ReadBinary: t_io.bin:checksum = " << frombinary.getChecksum() << std::endl;

    QDPIO::cout <<  "ReadBinary: now seek: current position = " << frombinary.currentPosition() << std::endl;
    frombinary.seek(pos);
    QDPIO::cout <<  "ReadBinary: after seeking: new current position = " << frombinary.currentPosition() << std::endl;
    read(frombinary, dd);
    QDPIO::cout <<  "ReadBinary: after reading: position = " << frombinary.currentPosition() << std::endl;
    pos = frombinary.currentPosition();
    frombinary.close();
  }

  // Test reader/writer
  {
    BinaryFileReaderWriter frombinary("t_io.bin");
    BinaryFileReaderWriter::pos_type pos = frombinary.currentPosition();
    QDPIO::cout <<  "ReadWriteBinary: t_io.bin: position = " << pos << std::endl;
    read(frombinary, aa); 
    QDPIO::cout <<  "ReadWriteBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl; 
    pos = frombinary.currentPosition();
    read(frombinary, dd);
    QDPIO::cout <<  "ReadWriteBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl;
    read(frombinary, astring, 100);
    QDPIO::cout <<  "ReadWriteBinary: t_io.bin: position = " << frombinary.currentPosition() << std::endl;
    QDPIO::cout <<  "ReadWriteBinary: t_io.bin:checksum = " << frombinary.getChecksum() << std::endl;

    QDPIO::cout <<  "ReadWriteBinary: now seek: current position = " << frombinary.currentPosition() << std::endl;
    frombinary.seek(pos);
    QDPIO::cout <<  "ReadWriteBinary: after seeking: new current position = " << frombinary.currentPosition() << std::endl;
    read(frombinary, dd);
    QDPIO::cout <<  "ReadWriteBinary: after reading: position = " << frombinary.currentPosition() << std::endl;
    pos = frombinary.currentPosition();
    frombinary.close();
  }

  {
    XMLFileWriter toxml("t_io.xml");
    push(toxml,"t_io");
    write(toxml,"a",a);
    write(toxml,"aa",aa);
    pop(toxml);
    toxml.flush();
    toxml.close();
  }

  Real x = 42.1;
  {
    QDPIO::cout << "Write some data to file t_io.txt\n";
    TextFileWriter totext("t_io.txt");
    totext << x;
    totext.flush();
    totext.close();
  }

  x = -1;
  {
    QDPIO::cout << "Read some data from file t_io.txt\n";
    TextFileReader fromtext("t_io.txt");
    fromtext >> x;
    fromtext.close();
  }

  QDPIO::cout << "The value :" << x << ": was read from t_io.txt" << std::endl;

  x = -1;
  {
    QDPIO::cout << "Enter a float for a test of reading stdin" << std::endl;
    QDPIO::cin >> x;
    QDP_info("The value :%g: was read from stdin", toFloat(x));
  }

  x = 17.1;
  {
    QDPIO::cout << "TextBufferWriter and Reader test: original float = " << x << std::endl;
    TextBufferWriter tobuftext;
    tobuftext << x;
    TextBufferReader frombuftext(tobuftext.str());
    frombuftext >> x;
    QDPIO::cout << "Read back a float from a TextBufferReader = " << x << std::endl;
  }

  x = 11.1;
  {
    QDPIO::cout << "BinaryBufferWriter and Reader test: original float = " << x << std::endl;
    BinaryBufferWriter tobufbinary;
    write(tobufbinary, x);

    BinaryBufferReader frombufbinary(tobufbinary.str());
    read(frombufbinary, x);
    QDPIO::cout << "Read back a float from a BinaryBufferReader = " << x << std::endl;
  }

  // Time to bolt
  QDP_finalize();

  exit(0);
}
