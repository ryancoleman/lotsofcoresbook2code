// $Id: t_qio.cc,v 1.29 2007-06-10 21:43:52 edwards Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("t_qio.xml");
  push(xml_out, "t_qio");

  QDP_serialparallel_t serpar = QDPIO_SERIAL;
  std::string test_file;

  // First test binary buffer writer in SINGLEFILE
  {
    test_file = "t_qio_buffer.lime";

    XMLBufferWriter file_xml, record_xml;
    Double d = 17.123456789010234567890123456;
    int rob = -5;

    {
      push(file_xml,"file_fred");
      write(file_xml,"d", d);
      push(file_xml,"file_sally");
      write(file_xml,"rob", rob);
      pop(file_xml);
      pop(file_xml);
    }
    
    {
      push(record_xml,"record_fred");
      write(record_xml,"d", d);
      push(record_xml,"record_sally");
      write(record_xml,"rob", rob);
      pop(record_xml);
      pop(record_xml);
    }

    {
      BinaryBufferWriter record_bin;
      write(record_bin, d);
      record_bin << rob;
      QDPUtil::n_uint32_t check = record_bin.getChecksum();
      record_bin << check;

      QDPFileWriter to_file(file_xml,test_file,QDPIO_SINGLEFILE,serpar,QDPIO_OPEN);
      write(to_file, record_xml, record_bin);

      QDPIO::cout << "BufferWriter test: d= " << d << "  rob=" << rob << "  checkSum= " << check << std::endl;
        }
  }


  // Read binary with a BufferReader
  {
    test_file = "t_qio_buffer.lime";

    XMLReader file_xml, record_xml;
    BinaryBufferReader record_bin;
    Double d;
    int rob;

    QDPFileReader from_file(file_xml,test_file,serpar);
    read(from_file, record_xml, record_bin);

    read(record_bin, d);
    record_bin >> rob;
    QDPUtil::n_uint32_t check;
    record_bin >> check;
    QDPIO::cout << "BufferReader test: d= " << d << "  rob=" << rob << "  checkSum= " << check << std::endl;
  }


  //
  // Loop over IO types for lattice fields
  //
  for(int i=0; i < 3; ++i)
  {
    QDP_volfmt_t volfmt;

    switch (i)
    {
    case 0:
      volfmt = QDPIO_SINGLEFILE;
      test_file = "t_qio_single.lime";
      QDPIO::cout << "\n\n\n\n***************SINGLEFILE tests*************\n" << std::endl;
      push(xml_out, "Singlefile");
      break;

    case 1:
      volfmt = QDPIO_MULTIFILE; 
      test_file = "t_qio_multi.lime";
      QDPIO::cout << "\n\n***************MULTIFILE tests*************\n" << std::endl;
      push(xml_out, "Multifile");
      break;

    case 2:
      volfmt = QDPIO_PARTFILE; 
      test_file = "t_qio_part.lime";
      QDPIO::cout << "\n\n***************PARTFILE tests*************\n" << std::endl;
      push(xml_out, "Partfile");
      break;

    default:
      QDPIO::cerr << "something wrong here" << std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "\n\n***************TEST WRITING*************\n" << std::endl;

    {
      push(xml_out, "Writing");

      XMLBufferWriter file_xml, record_xml;
      Double d = 17;
      int rob = -5;

      {
	push(file_xml,"file_fred");
	write(file_xml,"d", d);
	push(file_xml,"file_sally");
	write(file_xml,"rob", rob);
	pop(file_xml);
	pop(file_xml);
      }

      {
	push(record_xml,"record_fred");
	write(record_xml,"d", d);
	push(record_xml,"record_sally");
	write(record_xml,"rob", rob);
	pop(record_xml);
	pop(record_xml);
      }

      QDPFileWriter *to_ptr;
      to_ptr = new QDPFileWriter(file_xml,test_file,volfmt,serpar,QDPIO_OPEN);
      QDPFileWriter& to = *to_ptr;

      write(xml_out, "file_xml", file_xml);
      write(xml_out, "open_to.bad", to.bad());

#if 0
      if (volfmt == QDPIO_SINGLEFILE)   // there seems to be a bug in multifile here
      {
        multi1d<Complex> ff(5);
        Double fsum = 0;
        for(int i=0; i < ff.size(); ++i)
        {
    	  random(ff[i]);
	  fsum += norm2(ff[i]);
        }
        write(to,record_xml,ff);
	write(xml_out, "first_to.bad", to.bad());

        QDPIO::cout << "First record test: fsum=" << fsum << std::endl;
        push(xml_out, "Record1");
        write(xml_out, "record_xml", record_xml);
        write(xml_out, "fsum", fsum);
        write(xml_out, "ff", ff);
        pop(xml_out);
      }
#endif

      {
	LatticeComplex a;
	random(a);
	write(to,record_xml,a);
	write(xml_out, "second_to.bad", to.bad());

	Real atest = Real(innerProductReal(a,shift(a,FORWARD,0)));
	QDPIO::cout << "Second record test: innerProduct(a,shift(a,0))=" 
		    << atest << std::endl;
	push(xml_out, "Record2");
	write(xml_out, "record_xml", record_xml);
	write(xml_out, "atest", atest);
	pop(xml_out);
      }

#if 1
      {
	LatticeColorMatrix b;
	random(b);
	write(to,record_xml,b);
	write(xml_out, "third_to.bad", to.bad());
	
	Real btest = Real(innerProductReal(b,shift(b,FORWARD,0)));
	QDPIO::cout << "Third record test: innerProduct(b,shift(b,0))=" 
		    << btest << std::endl;
	push(xml_out, "Record3");
	write(xml_out, "record_xml", record_xml);
	write(xml_out, "btest", btest);
	pop(xml_out);
      }
#endif

      {
	multi1d<LatticeColorMatrix> c(Nd);
        Double fsum = 0;
        for(int i=0; i < c.size(); ++i)
        {
    	  random(c[i]);
	  fsum += innerProductReal(c[i],shift(c[i],FORWARD,0));
        }

	QDPIO::cout << "Fourth record test is a lattice Color Matrix" << std::endl;
	StopWatch swatch;
	swatch.reset();
	swatch.start();
        write(to,record_xml,c);
	swatch.stop();
	double secs = swatch.getTimeInSeconds();
	double bytes = Layout::vol()*Nd*Nc*Nc*2*sizeof(REAL);
	double Mbytes = bytes/(double)(1024*1024);
	QDPIO::cout << "Wrote " << Mbytes << " Mbytes in " << secs << " seconds " << std::endl;
	QDPIO::cout << "Transfer Rate: " << Mbytes/secs << " Mbytes/s" << std::endl;
	QDPIO::cout << std::flush;

	write(xml_out, "fourth_to.bad", to.bad());

        QDPIO::cout << "Fourth record test: fsum=" << fsum << std::endl;
        push(xml_out, "Record4");
        write(xml_out, "record_xml", record_xml);
        write(xml_out, "fsum", fsum);
//        write(xml_out, "c", c);
        pop(xml_out);
      }

      close(to);
      delete to_ptr;

      pop(xml_out);   // writing
    }

    QDPIO::cout << "\n\n***************TEST READING*******************\n" << std::endl;

    {
      push(xml_out, "Reading");

      XMLReader file_xml, record_xml;
      QDPFileReader* from_ptr;
      from_ptr = new QDPFileReader(file_xml,test_file,serpar);
      QDPFileReader& from = *from_ptr;

      QDPIO::cout << "Here is the contents of  file_xml" << std::endl;
      std::ostringstream outbuf;
      file_xml.print(outbuf);
      QDPIO::cout << outbuf.str() << std::endl;
      

      write(xml_out, "file_xml", file_xml);
      write(xml_out, "open_frombad", from.bad());

#if 0
      if (volfmt == QDPIO_SINGLEFILE)
      {
        multi1d<Complex> ff(5);   // needs to be free
        read(from,record_xml,ff);
	write(xml_out, "first_from.bad", from.bad());

        QDPIO::cout << "Here is the contents of first  record_xml" << std::endl;
	//        record_xml.print(cout);
	std::ostringstream ost;
	record_xml.print(ost);
	QDPIO::cout << ost.str() << std::endl;

        Double fsum = 0;
        for(int i=0; i < ff.size(); ++i)
  	  fsum += norm2(ff[i]);
        QDPIO::cout << "First record test: fsum=" << fsum << std::endl;
        push(xml_out, "Record1");
        write(xml_out, "record_xml", record_xml);
        write(xml_out, "fsum", fsum);
        write(xml_out, "ff", ff);
        pop(xml_out);
      }
#endif

      {
	LatticeComplex a;
	read(from,record_xml,a);
	write(xml_out, "second_from.bad", from.bad());

	QDPIO::cout << "Here is the contents of second  record_xml" << std::endl;
	std::ostringstream ost;
	record_xml.print(ost);
	QDPIO::cout << ost.str() << std::endl;


	Real atest = Real(innerProductReal(a,shift(a,FORWARD,0)));
	QDPIO::cout << "Second record check: innerProduct(a,shift(a,0))=" 
		    << atest << std::endl;
	push(xml_out, "Record2");
	write(xml_out, "record_xml", record_xml);
	write(xml_out, "atest", atest);
	pop(xml_out);
      }

#if 1
      {
	LatticeColorMatrix b;
	QDPIO::cout << "LatticeColorMatrix created " << std::endl;
	QDPIO::cout << "About to read from it " << std::endl;

	read(from,record_xml,b);
	write(xml_out, "third_from.bad", from.bad());

	QDPIO::cout << "Here is the contents of third  record_xml" << std::endl;
	std::ostringstream ost;
	record_xml.print(ost);
	QDPIO::cout << ost.str() << std::endl;


	Real btest = Real(innerProductReal(b,shift(b,FORWARD,0)));
	QDPIO::cout << "Third record check: innerProduct(b,shift(b,0))=" 
		    << btest << std::endl;
	push(xml_out, "Record3");
	write(xml_out, "record_xml", record_xml);
	write(xml_out, "btest", btest);
	pop(xml_out);
      }
#endif

      {
	multi1d<LatticeColorMatrix> c(Nd);   // array size should be free

	QDPIO::cout << "Multi1d LatticeColorMatrix created " << std::endl;
	QDPIO::cout << "About to try treading into it" << std::endl;
	StopWatch swatch;
	swatch.reset();
	swatch.start();
        read(from,record_xml,c);
	swatch.stop();

	double secs = swatch.getTimeInSeconds();
	double bytes = Layout::vol()*Nd*Nc*Nc*2*sizeof(REAL);
	double Mbytes = bytes/(double)(1024*1024);
	QDPIO::cout << "Read " << Mbytes << " Mbytes in " << secs << " seconds " << std::endl;
	QDPIO::cout << "Transfer Rate: " << Mbytes/secs << " Mbytes/s" << std::endl;
	
	write(xml_out, "fourth_from.bad", from.bad());

        Double fsum = 0;
        for(int i=0; i < c.size(); ++i)
	  fsum += innerProductReal(c[i],shift(c[i],FORWARD,0));

        QDPIO::cout << "Fourth record test: fsum=" << fsum << std::endl;
        push(xml_out, "Record4");
        write(xml_out, "record_xml", record_xml);
        write(xml_out, "fsum", fsum);
//        write(xml_out, "c", c);
        pop(xml_out);
      }

      close(from);   // reading
      //      delete from_ptr; // delete it to keep valgrind happy.
      pop(xml_out);
    }

    pop(xml_out);  // single or multifile
  }


  pop(xml_out);   // t_qio
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
