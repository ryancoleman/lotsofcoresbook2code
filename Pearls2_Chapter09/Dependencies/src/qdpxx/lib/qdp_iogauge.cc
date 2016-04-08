//
// QDP data parallel interface
/*!
 * @file
 * @brief  Various gauge readers/writers and propagator readers/writers.
 */

#include "qdp.h"
#include "qdp_iogauge.h"

#include "time.h"

#include <unistd.h>
#include <string>
using std::string;

// QCDOC HACK. QCDOC Does not have gethostname
// provided in qdp_util.cc
#ifndef HAVE_GETHOSTNAME
extern int gethostname(char *, size_t);
#endif

namespace QDP {


// Anonymous namespace
namespace
{
  // Float tolerance
  const Double tol = 1.e-5;  /* tolerance for floating point checks */

  // Grrh, I do not want to expose the plaquette code.
  void mesplq(Double& w_plaq, Double& link, const multi1d<LatticeColorMatrix>& u)
  {
    w_plaq = link = 0.0;

    // Compute the average plaquettes
    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
	/* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	LatticeColorMatrix tmp_0 = shift(u[nu],FORWARD,mu) * adj(shift(u[mu],FORWARD,nu));

	/* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);

	/* tmp = sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu))) */
	Double tmp = sum(real(trace(u[mu]*tmp_1)));

	w_plaq += tmp;
      }
    }
  
    // NERSC normalization
    w_plaq *= 2.0 / double(Layout::vol()*Nd*(Nd-1)*Nc);
  
    // Compute the average link
    for(int mu=0; mu < Nd; ++mu)
      link += sum(real(trace(u[mu])));

    link /= double(Layout::vol()*Nd*Nc);
  }

} // end anonymous namespace



//! Write a multi1d array
template<class T>
std::ostream& operator<<(std::ostream& s, const multi1d<T>& d)
{
  s << d[0];
  for(int i=1; i < d.size(); ++i)
    s << " " << d[i];

  return s;
}


ArchivGauge_t::ArchivGauge_t()
{
  mat_size = 2*Nc*(Nc-1);
  float_size = 4; // 32 bits
  nrow = Layout::lattSize();
  boundary.resize(Nd);
  boundary = 1;   // periodic
  sequence_number = 0;
  ensemble_label = "NERSC archive";
  creator = "QDP++";
  {
    const int namelen = 128;
    char name[namelen];
    gethostname(name, namelen);
    name[namelen-1] = '\0';
    creator_hardware = name;
  }
  checksum = 0;

  time_t now = time(NULL);
  {
    char *tmp = ctime(&now);
    int date_size = strlen(tmp);
    char *datetime = new(std::nothrow) char[date_size+1];
    if( datetime == 0x0 ) { 
      QDP_error_exit("Unable to allocate datetime in qdp_iogauge.cc\n");
    }

    strcpy(datetime,ctime(&now));

    for(int i=0; i < date_size; ++i)
      if ( datetime[i] == '\n' )
      {
	datetime[i] = '\0';
	date_size = i;
	break;
     }   

    creation_date = datetime;
    delete[] datetime;
  }
  archive_date  = creation_date;

  {
    std::ostringstream s;
    s << "X" << now;
    ensemble_id = s.str();
  }

  w_plaq = 0;   // WARNING: bogus
  link = 0;     // WARNING: bogus
}


//! Source header read
void read(XMLReader& xml, const std::string& path, ArchivGauge_t& header)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "mat_size", header.mat_size);
  read(paramtop, "float_size", header.float_size);
  read(paramtop, "nrow", header.nrow);
  read(paramtop, "boundary", header.boundary);
  read(paramtop, "w_plaq", header.w_plaq);
  read(paramtop, "link", header.link);
  read(paramtop, "ensemble_id", header.ensemble_id);
  read(paramtop, "ensemble_label", header.ensemble_label);
  read(paramtop, "creator", header.creator);
  read(paramtop, "creator_hardware", header.creator_hardware);
  read(paramtop, "creation_date", header.creation_date);
  read(paramtop, "archive_date", header.archive_date);
  read(paramtop, "sequence_number", header.sequence_number);

  // read a hex as a string and then convert
  {
    std::string chk;
    read(paramtop, "checksum", chk);
    std::istringstream s(chk);
    s >> header.checksum;
  }
}


//! Source header writer
void write(XMLWriter& xml, const std::string& path, const ArchivGauge_t& header)
{
  push(xml, path);

  write(xml, "mat_size", header.mat_size);
  write(xml, "float_size", header.float_size);
  write(xml, "nrow", header.nrow);
  write(xml, "boundary", header.boundary);
  write(xml, "w_plaq", header.w_plaq);
  write(xml, "link", header.link);
  write(xml, "ensemble_id", header.ensemble_id);
  write(xml, "ensemble_label", header.ensemble_label);
  write(xml, "creator", header.creator);
  write(xml, "creator_hardware", header.creator_hardware);
  write(xml, "creation_date", header.creation_date);
  write(xml, "archive_date", header.archive_date);
  write(xml, "sequence_number", header.sequence_number);

  // write as a hex 
  {
    std::ostringstream s;
    s.setf(std::ios_base::hex, std::ios_base::basefield);
    s << header.checksum;
    write(xml, "checksum", s.str());
  }

  pop(xml);
}


//-----------------------------------------------------------------------
// Read a QCD archive file header
//! Read a QCD (NERSC) Archive format gauge field header
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param cfg_in     binary writer object ( Modify )

 \note This can handle three-row format link matrices if the
 \c DATATYPE key has the value \c 4D_SU3_GAUGE_3x3
 or two-row format matrices if it has the value \c 4D_SU3_GAUGE

 The plaquette, link value, and checksum are read and compared against 
 computed values.
*/    

static void readArchivHeader(BinaryReader& cfg_in, ArchivGauge_t& header)
{
  if (Nd != 4)
  {
    QDPIO::cerr << "Expecting Nd == 4" << std::endl;
    QDP_abort(1);
  }

  const size_t max_line_length = 128;

  // The expected lattice size of the gauge field
  header.nrow.resize(Nd);

  /* For now, read and throw away the header */
  std::string line;

  QDPIO::cout << "Start of header" << std::endl;

  cfg_in.read(line, max_line_length);
  QDPIO::cout << line << std::endl;
  
  if (line != std::string("BEGIN_HEADER"))
  {
    QDPIO::cerr << "Missing BEGIN_HEADER" << std::endl;
    QDP_abort(1);
  }

  /* assume matrix size is 2*Nc*Nc (matrix is UNcompressed) 
     and change if we find out otherwise */
  header.mat_size=2*Nc*Nc;

  /* Begin loop on lines */
  int  lat_size_cnt = 0;

  while (1)
  {
    cfg_in.read(line, max_line_length);
    QDPIO::cout << line << std::endl;

    if (line == std::string("END_HEADER")) break;

    int itmp, dd;
    std::string::size_type off;

    // Snarf the first token
    char tokenn[max_line_length];
    if ( sscanf(line.c_str(), "%s", tokenn) != 1 ) 
    {
      QDPIO::cerr << __func__ 
		  << ": incorrectly parsed header line=XX" << line << "XX" << std::endl;
      QDP_abort(1);
    }
    std::string token = tokenn;

    // Scan for first non-space char after "="
    off = line.find('=');
    if ( off == std::string::npos )
    {
      QDPIO::cerr << __func__ 
		  << ": incorrectly parsed header line=XX" << line << "XX" << std::endl;
      QDP_abort(1);
    }
    off = line.find_first_not_of(' ', off+1);
    std::string value;
    if ( off == std::string::npos )
    {
//      QDPIO::cerr << __func__ 
//		  << ": incorrectly parsed header line=XX" << line << "XX" << std::endl;
//      QDP_abort(1);
      value = "";
    }
    else
    {
      value = line.substr(off, line.length()-off+1);
//    QDPIO::cout << "value = XX" << value << "XX" << std::endl;
    }


    // Scan for the datatype then scan for it
    if ( token == std::string("DATATYPE") )
    {
      /* Check if it is uncompressed */
      if ( value == std::string("4D_SU3_GAUGE_3x3") )
      {
	header.mat_size=18;   /* Uncompressed matrix */
	if (Nc != 3)
	{
	  QDPIO::cerr << __func__ << ": expecting Nc == 3" << std::endl;
	  QDP_abort(1);
	}
      }
      else if ( value == std::string("4D_SU3_GAUGE") )
      {
	header.mat_size=12;   /* Compressed matrix */
	if (Nc != 3)
	{
	  QDPIO::cerr << __func__ << ": expecting Nc == 3" << std::endl;
	  QDP_abort(1);
	}
      }
      else if ( value == std::string("4D_SU4_GAUGE") )
      {
	if (Nc != 4)
	{
	  QDPIO::cerr << __func__ << ": expecting Nc == 4" << std::endl;
	  QDP_abort(1);
	}
      }
      else
      {
	QDPIO::cerr << __func__ 
		    << ": unknown gauge type = XX" << value << "XX" << std::endl;
	QDP_abort(1);
      }
    }


    // Scan for the plaq
    double dtmp;
    if ( sscanf(line.c_str(), "PLAQUETTE = %lf", &dtmp) == 1 ) 
    {
      header.w_plaq = dtmp;
    }

    // Scan for the link
    if ( sscanf(line.c_str(), "LINK_TRACE = %lf", &dtmp) == 1 ) 
    {
      header.link = dtmp;
    }

    // Scan for the checksum
    unsigned long chk;
    if ( sscanf(line.c_str(), "CHECKSUM = %lx", &chk) == 1 ) 
    {
      header.checksum = chk;
    }

    // Scan for the sequence number
    if ( sscanf(line.c_str(), "SEQUENCE_NUMBER = %d", &itmp) == 1 ) 
    {
      header.sequence_number = itmp;
    }

    // Scan for the ensemble label
    if ( token == std::string("ENSEMBLE_LABEL") )
    {
      header.ensemble_label = value;
    }

    // Scan for the ensemble id
    if ( token == std::string("ENSEMBLE_ID") )
    {
      header.ensemble_id = value;
    }

    // Scan for the creator
    if ( token == std::string("CREATOR") )
    {
      header.creator = value;
    }

    // Scan for the creator machine
    if ( token == std::string("CREATOR_MACHINE") )
    {
      header.creator_hardware = value;
    }

    // Scan for the creator hardware
    if ( token == std::string("CREATOR_HARDWARE") )
    {
      header.creator_hardware = value;
    }

    // Scan for the creation date
    if ( token == std::string("CREATION_DATE") )
    {
      header.creation_date = value;
    }

    // Scan for the archive date
    if ( token == std::string("ARCHIVE_DATE") )
    {
      header.archive_date = value;
    }

    // Find the lattice size of the gauge field
    if ( sscanf(line.c_str(), "DIMENSION_%d = %d", &dd, &itmp) == 2 ) 
    {
      /* Found a lat size */
      if (dd < 1 || dd > Nd)
      {
	QDPIO::cerr << __func__ << ": dimension number out of bounds" << std::endl;
	QDP_abort(1);
      }

      header.nrow[dd-1] = itmp;
      ++lat_size_cnt;
    }
    
    // Find the boundary conditions
    if ( sscanf(line.c_str(), "BOUNDARY_%d = %d", &dd, &itmp) == 2 ) 
    {
      /* Found a lat size */
      if (dd < 1 || dd > Nd)
      {
	QDPIO::cerr << __func__ << ": dimension number out of bounds" << std::endl;
	QDP_abort(1);
      }
      
      header.boundary[dd-1] = itmp;
    }

    if( token == std::string("FLOATING_POINT") )
    {
      if( value == std::string("IEEE32BIG") || value == std::string("IEEE32") ) 
      {
	header.float_size=4;
      }
      else if( value == std::string("IEEE64BIG") )
      {
	header.float_size=8;
      }
      else 
      {
	QDPIO::cerr << __func__ 
		    << ": unknown floating point type = XX" << value << "XX" << std::endl;
	QDP_abort(1);
      }
    }
  }

  QDPIO::cout << "End of header" << std::endl;

  // Sanity check
  if (lat_size_cnt != Nd)
  {
    QDPIO::cerr << __func__ << ": did not find all the lattice sizes" << std::endl;
    QDP_abort(1);
  }

  for(int dd=0; dd < Nd; ++dd)
    if (header.nrow[dd] != Layout::lattSize()[dd])
    {
      QDPIO::cerr << __func__ << ": archive lattice size does not agree with current size" << std::endl;
      QDP_abort(1);
    }
}


//-----------------------------------------------------------------------
//! Read a NERSC Gauge Connection  Archive file
// See the corresponding  qdp_*_specific.cc files
//! Writes a NERSC Gauge Connection Archive gauge configuration file
/*!
 * \ingroup io
 An architecture-specific version of this routine is called by the generic
 readArchiv functions.

 The data is written in big-endian IEEE format to the file.
 If the host nachine is little-endian, the data is byte-swapped.

  \param cfg_in    A binary reader
  \param u          The gauge configuration 
  \param mat_size   The number of floating-point numbers per link matrix in
  the file. This should be 12 to write two-row format or 18 for three-row
  format. 
  \param float_size
  \param checksum   The 32bit parity checksum
  
  \pre The binary writer should have already opened the file, and should be
  pointing to the beginning of the binary data.
*/
void readArchiv(BinaryReader& cfg_in, multi1d<LatticeColorMatrix>& u, 
		n_uint32_t& checksum, int mat_size, int float_size);



// Read a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param file       path ( Read )
 */    
void readArchiv(ArchivGauge_t& header, multi1d<LatticeColorMatrix>& u, const std::string& file)
{
  BinaryFileReader cfg_in(file);

  readArchivHeader(cfg_in, header);   // read header
  n_uint32_t checksum;
  // expects to be positioned at the beginning of the binary payload
  readArchiv(cfg_in, u, checksum, header.mat_size, header.float_size);

  if (checksum != header.checksum)
  {
    QDPIO::cerr << __func__ << ": checksum mismatch: new=" << checksum 
		<< "  header value= " << header.checksum << std::endl;
    QDP_abort(1);
  }

  Double w_plaq, link;
  mesplq(w_plaq, link, u);
  if (toBool(fabs(header.w_plaq - w_plaq) > tol))
  {
    QDPIO::cerr << __func__ << ": plaquette out of bounds: new=" << w_plaq 
		<< "  header value= " << header.w_plaq << std::endl;
    QDP_abort(1);
  }

  if (toBool(fabs(header.link - link) > tol))
  {
    QDPIO::cerr << __func__ << ": link out of bounds: new=" << link 
		<< "  header value= " << header.link << std::endl;
    QDP_abort(1);
  }

  cfg_in.close();
}


//-----------------------------------------------------------------------
// Read a Archive configuration file
/*
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readArchiv(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const std::string& cfg_file)
{
  ArchivGauge_t header;

  // Read the config and its binary header
  readArchiv(header, u, cfg_file);

  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  write(xml_buf, "NERSC", header);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const std::string& e)
  { 
    QDPIO:: cerr << "Error in readArchiv: " << e << std::endl;
    QDP_abort(1);
  }
}



//-----------------------------------------------------------------------
// Read a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    
void readArchiv(multi1d<LatticeColorMatrix>& u, const std::string& cfg_file)
{
  ArchivGauge_t header;
  readArchiv(header, u, cfg_file); // throw away the header
}




//-----------------------------------------------------------------------
// Write a QCD archive file
//! Write a QCD (NERSC) Archive format gauge field
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Read )
 * \param cfg_out    binary writer object ( Modify )

 \pre The information in the header should be filled in.
 
 \note The value 0 is written as checksum.
 \note The token \c FLOATING_POINT is always given the value \c IEEE32BIG
 */    
static void writeArchivHeader(BinaryWriter& cfg_out, const ArchivGauge_t& header)
{
  if (Nd != 4)
  {
    QDPIO::cerr << "Expecting Nd == 4" << std::endl;
    QDP_abort(1);
  }

  if (Nc != 3)
  {
    QDPIO::cerr << "Expecting Nc == 3" << std::endl;
    QDP_abort(1);
  }

  std::ostringstream head;

  head << "BEGIN_HEADER\n";

  head << "CHECKSUM = ";
  head.setf(std::ios_base::hex, std::ios_base::basefield);
  head << header.checksum << std::endl;
  head.setf(std::ios_base::dec, std::ios_base::basefield);
  head << "LINK_TRACE = " << header.link << "\n"
       << "PLAQUETTE = " << header.w_plaq << "\n";

  head << "DATATYPE = 4D_SU3_GAUGE\n"
       << "HDR_VERSION = 1.0\n"
       << "STORAGE_FORMAT = 1.0\n";

  for(int i=1; i <= Nd; ++i)
    head << "DIMENSION_" << i << " = " << Layout::lattSize()[i-1] << "\n";

  for(int i=0; i < Nd; ++i)
    if (header.boundary[i] == 1)
      head << "BOUNDARY_" << (i+1) << " = PERIODIC\n";
    else if (header.boundary[i] == -1)
      head << "BOUNDARY_" << (i+1) << " = ANTIPERIODIC\n";
    else
    {
      QDPIO::cerr << "writeArchiv: unknown boundary type";
      QDP_abort(1);
    }

  head << "ENSEMBLE_ID = " << header.ensemble_id << "\n"
       << "ENSEMBLE_LABEL = " << header.ensemble_label << "\n"
       << "SEQUENCE_NUMBER = " << header.sequence_number << "\n"
       << "CREATOR = " << header.creator << "\n"
       << "CREATOR_HARDWARE = " << header.creator_hardware << "\n"
       << "CREATION_DATE = " << header.creation_date << "\n"
       << "ARCHIVE_DATE = " << header.archive_date << "\n"
       << "FLOATING_POINT = IEEE32BIG\n";

  head << "END_HEADER\n";

  cfg_out.writeArray(head.str().c_str(), 1, head.str().size());
}


// Write a QCD archive file
// See the corresponding  qdp_*_specific.cc files

//! Writes a NERSC Gauge Connection Archive gauge configuration file
/*!
 * \ingroup io
 An architecture-specific version of this routine is called by the generic
 readArchiv functions.

 The data is written in big-endian IEEE format to the file.
 If the host nachine is little-endian, the data is byte-swapped.

  \param cfg_out    A binary writer
  \param u          The gauge configuration 
  \param mat_size   The number of floating-point numbers per link matrix to
  write. This should be 12 to write two-row format or 18 for three-row format.

  \pre The binary writer should have already opened the file.
*/
  
void writeArchiv(BinaryWriter& cfg_out, const multi1d<LatticeColorMatrix>& u,
		 int mat_size);


//-----------------------------------------------------------------------
// Write a QCD archive file
// Write a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param xml        xml writer holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param file       path ( Read )
 */    
void writeArchiv(ArchivGauge_t& header, const multi1d<LatticeColorMatrix>& u, const std::string& file)
{
  Double w_plaq, link;
  mesplq(w_plaq, link, u);
  header.w_plaq = w_plaq;
  header.link = link;
  header.checksum = computeChecksum(u, header.mat_size);

  BinaryFileWriter cfg_out(file);

  writeArchivHeader(cfg_out, header);   // write header
  writeArchiv(cfg_out, u, header.mat_size);  // continuing writing after header

  cfg_out.close();
}


// Write a Archive configuration file
/*
 * \ingroup io
 *
 * \param xml        xml writer holding config info ( Read )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeArchiv(XMLBufferWriter& xml, const multi1d<LatticeColorMatrix>& u, 
		 const std::string& cfg_file)
{
  ArchivGauge_t header;
  Double w_plaq, link;
  mesplq(w_plaq, link, u);
  header.w_plaq = w_plaq;
  header.link = link;
  header.checksum = computeChecksum(u, header.mat_size);

  XMLReader  xml_in(xml);   // use the buffer writer to instantiate a reader
  read(xml_in, "/NERSC", header);

  writeArchiv(header, u, cfg_file);
}


// Write a Archive configuration file
/*
 * \ingroup io
 *
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeArchiv(const multi1d<LatticeColorMatrix>& u, 
		 const std::string& cfg_file)
{
  ArchivGauge_t header;
  writeArchiv(header, u, cfg_file);
}


} // namespace QDP;
