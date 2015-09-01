// -*- C++ -*-

/*! @file
 * @brief NERSC Gauge Connection Archive gauge support
 *
 * Functions for reading and writing gauge fields in NERSC Gauge Connection
   Archive format.

   See http://qcd.nersc.gov/ for details.
 */

#ifndef QDP_IOGAUGE_H
#define QDP_IOGAUGE_H


namespace QDP {

//! NERSC Archive gauge field header
/*!
  See http://qcd.nersc.gov/ for details.
 
  \note This does not contain the checksum value.
 */
struct ArchivGauge_t
{

  //! Initializes a NERSC Archive header with default values
  /*!
   * \ingroup io
 
   The defaults are:
   - Two-row matrix storage
   - Floating-point precision is 32 bits
   - Periodic boundary conditions
   - Sequence number 1
   - Ensemble label is "NERSC archive"
   - Creator is "QDP++"
   - Creator hardware is "QDP++"
   -  The creation date is obtained from the system clock as the time when this
   function is called.
   - The archival date is the creation date.
   - The average plaquette and link are 0.
   - The ensemble ID is "X" followed by the creation date.
   .
   */    
  ArchivGauge_t();

  multi1d<int> nrow;            /*!< Lattice size */
  multi1d<int> boundary;        /*!< Boundary conditions */

  Real  w_plaq;                 /*!< Mean normalised plaquette */
  Real  link;                   /*!< Mean link trace */

  /* assume matrix size is 12 (matrix is compressed) 
     and change if we find out otherwise */
  size_t      mat_size;         /*!< Number of floating point numbers stored
				  per link matrix. This effectively specifies
				  whether the matrix is stored in two or
				  three-row format.
				*/

  /* Our Columbia friends have sneakily defined IEEE64BIG  */
  size_t      float_size;       /*!< Floating-point precision */

  n_uint32_t  checksum;         /*!< Checksum */
  int         sequence_number;  /*!< Sequence number */
  std::string ensemble_id;      /*!< Ensemble ID */
  std::string ensemble_label;   /*!< Ensemble label */
  std::string creator;          /*!< Creator */		
  std::string creator_hardware; /*!< Creator hardware */
  std::string creation_date;	  /*!< Creation date */
  std::string archive_date;     /*!< Archive date */     
};


//! Reads a Gauge Connection header from XML into a header container
/*!
  \pre The XMLReader should contain the header information in the following
  tags:
  \verbatim
     <mat_size>...        </mat_size>
     <float_size>...      </float_size>
     <nrow>...            </nrow>
     <boundary>...        </boundary>
     <ensemble_id>...     </ensemble_id>
     <ensemble_label>...  </ensemble_label>
     <creator>...         </creator>
     <creator_hardware>...</creator_hardware>
     <creation_date>...   </creation_date>
     <archive_date>...    </archive_date>
  \endverbatim

  \param xml The container of the XML metadata
  \param path The Xpath to the tag containing the NERSC tags
  \param header The header to which the NERSC header information is written.
*/
void read(XMLReader& xml, const std::string& path, ArchivGauge_t& header);

//! Writes a Gauge Connection header from a header container into XML
/*!
  \param xml The XML container to which the metadata is written
  \param path The Xpath to the tag containing the NERSC tags
  \param header The header from which the NERSC header information is read..
  
\post The XMLWriter will contain the header information in the following
  tags:
  \verbatim
     <mat_size>...	  </mat_size>
     <nrow>...    	  </nrow>
     <float_size>...      </float_size>
     <boundary>...	  </boundary>
     <ensemble_id>...     </ensemble_id>
     <ensemble_label>...  </ensemble_label>
     <creator>...         </creator>
     <creator_hardware>...</creator_hardware>
     <creation_date>...   </creation_date>
     <archive_date>...    </archive_date>
  \endverbatim
*/
void write(XMLWriter& xml, const std::string& path, const ArchivGauge_t& header);


//! Compute simple NERSC-like checksum of a gauge field
/*
  \ingroup io
 Compute the checksum of a gauge field

  \param u          gauge configuration ( Read )
  \return checksum
*/    
n_uint32_t computeChecksum(const multi1d<LatticeColorMatrix>& u, int mat_size);


//! Reads a NERSC Gauge Connection Archive format gauge field
/*!
  \ingroup io
 The data in the file is assumed to be big-endian.
 If the host nachine is little-endian, the data is byte-swapped.
 Based on the header information, the precision of the data can be converted.
 
  \param header     A container for the Gauge Connection header metadata ( Modify )
  \param u          The gauge configuration ( Modify )
  \param file       The file name ( Read )

 \note This can handle three-row format link matrices if the
 \c DATATYPE header token has the value \c 4D_SU3_GAUGE_3x3
 or two-row format matrices if it has the value \c 4D_SU3_GAUGE
 
 The plaquette, link and checksum values are ignored.
  
 */    
void readArchiv(ArchivGauge_t& header, multi1d<LatticeColorMatrix>& u, const std::string& file);

//! Reads a NERSC Gauge Connection Archive gauge configuration file
/*!
  \ingroup io
 The data in the file is assumed to be big-endian.
 If the host nachine is little-endian, the data is byte-swapped.
 Based on the header information, the precision of the data can be converted.
 
  \param xml        A container for the Gauge Connection header metadata as XML  ( Modify )
  \param u          The gauge configuration ( Modify )
  \param cfg_file   The file name ( Read )

  \note This can handle three-row format link matrices if the
 \c DATATYPE header token has the value \c 4D_SU3_GAUGE_3x3
 or two-row format matrices if it has the value \c 4D_SU3_GAUGE

 The plaquette, link and checksum values are ignored.

*/    
void readArchiv(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const std::string& cfg_file);


//! Reads a NERSC Gauge Connection Archive gauge configuration file
/*!
 \ingroup io
 The data in the file is assumed to be big-endian.
 If the host nachine is little-endian, the data is byte-swapped.

 \param u          The gauge configuration ( Modify )
 \param cfg_file   The file name ( Read )

  \note This can handle three-row format link matrices if the
 \c DATATYPE header token has the value \c 4D_SU3_GAUGE_3x3
 or two-row format matrices if it has the value \c 4D_SU3_GAUGE

 The plaquette, link and checksum values are ignored.
*/

void readArchiv(multi1d<LatticeColorMatrix>& u, const std::string& cfg_file);


//! Writes a NERSC Gauge Connection Archive gauge configuration file
/*!
 * \ingroup io
 The data is written in 32-bit big-endian IEEE format to the file.
 If the host nachine is little-endian, the data is byte-swapped.

  \param header     A container for the Gauge Connection header metadata
  \param u          The gauge configuration 
  \param file       The file name 

  \pre The information in the header should be filled in.
 
 \note The value 0 is written as checksum.
 */    
void writeArchiv(ArchivGauge_t& header, const multi1d<LatticeColorMatrix>& u, const std::string& file);


//! Writes a NERSC Gauge Connection Archive gauge configuration file
/*!
 * \ingroup io
 The data is written in 32-bit big-endian IEEE format to the file.
 If the host nachine is little-endian, the data is byte-swapped.

  \param header     A container for the Gauge Connection header metadata as XML
  \param u          The gauge configuration 
  \param cfg_file       The file name

 \pre The information in the header should be filled in.
 
 \note The value 0 is written as checksum.
 \note The token \c FLOATING_POINT is always given the value \c IEEE32BIG
 */    
void writeArchiv(XMLBufferWriter& xml, const multi1d<LatticeColorMatrix>& u, 
		 const std::string& cfg_file);

//! Writes a NERSC Gauge Connection Archive gauge configuration file
/*!
 * \ingroup io
 The data is written in 32-bit big-endian IEEE format to the file.
 If the host nachine is little-endian, the data is byte-swapped.

  \param u          The gauge configuration 
  \param cfg_file       The file name

  \note Since no header information is supplied, the default ArchivGauge_t
  values are used.
*/    

void writeArchiv(const multi1d<LatticeColorMatrix>& u, 
		 const std::string& cfg_file);


} // namespace QDP

#endif
