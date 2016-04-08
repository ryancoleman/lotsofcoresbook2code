
/*! @file
 * @brief Scalar specific routines
 * 
 * Routines for scalar implementation
 */

#include "qdp.h"
#include "qdp_util.h"

namespace QDP {

//-----------------------------------------------------------------------------
//! Initializer for generic map constructor
void Map::make(const MapFunc& func)
{
#if QDP_DEBUG >= 3
  QDP_info("Map::make");
#endif

  //--------------------------------------
  // Setup the communication index arrays
  goffsets.resize(Layout::vol());

  /* Get the offsets needed for neighbour comm.
     * goffsets(position)
     * the offsets contain the current site, i.e the neighbour for site i
     * is  goffsets(i,dir,mu) and NOT  i + goffset(..) 
     */
  // nrow is not used here. Comment it out to satisfy -Wall 
  // const multi1d<int>& nrow = Layout::lattSize();

  // Loop over the sites on this node
#pragma omp parallel for
  for(int linear=0; linear < Layout::vol(); ++linear)
  {
    // Get the true lattice coord of this linear site index
    multi1d<int> coord = Layout::siteCoords(0, linear);

    // Source neighbor for this destination site
    multi1d<int> fcoord = func(coord,+1);

    // Source linear site and node
    goffsets[linear] = Layout::linearSiteIndex(fcoord);
  }

#if 0
  for(int ipos=0; ipos < Layout::vol(); ++ipos)
    fprintf(stderr,"goffsets(%d,%d,%d) = %d\n",ipos,goffsets(ipos));
#endif
}



//-----------------------------------------------------------------------
// Compute simple NERSC-like checksum of a gauge field
/*
 * \ingroup io
 *
 * \param u          gauge configuration ( Read )
 *
 * \return checksum
 */    

uint32_t computeChecksum(const multi1d<LatticeColorMatrix>& u,
			   int mat_size)
{
  size_t size = sizeof(REAL32);
  size_t su3_size = size*mat_size;
  uint32_t checksum = 0;   // checksum
  const int nodeSites = Layout::sitesOnNode();

  multi1d<multi1d<ColorMatrix> > sa(Nd);   // extract gauge fields

  for(int dd=0; dd<Nd; dd++)        /* dir */
  {
    sa[dd].resize(nodeSites);
    QDP_extract(sa[dd], u[dd], all);
  }

  char  *chk_buf = new(std::nothrow) char[su3_size];
  if( chk_buf == 0x0 ) { 
    QDP_error_exit("Unable to allocate chk_buf\n");
  }

  for(int linear=0; linear < nodeSites; ++linear)
  {
    for(int dd=0; dd<Nd; dd++)        /* dir */
    {
      switch (mat_size)
      {
      case 12:
      {
	REAL32 su3[2][3][2];

	for(int kk=0; kk<Nc; kk++)      /* color */
	  for(int ii=0; ii<2; ii++)    /* color */
	  {
	    Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
	    su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	    su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	  }

	memcpy(chk_buf, &(su3[0][0][0]), su3_size);
      }
      break;

      case 18:
      {
	REAL32 su3[3][3][2];

	for(int kk=0; kk<Nc; kk++)      /* color */
	  for(int ii=0; ii<Nc; ii++)    /* color */
	  {
	    Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
	    su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	    su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	  }

	memcpy(chk_buf, &(su3[0][0][0]), su3_size);
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": unexpected size" << std::endl;
	QDP_abort(1);
      }

      // Compute checksum
      uint32_t* chk_ptr = (uint32_t*)chk_buf;
      for(int i=0; i < mat_size*size/sizeof(uint32_t); ++i)
	checksum += chk_ptr[i];
    }
  }

  delete[] chk_buf;

  // Get all nodes to contribute
  QDPInternal::globalSumArray((unsigned int*)&checksum, 1);   // g++ requires me to narrow the type to unsigned int

  return checksum;
}


//-----------------------------------------------------------------------
// Read a QCD archive file
// Read a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param cfg_in     binary writer object ( Modify )
 * \param u          gauge configuration ( Modify )
 */    

void readArchiv(BinaryReader& cfg_in, multi1d<LatticeColorMatrix>& u, 
		uint32_t& checksum, int mat_size, int float_size)
{
  ColorMatrix  sitefield;
  char *su3_buffer;

  REAL su3[Nc][Nc][2];
  checksum = 0;

  su3_buffer = new char[ Nc*Nc*2*float_size ];
  if( su3_buffer == 0x0 ) { 
    QDP_error_exit("Unable to allocate input buffer\n");
  }

  // Find the location of each site and send to primary node
  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize());
  
    for(int dd=0; dd<Nd; dd++)        /* dir */
    {
      /* Read an fe variable and write it to the BE */
      cfg_in.readArray(su3_buffer, float_size, mat_size);

      if (cfg_in.fail()) {
	QDP_error_exit("Error reading configuration");
      }


      // Compute checksum
      uint32_t* chk_ptr = (uint32_t*)su3_buffer;
      for(int i=0; i < mat_size*float_size/sizeof(uint32_t); ++i)
	checksum += chk_ptr[i];


      /* Transfer from input buffer to the actual su3 buffer, 
	 downcasting it to float if necessary */
      if ( float_size == 4 ) 
      { 
	REAL32 *su3_bufp = (REAL32 *)su3_buffer;
	REAL *su3_p = (REAL *)su3;

	for(int cp_index=0; cp_index < mat_size; cp_index++) { 
	  su3_p[cp_index] = (REAL)su3_bufp[cp_index];
	}
      }
      else if ( float_size == 8 ) 
      {
	REAL64 *su3_bufp = (REAL64 *)su3_buffer;
	REAL  *su3_p = (REAL *)su3;

	for(int cp_index =0; cp_index < mat_size; cp_index++) { 
	  
	  su3_p[cp_index] = (REAL)su3_bufp[cp_index];
	}
      }

      /* Reconstruct the third column  if necessary */
      if (mat_size == 12) 
      {
	su3[2][0][0] = su3[0][1][0]*su3[1][2][0] - su3[0][1][1]*su3[1][2][1]
	   	     - su3[0][2][0]*su3[1][1][0] + su3[0][2][1]*su3[1][1][1];
	su3[2][0][1] = su3[0][2][0]*su3[1][1][1] + su3[0][2][1]*su3[1][1][0]
		     - su3[0][1][0]*su3[1][2][1] - su3[0][1][1]*su3[1][2][0];

	su3[2][1][0] = su3[0][2][0]*su3[1][0][0] - su3[0][2][1]*su3[1][0][1]
		     - su3[0][0][0]*su3[1][2][0] + su3[0][0][1]*su3[1][2][1];
	su3[2][1][1] = su3[0][0][0]*su3[1][2][1] + su3[0][0][1]*su3[1][2][0]
		     - su3[0][2][0]*su3[1][0][1] - su3[0][2][1]*su3[1][0][0];
          
	su3[2][2][0] = su3[0][0][0]*su3[1][1][0] - su3[0][0][1]*su3[1][1][1]
		     - su3[0][1][0]*su3[1][0][0] + su3[0][1][1]*su3[1][0][1];
	su3[2][2][1] = su3[0][1][0]*su3[1][0][1] + su3[0][1][1]*su3[1][0][0]
		     - su3[0][0][0]*su3[1][1][1] - su3[0][0][1]*su3[1][1][0];
      }

      /* Copy into the big array */
      for(int kk=0; kk<Nc; kk++)      /* color */
      {
	for(int ii=0; ii<Nc; ii++)    /* color */
	{
	  Real re = su3[ii][kk][0];
	  Real im = su3[ii][kk][1];
	  Complex sitecomp = cmplx(re,im);
	  pokeColor(sitefield,sitecomp,ii,kk);
	}
      }

      pokeSite(u[dd], sitefield, coord);
    }
  }
  delete [] su3_buffer;
}



//-----------------------------------------------------------------------
// Write a QCD archive file
// Write a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param cfg_out    binary writer object ( Modify )
 * \param u          gauge configuration ( Read )
 */    
void writeArchiv(BinaryWriter& cfg_out, const multi1d<LatticeColorMatrix>& u,
		 int mat_size)
{
  ColorMatrix  sitefield;
  float su3[3][3][2];

  // Find the location of each site and send to primary node
  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize());

    for(int dd=0; dd<Nd; dd++)        /* dir */
    {
      sitefield = peekSite(u[dd], coord);

      if ( mat_size == 12 ) 
      {
	for(int kk=0; kk < Nc; kk++)      /* color */
	  for(int ii=0; ii < Nc-1; ii++)    /* color */
	  {
	    Complex sitecomp = peekColor(sitefield,ii,kk);
	    su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	    su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	  }
      }
      else
      {
	for(int kk=0; kk < Nc; kk++)      /* color */
	  for(int ii=0; ii < Nc; ii++)    /* color */
	  {
	    Complex sitecomp = peekColor(sitefield,ii,kk);
	    su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	    su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	  }
      }

      // Write a site variable
      cfg_out.writeArray((char *)&(su3[0][0][0]),sizeof(float), mat_size);
    }
  }

  if (cfg_out.fail())
    QDP_error_exit("Error writing configuration");
}


} // namespace QDP;
