#ifndef QPHIX_QDP_PACKER_H
#define QPHIX_QDP_PACKER_H



#include "qdp.h"
#include "qphix/geometry.h"

#include "qphix/dslash_def.h"
#include "qphix/qphix_config.h"

#if defined(QPHIX_MIC_SOURCE)
#include <immintrin.h>
#endif

#ifdef QPHIX_BUILD_CLOVER
#include "qphix/clover_dslash_def.h"
#endif

using namespace QDP;

namespace QPhiX { 


  template<typename FT, int veclen, int soalen, bool compress, typename QDPGauge>
    void qdp_pack_gauge(const QDPGauge& u,  
			typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb0, 
			typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb1, 
			Geometry<FT,veclen,soalen,compress>& s)
  {
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int nvecs = s.nVecs();
    int nyg = s.nGY();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();
    

    // Shift the lattice to get U(x-mu)
    QDPGauge u_minus(4);
    for(int mu=0; mu < 4; mu++) {
      u_minus[mu] = shift(u[mu], BACKWARD, mu);
    }

    
#pragma omp parallel for collapse(4)
    for(int t = 0; t < Nt; t++) {
      for(int z = 0; z < Nz; z++) {
	for(int y = 0; y < Ny; y++) {
	  for(int s = 0; s < nvecs; s++) {
	    for(int mu = 0; mu < 4; mu++) {
	      int outer_c = 3;
	      if ( compress ) {
		outer_c = 2;
	      }
	      for(int c = 0; c < outer_c; c++) {
		for(int c2 = 0; c2 < 3; c2++) {
		  for(int x = 0; x < soalen; x++) {

		    //#ifndef USE_PACKED_GAUGES
		    //int xx = x;
		    //int block = ((t*Nz+z)*Ny+y)*nvecs+s;

		    //#endif
		    //#else // USE_PACKED_GAUGES
		    int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
		    int xx = (y%nyg)*soalen+x;
		    // #endif // USE_PACKED_GAUGES

		    int qdpsite = x + soalen*(s + nvecs*(y + Ny*(z + Nz*t)));
		    u_cb0[block][2*mu][c][c2][0][xx] = u_minus[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb0[block][2*mu][c][c2][1][xx] = u_minus[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).imag();
		    u_cb0[block][2*mu+1][c][c2][0][xx] = u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb0[block][2*mu+1][c][c2][1][xx] = u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).imag();
		    
		    
		    u_cb1[block][2*mu][c][c2][0][xx] = u_minus[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb1[block][2*mu][c][c2][1][xx] = u_minus[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).imag();
		    u_cb1[block][2*mu+1][c][c2][0][xx] = u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb1[block][2*mu+1][c][c2][1][xx] = u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).imag();
		  }   
		}   
	      }
	    }
	  }
	}
      }
    }
  }

#ifdef QPHIX_BUILD_CLOVER
  template<typename FT, int veclen, int soalen, bool compress, typename ClovTerm>
    void qdp_pack_clover(const ClovTerm& qdp_clov_in,
			 typename ClovDslash<FT,veclen,soalen,compress>::CloverBlock* cl_out,Geometry<FT,veclen,soalen,compress>& s, int cb)
  {
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int nvecs = s.nVecs();
    int nyg = s.nGY();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();
    
#pragma omp parallel for collapse(4)
    for(int t = 0; t < Nt; t++) {
      for(int z = 0; z < Nz; z++) {
	for(int y = 0; y < Ny; y++) {
	  for(int s = 0; s < nvecs; s++) {
	    for(int x = 0; x < soalen; x++) {
	      
	      int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
	      int xx = (y%nyg)*soalen+x;
	      int qdpsite = x + soalen*(s + nvecs*(y + Ny*(z + Nz*t)))+rb[cb].start();
	      
	      for(int d=0; d < 6; d++) { 
		cl_out[block].diag1[d][xx]=qdp_clov_in[qdpsite].diag[0][d].elem();
	      }
	      for(int od=0; od < 15; od++) { 
		cl_out[block].off_diag1[od][RE][xx]=qdp_clov_in[qdpsite].offd[0][od].real();
		cl_out[block].off_diag1[od][IM][xx]=qdp_clov_in[qdpsite].offd[0][od].imag();
	      }

	      for(int d=0; d < 6; d++) { 
		cl_out[block].diag2[d][xx]=qdp_clov_in[qdpsite].diag[1][d].elem();
	      }
	      for(int od=0; od < 15; od++) { 
		cl_out[block].off_diag2[od][RE][xx]=qdp_clov_in[qdpsite].offd[1][od].real();
		cl_out[block].off_diag2[od][IM][xx]=qdp_clov_in[qdpsite].offd[1][od].imag();
	      }
	    }
	  }
	}
      }
    }
  }
  
#endif  // IFDEF BUILD CLOVER


  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
  void qdp_pack_cb_spinor(const QDPSpinor& psi_in, 
			  typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi,
			  Geometry<FT,veclen,soalen,compress>& s,
			  int cb) 
  { 
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();

#pragma omp parallel for collapse(4)
      for(int t=0; t < Nt; t++) {
	for(int z=0; z < Nz; z++) {
	  for(int y=0; y < Ny; y++) {
	    for(int s=0; s < nvecs; s++) { 
	      for(int col=0; col < 3; col++)  {
		for(int spin=0; spin < 4; spin++) { 
		  for(int x=0; x < soalen; x++) { 

		    int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		    int x_coord = s*soalen + x;
		    
		    int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;
		    psi[ind][col][spin][0][x] = psi_in.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).real();
		    psi[ind][col][spin][1][x] = psi_in.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).imag();
		  
		  }
		}
	      }
	    }
	  }
	}
      }
      
  }

  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
  void qdp_pack_spinor(const QDPSpinor& psi_in, 
		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_even, 
		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_odd,
		       Geometry<FT,veclen,soalen,compress>& s) 
  {
    qdp_pack_cb_spinor(psi_in,psi_even,s,0);
    qdp_pack_cb_spinor(psi_in,psi_odd,s,1);
  }

#if 0
  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
  void qdp_pack_spinor(const QDPSpinor& psi_in, 
		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_even, 
		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_odd,
		       Geometry<FT,veclen,soalen,compress>& s) 
  { 
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();

#pragma omp parallel for collapse(4)
      for(int t=0; t < Nt; t++) {
	for(int z=0; z < Nz; z++) {
	  for(int y=0; y < Ny; y++) {
	    for(int s=0; s < nvecs; s++) { 
	      for(int col=0; col < 3; col++)  {
		for(int spin=0; spin < 4; spin++) { 
		  for(int x=0; x < soalen; x++) { 

		    int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		    int x_coord = s*soalen + x;
		    
		    int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;
		    psi_even[ind][col][spin][0][x] = psi_in.elem(rb[0].start()+qdp_ind).elem(spin).elem(col).real();
		    psi_even[ind][col][spin][1][x] = psi_in.elem(rb[0].start()+qdp_ind).elem(spin).elem(col).imag();
		    
		    psi_odd[ind][col][spin][0][x] = psi_in.elem(rb[1].start()+qdp_ind).elem(spin).elem(col).real();
		    psi_odd[ind][col][spin][1][x] = psi_in.elem(rb[1].start()+qdp_ind).elem(spin).elem(col).imag();
		  
		  }
		}
	      }
	    }
	  }
	}
      }
      
  }
#endif
  

  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
    void qdp_unpack_cb_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_packed, 
			      QDPSpinor& chi,
			      Geometry<FT,veclen,soalen,compress>& s,
			      int cb) 
  { 
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();


#pragma omp parallel for collapse(4)    
    for(int t=0; t < Nt; t++) {
      for(int z=0; z < Nz; z++) {
	for(int y=0; y < Ny; y++) {
	  for(int s=0; s < nvecs; s++) { 
	    for(int spin=0; spin < 4; spin++) { 
	      for(int col=0; col < 3; col++)  {
		for(int x=0; x < soalen; x++) { 

		  int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		  int x_coord = s*soalen + x;
		  
		  int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;
		  
		  chi.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).real() =  chi_packed[ind][col][spin][0][x];
		  chi.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).imag() =  chi_packed[ind][col][spin][1][x];

		}
	      }
	    }
	  }
	}
      }
    }
  }

  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
    void qdp_unpack_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_even, 
			   typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_odd,
			   QDPSpinor& chi,
			   Geometry<FT,veclen,soalen,compress>& s) 
  {
    qdp_unpack_cb_spinor(chi_even,chi,s,0);
    qdp_unpack_cb_spinor(chi_odd,chi,s,1);
  }

#if 0
  template<typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
    void qdp_unpack_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_even, 
			   typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_odd,
			   QDPSpinor& chi,
			   Geometry<FT,veclen,soalen,compress>& s) 
  { 

    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();


#pragma omp parallel for collapse(4)    
    for(int t=0; t < Nt; t++) {
      for(int z=0; z < Nz; z++) {
	for(int y=0; y < Ny; y++) {
	  for(int s=0; s < nvecs; s++) { 
	    for(int spin=0; spin < 4; spin++) { 
	      for(int col=0; col < 3; col++)  {
		for(int x=0; x < soalen; x++) { 

		  int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		  int x_coord = s*soalen + x;
		  
		  int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;
		  
		  chi.elem(rb[0].start()+qdp_ind).elem(spin).elem(col).real() =  chi_even[ind][col][spin][0][x];
		  chi.elem(rb[0].start()+qdp_ind).elem(spin).elem(col).imag() =  chi_even[ind][col][spin][1][x];

		  chi.elem(rb[1].start()+qdp_ind).elem(spin).elem(col).real() =  chi_odd[ind][col][spin][0][x];
		  chi.elem(rb[1].start()+qdp_ind).elem(spin).elem(col).imag() =  chi_odd[ind][col][spin][1][x];

		}
	      }
	    }
	  }
	}
      }
    }
  }
#endif

#if defined(QPHIX_MIC_SOURCE)

  // Downconvert an array of float-vecs to an array of float 16 vecs
  void downconvert_array(const float *from, half *to, const unsigned int nvecs)
  {
#pragma omp parallel for shared(from,to)
    for(int i=0; i < nvecs; i++) { 
      __m512 in = _mm512_load_ps((void *)(from+16*i));
      _mm512_extstore_ps((void *)(to+16*i),in,_MM_DOWNCONV_PS_FLOAT16,_MM_HINT_NT);
    }
  }

  // Upconvert an araray of float16 vecs to an array of float 32
  void upconvert_array(const half* from, float *to, const unsigned int nvecs)
  {
#pragma omp parallel for shared(from, to)
    for(int i=0; i < nvecs; i++) { 
      __m512 in = _mm512_extload_ps((void *)(from+16*i),_MM_UPCONV_PS_FLOAT16, _MM_BROADCAST32_NONE, _MM_HINT_T0);
      _mm512_storenrngo_ps((void *)(to+16*i),in);
    }
  }


  // Half precision packers... 
 template<int soalen, bool compress, typename QDPGauge>
   void qdp_pack_gauge(const QDPGauge& u,  
		       typename Geometry<half,16,soalen,compress>::SU3MatrixBlock *u_cb0, 
		       typename Geometry<half,16,soalen,compress>::SU3MatrixBlock *u_cb1, 
		       Geometry<half,16,soalen,compress>& s)
 {
   
   typedef typename Geometry<float,16,soalen,compress>::SU3MatrixBlock GaugeF;
   
   int latt_size[4];
   int By,Bz,NCores,Sy,Sz,PadXY,PatXYZ,MinCt;
   latt_size[0]=s.Nx();
   latt_size[1]=s.Ny();
   latt_size[2]=s.Nz();
   latt_size[3]=s.Nt();
   // Make a geometry for allocating
   Geometry<float,16,soalen,compress> g_float(latt_size,
					      s.getBy(),
					      s.getBz(),
					      s.getNumCores(),
					      s.getSy(),
					      s.getSz(),
					      s.getPadXY(),
					      s.getPadXYZ(),
					      s.getMinCt());
					    
   GaugeF* tmp_cb0 = (GaugeF *)g_float.allocCBGauge();
   GaugeF* tmp_cb1 = (GaugeF *)g_float.allocCBGauge();
   // CHECK THESE ALLOCATIONS

   // OK Now pack the float
   qdp_pack_gauge<float,16,soalen,compress,QDPGauge>(u, tmp_cb0, tmp_cb1, g_float);

   // This is copied out of the allocation routine. 
   // Basically we take all that we have allocated
   // and divide by veclen*sizeof(float) to get the number 
   // of vectors to downconvert
   unsigned int n_f_vecs = (((s.getPxyz()*s.Nt()*soalen)/16)*sizeof(GaugeF))/(16*sizeof(float));

   downconvert_array((float *)tmp_cb0, (half *)u_cb0, n_f_vecs);
   downconvert_array((float *)tmp_cb1, (half *)u_cb1, n_f_vecs);

   g_float.free(tmp_cb0);
   g_float.free(tmp_cb1);
 }

 template<int soalen, bool compress, typename QDPSpinor>
   void qdp_pack_spinor(const QDPSpinor& psi_in, 
			typename Geometry<half,16,soalen, compress>::FourSpinorBlock* psi_even, 
			typename Geometry<half,16,soalen, compress>::FourSpinorBlock* psi_odd,
			Geometry<half,16,soalen,compress>& s) 
   {
     typedef typename Geometry<float,16,soalen,compress>::FourSpinorBlock SpinorF;

     int latt_size[4];
     int By,Bz,NCores,Sy,Sz,PadXY,PatXYZ,MinCt;
     latt_size[0]=s.Nx();
     latt_size[1]=s.Ny();
     latt_size[2]=s.Nz();
     latt_size[3]=s.Nt();
     // Make a geometry for allocating
     Geometry<float,16,soalen,compress> g_float(latt_size,
						s.getBy(),
						s.getBz(),
						s.getNumCores(),
						s.getSy(),
						s.getSz(),
						s.getPadXY(),
						s.getPadXYZ(),
						s.getMinCt());
     
     
     SpinorF* tmp_cb0_alloc = (SpinorF *)g_float.allocCBFourSpinor();
     SpinorF* tmp_cb1_alloc = (SpinorF *)g_float.allocCBFourSpinor();

     SpinorF* tmp_cb0 = tmp_cb0_alloc+1;
     SpinorF* tmp_cb1 = tmp_cb1_alloc+1;

     // CHECK THESE ALLOCATIONS
     
     // OK Now pack the float
     qdp_pack_spinor<float,16,soalen,compress,QDPSpinor>(psi_in, tmp_cb0, tmp_cb1, g_float);
 
     // This is copied out of the allocation routine. 
     // Basically we take all that we have allocated
     // and divide by veclen*sizeof(float) to get the number 
     // of vectors to downconvert
 
     unsigned int n_f_vecs = ((s.getPxyz()*s.Nt())*sizeof(SpinorF))/(16*sizeof(float));


     downconvert_array((float *)tmp_cb0, (half *)psi_even, n_f_vecs);
     downconvert_array((float *)tmp_cb1, (half *)psi_odd,  n_f_vecs);
     
     g_float.free(tmp_cb0_alloc);
     g_float.free(tmp_cb1_alloc);
     
   }
 
 template<int soalen, bool compress, typename QDPSpinor>
   void qdp_unpack_spinor(typename Geometry<half,16,soalen,compress>::FourSpinorBlock* chi_even, 
			  typename Geometry<half,16,soalen,compress>::FourSpinorBlock* chi_odd,
			  QDPSpinor& chi,
			  Geometry<half,16,soalen,compress>& s) 
   {
     
     typedef typename Geometry<float,16,soalen,compress>::FourSpinorBlock SpinorF;
     
     int latt_size[4];
     int By,Bz,NCores,Sy,Sz,PadXY,PatXYZ,MinCt;
     latt_size[0]=s.Nx();
     latt_size[1]=s.Ny();
     latt_size[2]=s.Nz();
     latt_size[3]=s.Nt();
     // Make a geometry for allocating
     Geometry<float,16,soalen,compress> g_float(latt_size,
						s.getBy(),
						s.getBz(),
						s.getNumCores(),
						s.getSy(),
						s.getSz(),
						s.getPadXY(),
						s.getPadXYZ(),
						s.getMinCt());
     
     
     // FIXME: CHECK THESE ALLOCATIONS
     SpinorF* tmp_cb0_alloc = (SpinorF *)g_float.allocCBFourSpinor();
     SpinorF* tmp_cb1_alloc = (SpinorF *)g_float.allocCBFourSpinor();

     SpinorF* tmp_cb0 = tmp_cb0_alloc+1;
     SpinorF* tmp_cb1 = tmp_cb1_alloc+1;

     // This is copied out of the allocation routine. 
     // Basically we take all that we have allocated
     // and divide by veclen*sizeof(float) to get the number 
     // of vectors to downconvert
 
     unsigned int n_f_vecs = ((s.getPxyz()*s.Nt())*sizeof(SpinorF))/(16*sizeof(float));

     upconvert_array((half *)chi_even, (float *)tmp_cb0, n_f_vecs);
     upconvert_array((half *)chi_odd, (float *)tmp_cb1, n_f_vecs);
     
     // OK Now pack the float
     qdp_unpack_spinor<float,16,soalen,compress,QDPSpinor>(tmp_cb0, tmp_cb1,chi, g_float);
     
     
     g_float.free(tmp_cb0_alloc);
     g_float.free(tmp_cb1_alloc);
     

   }

#ifdef QPHIX_BUILD_CLOVER
  template<int soalen, bool compress, typename ClovTerm>
    void qdp_pack_clover(const ClovTerm& qdp_clov_in,
			 typename Geometry<half,16,soalen,compress>::CloverBlock* cl_out,
			 Geometry<half,16,soalen,compress>& s, int cb)
  {
     typedef typename Geometry<float,16,soalen,compress>::CloverBlock ClovF;
     
     int latt_size[4];
     int By,Bz,NCores,Sy,Sz,PadXY,PatXYZ,MinCt;
     latt_size[0]=s.Nx();
     latt_size[1]=s.Ny();
     latt_size[2]=s.Nz();
     latt_size[3]=s.Nt();
     // Make a geometry for allocating
     Geometry<float,16,soalen,compress> g_float(latt_size,
						s.getBy(),
						s.getBz(),
						s.getNumCores(),
						s.getSy(),
						s.getSz(),
						s.getPadXY(),
						s.getPadXYZ(),
						s.getMinCt());
				
     
     
     
     ClovF* tmp_clov = (ClovF *)g_float.allocCBClov();

     // This is copied out of the allocation routine. 
     // Basically we take all that we have allocated
     // and divide by veclen*sizeof(float) to get the number 
     // of vectors to downconvert

     unsigned int n_f_vecs = (( (s.getPxyz()*s.Nt()*soalen) / 16  )*sizeof(ClovF))/(16*sizeof(float));
	  

     // OK Now pack the float
     qdp_pack_clover<float,16,soalen,compress,ClovTerm>(qdp_clov_in, tmp_clov, g_float, cb);
     
     downconvert_array((float *)tmp_clov,(half *)cl_out,n_f_vecs);
     g_float.free(tmp_clov);

  }
#endif // Build Clover

#endif  // if defined(QPHIX_MIC_SOURCE)

};


#endif
