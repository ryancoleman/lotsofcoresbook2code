#ifndef QPHIX_BLAS_UTILS_H
#define QPHIX_BLAS_UTILS_H

namespace QPhiX {

  namespace BLASUtils { 
    
    //  res = alpha x + y
    //  alpha is complex
    template<typename FT, int S>
    inline void 
    cmadd(FT res[2][S], FT alpha[2], FT x[2][S], FT y[2][S]) 
    {
      //  (a[RE] x[RE] - a[IM] y[IM])  + res[RE]
      //  (a[RE] y[IM] + a[IM] y[RE])  + res[IM]
#pragma simd
      for(int s=0; s < S; s++) { 
	res[0][s] = alpha[0]*x[0][s] - alpha[1]*x[1][s]  + y[0][s];
	res[1][s] = alpha[0]*x[1][s] + alpha[1]*x[0][s]  + y[1][s];
      }
    }
    
    // res = -alpha x + y
    // res = y - alpha x
    // alpha is complex
    template<typename FT, int S>
    inline void 
    cnmadd(FT res[2][S], FT alpha[2], FT x[2][S], FT y[2][S]) 
    {
      //  res[RE] -(a[RE] x[RE] - a[IM] y[IM])
      // =res[RE] - a[RE] x[RE] + a[IM] y[IM]
      
      //  res[IM] -(a[RE] y[IM] + a[IM] y[RE])
      // =res[IM] -a[RE]y[IM] - a[IM] y[RE]
#pragma simd 
      for(int s=0; s < S; s++) { 
	res[0][s] = y[0][s] - alpha[0]*x[0][s] + alpha[1]*x[1][s];
	res[1][s] = y[1][s] - alpha[0]*x[1][s] - alpha[1]*x[0][s];
      }
    }
    
   
    // Generic stream in
    template<typename FT, int V>
    inline void
      streamInSpinor(FT* restrict dst, const FT* restrict src, int numvec) { 
      
#if defined(__MIC__)
      const int prefdist1 = 12;
      const int prefdist2 = 64;
      
      const char* prefl1base = (const char *)src+prefdist1*64;
      const char* prefl2base = (const char *)src+prefdist2*64;
#endif
#pragma vector aligned(src)
#pragma vector aligned(dst)
      for(int v=0; v < numvec; v++) { 
#if defined(__MIC__)
#endif
	
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }
    
    // Generic write  out
    template<typename FT, int V>
    inline void
      writeSpinor(FT* restrict dst, const FT* restrict src, int numvec) { 
      
#pragma vector aligned(src)
#pragma vector aligned(dst)
      for(int v=0; v < numvec; v++) { 
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }
    

    // Generic stream out
    template<typename FT, int V>
    inline void
      streamOutSpinor(FT* restrict dst, const FT* restrict src, int numvec) { 
      
#pragma vector aligned(src)
#pragma vector aligned(dst)
#pragma vector nontemporal(dst)
      for(int v=0; v < numvec; v++) { 
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }
    

    // Stream In to a different type
    template<typename FT, int V>
    inline void
      streamInSpinor(typename ArithType<FT>::Type* restrict dst, const FT* restrict src, int numvec) { 
      
#if defined(__MIC__)
      const int prefdist1 = 12;
      const int prefdist2 = 64;
      
      const char* prefl1base = (const char *)src+prefdist1*64;
      const char* prefl2base = (const char *)src+prefdist2*64;
#endif
#pragma vector aligned(src)
#pragma vector aligned(dst)
      for(int v=0; v < numvec; v++) { 
#if defined(__MIC__)
#endif
	
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }
    
    // Write out to a different type
    template<typename FT, int V>
    inline void
      writeSpinor(FT* restrict dst, const typename ArithType<FT>::Type* restrict src, int numvec) { 
      
#pragma vector aligned(src)
#pragma vector aligned(dst)
      for(int v=0; v < numvec; v++) { 
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }
    

    // Stream out to a different type
    template<typename FT, int V>
    inline void
      streamOutSpinor(FT* restrict dst, const typename ArithType<FT>::Type* restrict src, int numvec) { 
      
#pragma vector aligned(src)
#pragma vector aligned(dst)
#pragma vector nontemporal(dst)
      for(int v=0; v < numvec; v++) { 
#pragma simd
	for(int s=0; s < V; s++) {
	  dst[v*V+s]=src[v*V+s];
	}
      }
    }


#if defined(__MIC__)
#include <immintrin.h>

    // Half prec specicialize 
    template<>
    inline void
      streamInSpinor<half,16>(typename ArithType<half>::Type* restrict dst, const half* restrict src, int numvec) { 
     
      const int prefdist1 = 12;
      const int prefdist2 = 64;
      
      const char* prefl1base = (const char *)src+prefdist1*64;
      const char* prefl2base = (const char *)src+prefdist2*64;

      for(int v=0; v < numvec; v++) { 
	
	__m512 r = _mm512_extload_ps((void *)&src[v*16], _MM_UPCONV_PS_FLOAT16, _MM_BROADCAST32_NONE, _MM_HINT_T0);
	_mm512_store_ps((void *)&dst[v*16], r);

      }
    }
    
   
    template<>
    inline void
      writeSpinor<half,16>(half* restrict dst, const typename ArithType<half>::Type* restrict src, int numvec) { 
      const int prefdist1 = 12;
      const int prefdist2 = 64;
      
      const char* prefl1base = (const char *)src+prefdist1*64;
      const char* prefl2base = (const char *)src+prefdist2*64;

      const char* prefl1baseo = (const char *)dst+prefdist1*64;
      const char* prefl2baseo = (const char *)dst+prefdist2*64;
      
      for(int v=0; v < numvec; v++) { 

	__m512 r = _mm512_load_ps((void *)&src[v*16]);
	_mm512_extstore_ps((void *)&dst[v*16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_T0);
      }
    }


    template<>
    inline void
      streamOutSpinor<half,16>(half* restrict dst, const typename ArithType<half>::Type* restrict src, int numvec) { 
      const int prefdist1 = 12;
      const int prefdist2 = 64;
      
      const char* prefl1base = (const char *)src+prefdist1*64;
      const char* prefl2base = (const char *)src+prefdist2*64;
      
      for(int v=0; v < numvec; v++) { 

	__m512 r = _mm512_load_ps((void *)&src[v*16]);
	_mm512_extstore_ps((void *)&dst[v*16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NT);
      }
    }

#endif // defined MIC    

    
  }; // Namespace BLAS UTILS
  
};



#endif
