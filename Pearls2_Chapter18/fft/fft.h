/** @file fft.h
 * @brief A master header file to define fft interfaces
 */
#ifndef IHPC_FFT_MASTER_H
#define IHPC_FFT_MASTER_H
#include <complex>

namespace fft
{
///enumeration to determine what FFT  engine is used
  enum {FFTW_ENG=0  /*< FFTW  */
      , FFTMKL_ENG=1  /* FFT MKL */
      , FFTW_DM_ENG=2  /* FFT MKL */
      , FFTMKL_DM_ENG=3  /* FFT MKL */
      , MAX_FFT_ENG /* size of this enum */
  };

  enum { FFT3D_SLAB=0, 
    FFT3D_SLAB2,
    FFT3D_SLAB3,
    FFT3D_PENCIL
  };

#if !defined(HAVE_LIBFFTW)
const int FFTW_ESTIMATE=0;
const int FFTW_MEASURE=0;
const int FFTW_PATIENT=0;
#endif

enum { FFT_INPLACE=0
    , FFT_COMPLEX=1
    , FFT_NUMBER_OF_TRANSFORMS=2
    , FFT_LENGTH=3
    , FFT_IN_DISTANCE=4
    , FFT_OUT_DISTANCE=5
    , FFT_IN_STRIDE=6
    , FFT_OUT_STRIDE=7
    , FFT_MAX
};

template<typename T1, typename T2>
struct is_complex2complex
{
  static const int value=1;
};

template<typename T>
struct is_complex2complex<T,std::complex<T> >
{
  static const int value=0;
};

/** dummy declaration of fft_engine_base */
template<typename T, int N, unsigned ENG> 
struct fft_engine_base { };

template<typename T, unsigned ENG, unsigned DIST> 
struct fft3d_engine_base { };

}
#include <fft/mkl_engine.h>
#endif
