//////////////////////////////////////////////////////////////////
// (c) Copyright 2014-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file mkl_engine.h
 * @brief definition of fft_engine_base<fftw_real,FFTMKL_ENG>
 */
#ifndef IHPC_MKL_ENGINE_BASE_H
#define IHPC_MKL_ENGINE_BASE_H
#include <mkl_dfti.h>

namespace fft
{
  inline _MKL_Complex16* mkl_mangle(std::complex<double>* in)
  {
    return reinterpret_cast<_MKL_Complex16*>(in);
  }

  inline _MKL_Complex8* mkl_mangle(std::complex<float>* in)
  {
    return reinterpret_cast<_MKL_Complex8*>(in);
  }

  inline DFTI_CONFIG_VALUE dfti_get_precision(double a)
  {
    return DFTI_DOUBLE;
  }

  inline DFTI_CONFIG_VALUE dfti_get_precision(float a)
  {
    return DFTI_SINGLE;
  }

  /** fft_engine_base<T,N,FFTMKL_ENG> to be specialized for mkl
   */
  template<typename T, int N>
    class fft_engine_base<T,N,FFTMKL_ENG>
    {
      public:
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef DFTI_DESCRIPTOR_HANDLE fft_plan_type;
        int num_user_threads;
        ///the plan for fft & ifft for c2c
        fft_plan_type my_handle; 
        ///the plan for ifft
        fft_plan_type c2r_handle; 
        ///stride of real data
        MKL_LONG stride_r[N+1];
        ///stride of complex data
        MKL_LONG stride_c[N+1];

        ///default constructor
        fft_engine_base(int nut=1):num_user_threads(nut),my_handle(0),c2r_handle(0) { }

        ///virtual destructor to clean up the plans
        virtual ~fft_engine_base()
        {
          if(my_handle)
            DftiFreeDescriptor(&my_handle);
          if(c2r_handle)
            DftiFreeDescriptor(&c2r_handle);
        }

        /** plan for outplace, complex-to-complex  transform
        */
        void create_plan(int* ngrid, int howmany, complex_type* in, size_t in_dist)
        {
          if(create_plan_mkl(ngrid,howmany,DFTI_COMPLEX,false))
          {
            if(in_dist==size_t(0))
            {
              in_dist=ngrid[0];
              for(int i=1; i<N; ++i) in_dist*=ngrid[i];
            }
            MKL_LONG status=  DftiSetValue(my_handle,DFTI_INPUT_DISTANCE,in_dist);
            DftiCommitDescriptor(my_handle);
          }
        }


        /** plan for outplace, real-to-complex  transform
        */
        void create_plan(int* ngrid, int howmany, complex_type* in, complex_type* out
            , size_t in_dist, size_t out_dist)
        {
          if(create_plan_mkl(ngrid,howmany,DFTI_COMPLEX,true))
          {
            if(in_dist==size_t(0))
            {
              in_dist=ngrid[0];
              for(int i=1; i<N; ++i) in_dist*=ngrid[i];
            }
            MKL_LONG status=  DftiSetValue(my_handle,DFTI_INPUT_DISTANCE,in_dist);

            if(out_dist==size_t(0))
            {
              out_dist=ngrid[0];
              for(int i=1; i<N; ++i) out_dist*=ngrid[i];
            }
            status=  DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE,out_dist);
            DftiCommitDescriptor(my_handle);
          }
        }

        /** plan for real-to-complex, always outplace */
        void create_plan(int* desc, int* ngrid, real_type* in, complex_type* out, size_t in_dist=0, size_t out_dist=0)
        {
          MKL_LONG ngrid_0[N],status;
          for(int i=0; i<N; ++i) ngrid_0[i]=static_cast<MKL_LONG>(ngrid[i]);
          if(my_handle) //not too solid
          {
            MKL_LONG ngrid_1[N];
            DftiGetValue(my_handle,DFTI_LENGTHS,ngrid_1);
            if(std::equal(ngrid_1,ngrid_1+N,ngrid_0)) return;
            DftiFreeDescriptor(&my_handle);
            DftiFreeDescriptor(&c2r_handle);
          }

          if(1 == N) //this is really bad!!!
          {
            status = DftiCreateDescriptor(&my_handle,dfti_get_precision(T()),DFTI_REAL,N,ngrid_0[0]);
            status = DftiCreateDescriptor(&c2r_handle,dfti_get_precision(T()),DFTI_REAL,N,ngrid_0[0]);
          }
          else
          {
            status = DftiCreateDescriptor(&my_handle,dfti_get_precision(T()),DFTI_REAL,N,ngrid_0);
            status = DftiCreateDescriptor(&c2r_handle,dfti_get_precision(T()),DFTI_REAL,N,ngrid_0);
          }

          if(0 != status) return;

          status = DftiSetValue(my_handle,  DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
          status = DftiSetValue(c2r_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
          status = DftiSetValue(my_handle, DFTI_PLACEMENT,DFTI_NOT_INPLACE);
          status = DftiSetValue(c2r_handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
          /*
           *  rs[3] = 1;             cs[3] = 1;
           rs[2] = (N3/2+1)*2;    cs[2] = (N3/2+1);
           rs[1] = N2*(N3/2+1)*2; cs[1] = N2*(N3/2+1);
           rs[0] = 0;             cs[0] = 0;
           */
          stride_r[N]=1; stride_c[N]=1;
          int i=N-1;
          stride_r[i]=(ngrid[i]/2+1)*2; stride_c[i]=(ngrid[i]/2+1);
          while(i>0)
          {
            --i;
            stride_r[i]=stride_r[i+1]*ngrid[i]; stride_c[i]=stride_c[i+1]*ngrid[i];
          }
          stride_r[0]=0; stride_c[0]=0;

          status = DftiSetValue(my_handle,  DFTI_INPUT_STRIDES, stride_r); if (0 != status) return;
          status = DftiSetValue(my_handle,  DFTI_OUTPUT_STRIDES, stride_c); if (0 != status) return;
          status = DftiSetValue(c2r_handle, DFTI_INPUT_STRIDES, stride_c); if (0 != status) return;
          status = DftiSetValue(c2r_handle, DFTI_OUTPUT_STRIDES, stride_r); if (0 != status) return;

          status=DftiSetValue(my_handle,  DFTI_NUMBER_OF_USER_THREADS,num_user_threads);
          status=DftiSetValue(c2r_handle, DFTI_NUMBER_OF_USER_THREADS,num_user_threads);

          DftiCommitDescriptor(my_handle);
          DftiCommitDescriptor(c2r_handle);
        }

        inline void execute_fft(complex_type* inout) const
        {
          DftiComputeForward(my_handle,inout);
        }

        inline void execute_ifft(complex_type* inout) const
        {
          DftiComputeBackward(my_handle,inout);
        }

        inline void execute_fft(complex_type* in, complex_type* out) const
        {
          DftiComputeForward(my_handle,in,out);
        }

        inline void execute_ifft(complex_type* in, complex_type* out) const
        {
          DftiComputeBackward(my_handle,in,out);
        }

        inline void execute_fft(real_type* in, complex_type* out) const
        {
          DftiComputeForward(my_handle,in,out);
        }
        inline void execute_ifft(complex_type* in, real_type* out) const
        {
          DftiComputeBackward(c2r_handle,in,out);
        }

      private:

        bool create_plan_mkl(int* ngrid, int howmany, DFTI_CONFIG_VALUE c2c, bool noinplace)
        {
          MKL_LONG ngrid_0[N],status;
          for(int i=0; i<N; ++i) ngrid_0[i]=static_cast<MKL_LONG>(ngrid[i]);
          if(my_handle) //not too solid
          {
            MKL_LONG ngrid_1[N];
            DftiGetValue(my_handle,DFTI_LENGTHS,ngrid_1);
            if(std::equal(ngrid_1,ngrid_1+N,ngrid_0)) return false;
            DftiFreeDescriptor(&my_handle);
          }

          if(1 == N) //this is really bad!!!
          {
            status = DftiCreateDescriptor(&my_handle,dfti_get_precision(T()),c2c,N,ngrid_0[0]);
            //status=  DftiSetValue(my_handle,DFTI_INPUT_DISTANCE,ngrid_0[0]);
          }
          else
          {
            status = DftiCreateDescriptor(&my_handle,dfti_get_precision(T()),c2c,N,ngrid_0);
            //status=  DftiSetValue(my_handle,DFTI_ORDERING,DFTI_BACKWARD_SCRAMBLED);
          }
          if(0 != status) return false; // exception

          status=DftiSetValue(my_handle,DFTI_NUMBER_OF_TRANSFORMS,howmany);
          status=DftiSetValue(my_handle,DFTI_NUMBER_OF_USER_THREADS,num_user_threads);

          if(noinplace)
            status=DftiSetValue(my_handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
          else
            status=DftiSetValue(my_handle,DFTI_PLACEMENT,DFTI_INPLACE);

          //commited by the caller. It may need to set more values, e.g., stride, offset etc
          //DftiCommitDescriptor(my_handle);
          //
          return true;
        }

        ///** create DFFI_DESCRIPTOR for complex-to-complex transformations */
        //void create_dfti_desc(int* desc)
        //{
        //  if(my_handle)
        //    return;
        //  //if(my_handle) DftiFreeDescriptor(&my_handle);
        //  DFTI_CONFIG_VALUE my_precision=dfti_get_precision(real_type());
        //  if(desc[FFT_COMPLEX])
        //  {
        //    DftiCreateDescriptor(&my_handle,my_precision,DFTI_COMPLEX,1,desc[FFT_LENGTH]);
        //    DftiSetValue(my_handle,DFTI_INPUT_DISTANCE, desc[FFT_IN_DISTANCE]);
        //    DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE, desc[FFT_OUT_DISTANCE]);
        //  }
        //  else
        //  {
        //    DftiCreateDescriptor(&my_handle,my_precision,DFTI_REAL,1,desc[FFT_LENGTH]);
        //    DftiSetValue(my_handle,DFTI_INPUT_DISTANCE, desc[FFT_IN_DISTANCE]);
        //    DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE, 2*desc[FFT_OUT_DISTANCE]);//2 for complex
        //  }
        //  if(!desc[FFT_INPLACE])
        //    DftiSetValue(my_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        //  DftiSetValue(my_handle,DFTI_NUMBER_OF_TRANSFORMS, desc[FFT_NUMBER_OF_TRANSFORMS]);
        //  DftiCommitDescriptor(my_handle);
        //}

        ///** create DFFI_DESCRIPTOR for real-to-complex/complex-to-real transformations */
        //void create_r2c_desc(int dims, int howmany, bool outplace)
        //{
        //  if(my_handle) return;
        //  //if(my_handle) DftiFreeDescriptor(&my_handle);
        //  DFTI_CONFIG_VALUE my_precision=dfti_get_precision(real_type());
        //  DftiCreateDescriptor(&my_handle,my_precision,DFTI_REAL,1,dims);
        //  if(outplace) DftiSetValue(my_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        //  if(howmany>1) DftiSetValue(my_handle,DFTI_NUMBER_OF_TRANSFORMS, howmany);
        //  DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE, dims+2);
        //  DftiCommitDescriptor(my_handle);
        //}
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $
 ***************************************************************************/
