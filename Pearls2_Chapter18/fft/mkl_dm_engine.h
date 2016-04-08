//////////////////////////////////////////////////////////////////
// (c) Copyright 2014-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file mkl_dm_engine.h
 * @brief definition of fft_engine_base<fftw_real,FFTMKL_DM_ENG>
 */
#ifndef IHPC_MKL_DM_ENGINE_BASE_H
#define IHPC_MKL_DM_ENGINE_BASE_H
#include <mkl_cdft.h>
#include <fft/mkl_engine.h>

namespace fft
{
  /** fft_engine_base<T,N,FFTMKL_ENG> to be specialized for mkl
   */
  template<typename T, int N>
    class fft_engine_base<T,N,FFTMKL_DM_ENG>
    {
      public:
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef DFTI_DESCRIPTOR_DM_HANDLE fft_plan_type;
        MPI_Comm my_comm;
        int num_user_threads;
        long my_size;
        long my_N[2],my_start[2];
        fft_plan_type my_handle;
        complex_type* work;
        

        ///default constructor
        explicit fft_engine_base(MPI_Comm comm, int nut=1)
          :my_comm(comm),num_user_threads(nut),my_handle(0),work(0) { }

        ///virtual destructor to clean up the plans
        virtual ~fft_engine_base()
        {
          if(my_handle)
            DftiFreeDescriptorDM(&my_handle);
          if(work) delete [] work;
        }

          //if(noinplace)
          //  status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
          //else
          //  status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_INPLACE);

        /** plan for outplace, complex-to-complex  transform
        */
        void create_plan(int* ngrid, int howmany, complex_type* in)
        {
          if(create_plan_mkl(ngrid,howmany,DFTI_COMPLEX))
          {
            if(work) delete [] work;
            work=new complex_type[my_size];
            MKL_LONG status= DftiSetValueDM(my_handle,CDFT_WORKSPACE,work);
            status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_INPLACE);
            DftiCommitDescriptorDM(my_handle);
          }
        }


        /** plan for outplace, complex-to-complex  transform
         *
         * This is not safe, since out needs to be sized properly
        */
        void create_plan(int* ngrid, int howmany, complex_type* in, complex_type* out)
        {
          if(create_plan_mkl(ngrid,howmany,DFTI_COMPLEX))
          {
            //no need to use workspace
            MKL_LONG status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
            DftiCommitDescriptorDM(my_handle);
          }
        }

        /** plan for outplace, real-to-complex */
        void create_plan(int* desc, int* ngrid, real_type* in, complex_type* out=0, size_t in_dist=-1, size_t out_dist=-1)
        {
          //if(create_plan_mkl(ngrid,howmany,DFTI_REAL,out && in != out))
          //  DftiCommitDescriptor(my_handle);
        }

        inline void execute_fft(complex_type* inout) const
        {
          DftiComputeForwardDM(my_handle,inout);
        }

        inline void execute_ifft(complex_type* inout) const
        {
          DftiComputeBackwardDM(my_handle,inout);
        }

        inline void execute_fft(complex_type* in, complex_type* out) const
        {
          DftiComputeForwardDM(my_handle,in,out);
        }

        inline void execute_ifft(complex_type* in, complex_type* out) const
        {
          DftiComputeBackwardDM(my_handle,in,out);
        }

        inline void execute_fft(real_type* in, complex_type* out) const
        {
          DftiComputeForwardDM(my_handle,in,out);
        }
        inline void execute_ifft(complex_type* in, real_type* out) const
        {
          DftiComputeBackwardDM(my_handle,in,out);
        }

      private:

        bool create_plan_mkl(int* ngrid, int howmany, DFTI_CONFIG_VALUE c2c)
        {
          //cleanup first
          if(my_handle) DftiFreeDescriptorDM(&my_handle);

          MKL_LONG ngrid_0[N],status;
          for(int i=0; i<N; ++i) ngrid_0[i]=static_cast<MKL_LONG>(ngrid[i]);
          if(1 == N) //this is really bad!!!
            status = DftiCreateDescriptorDM(my_comm,&my_handle,dfti_get_precision(T()),c2c,N,ngrid_0[0]);
          else
            status = DftiCreateDescriptorDM(my_comm,&my_handle,dfti_get_precision(T()),c2c,N,ngrid_0);

          if(0 != status) return false; // exception

          status=DftiSetValueDM(my_handle,DFTI_NUMBER_OF_TRANSFORMS,howmany);
          status=DftiSetValueDM(my_handle,DFTI_NUMBER_OF_USER_THREADS,num_user_threads);
          //skip the final transpose
          status=DftiSetValueDM(my_handle,DFTI_TRANSPOSE,DFTI_ALLOW);

          //if(noinplace)
          //  status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
          //else
          //  status=DftiSetValueDM(my_handle,DFTI_PLACEMENT,DFTI_INPLACE);

          /** get the local values so that local data can be created */
          status = DftiGetValueDM(my_handle,CDFT_LOCAL_SIZE,&my_size);

          status = DftiGetValueDM(my_handle,CDFT_LOCAL_NX,&my_N[0]);
          status = DftiGetValueDM(my_handle,CDFT_LOCAL_OUT_NX,&my_N[1]);

          status = DftiGetValueDM(my_handle,CDFT_LOCAL_X_START,&my_start[0]);
          status = DftiGetValueDM(my_handle,CDFT_LOCAL_OUT_X_START,&my_start[1]);

          return true;
        }

    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $
 ***************************************************************************/
