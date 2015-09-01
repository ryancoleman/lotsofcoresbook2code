#include "qdp.h"
#include <string.h>

namespace QDPGlobalSums {

  // Given an array x of length elements
  // This routine will compute:
  //        y[i] = sum x[i]  for each individual i
  // ie the return value will also be an array of length element.
  //
  // Implementation note: I create a table all_data[nproc][length]
  // I do all my message passing first, then I sum.
  // 
  // Consider when no_of_procs in dim=3, length=1
  //
  // after the two communications the procs will have teh following all
  // data 
  //
  // proc0             proc1          proc2
  //  proc0_data       proc1_data     proc2_data
  //  proc2_data       proc0_data     proc1_data
  //  proc1_data       proc2_data     proc0_data
  //
  // where proci_data is the data originating from processor i.
  // (receives are always done from the minus direction, sends to plus)
  // 
  // I then sum on each processor, starting the sum from element i,
  // where i is the index of the processor 
  // so on proc0 I start the sum from all_data element 0
  //    on proc1                 from all_data element 1 (proc0_data)
  //    on proc2                 from all_data element 2 (proc0_data)
  // as the summation index becomes too big, I wrap it around.
  //
  // This way, all the processors carry out their sums in exactly the
  // same order, which should yield binary exactness

  template<typename T>
  QMP_status_t sumTDirection(T* x, int length, int dim)
  {
    int blocksize = sizeof(T)*length;
    if (blocksize % 8 != 0) { 
	blocksize += (8 - blocksize%8);
    }
    // Allocate space for message to send
    // Communicate all the data at once (ie sizeof(T)*length bytes)
    QMP_mem_t* send_mem = QMP_allocate_aligned_memory(blocksize,
						      8, 
						      (QMP_MEM_COMMS|QMP_MEM_FAST));
    if( send_mem == 0x0 ) {
      return QMP_NOMEM_ERR;
    }

    QMP_mem_t* recv_mem = QMP_allocate_aligned_memory(blocksize,
						      8, 
						      (QMP_MEM_COMMS|QMP_MEM_FAST));

    if( recv_mem == 0x0 ) { 
      return QMP_NOMEM_ERR;
    }

    // In addition I need to send to PLUS dir and receive from 
    // Minus dir
    void *sendmem_pointer = QMP_get_memory_pointer(send_mem);
    void *recvmem_pointer = QMP_get_memory_pointer(recv_mem);

    // (I leave the trailers full of junk)
    // memset(sendmem_pointer, 0x0, blocksize);
    // memset(recvmem_pointer, 0x0, blocksize);

    QMP_msgmem_t send_msgmem = QMP_declare_msgmem( sendmem_pointer, 
						    blocksize) ;

    QMP_msgmem_t recv_msgmem = QMP_declare_msgmem( recvmem_pointer, 
						     blocksize) ;
    // Send to + dir
    QMP_msghandle_t send_handle = QMP_declare_send_relative(send_msgmem, dim, +1, 0);

    // Recv from -dir
    QMP_msghandle_t recv_handle = QMP_declare_receive_relative(recv_msgmem, dim, -1, 0);

    // Get the number of CPU-s in this direction
    // Do I need to free these?
    const int* logical_dimensions = QMP_get_logical_dimensions();
    const int* logical_coordinates = QMP_get_logical_coordinates();
    
    int procs_in_dimension = logical_dimensions[dim];

    multi2d<T> all_data(length, procs_in_dimension);
    

    // Copy my data to send buffer to stat with
    memcpy(sendmem_pointer, x, sizeof(T)*length);

    // Copy my data to all_data 
    for(int j=0; j < length; j++) { 
      all_data(j,0) = x[j];
    }

    for(int i=0; i < procs_in_dimension-1; i++) { 
      QMP_status_t status;

      // Start receiving from -dir
      status  = QMP_start(recv_handle);
      if( status != QMP_SUCCESS ) { 
	return status;
      }

      // Start sending to +dir
      status = QMP_start(send_handle);
      if( status != QMP_SUCCESS ) { 
	return status;
      }

      // Finish the send
      status = QMP_wait(send_handle);
      if( status != QMP_SUCCESS ) { 
	return status;
      }

      status = QMP_wait(recv_handle);
      if( status != QMP_SUCCESS ) { 
	return status;
      }

      // Copy what I have received so I can send it on
      memcpy(sendmem_pointer, recvmem_pointer, sizeof(T)*length);

      for(int j=0; j < length; j++) { 
	all_data(j,i+1) = ((T *)recvmem_pointer)[j];
      }
    }


    for(int j=0; j < length; j++) { 
    // The index of the data received from 0 in this dimension
      int my_index = logical_coordinates[dim];
      x[j] = all_data(j, my_index);
      for(int i=0; i < procs_in_dimension-1; i++) {
	// Increment pointer with wraparound
	my_index = (my_index + 1) % procs_in_dimension;
	x[j] += all_data(j, my_index);
      }
    }

    // Free up the comms
    QMP_free_msghandle(recv_handle);
    QMP_free_msghandle(send_handle);
    QMP_free_msgmem(recv_msgmem);
    QMP_free_msgmem(send_msgmem);

    QMP_free_memory(recv_mem);
    QMP_free_memory(send_mem);
    // free(logical_dimensions);
    // free(logical_coordinates);
    return QMP_SUCCESS;
  }

  // Global sum: call sumTDirection in all available directions
  template<typename T>
  QMP_status_t sumT(T* x, int length)
  {
    // Get the number of dimensions
    int ndim = QMP_get_logical_number_of_dimensions();

    for(int dim=0; dim < ndim; dim++) { 
      QMP_status_t  status = sumTDirection<T>(x, length, dim);
      if (status != QMP_SUCCESS ) { 
	return status;
      }
    }

    return QMP_SUCCESS;
  }

  QMP_status_t QDP_sum_int(int *i) { 
#ifndef USE_QDP_QMP_GLOBAL_SUM
    return QMP_sum_int(i);
#else 
    return sumT<int>(i, 1);
#endif
  }

  QMP_status_t QDP_sum_float_array(float *x, int length) {
#ifndef USE_QDP_QMP_GLOBAL_SUM
    return QMP_sum_float_array(x,length);
#else
    return sumT<float>(x,length);
#endif
  }

  QMP_status_t QDP_sum_double_array(double *x, int length) { 
#ifndef USE_QDP_QMP_GLOBAL_SUM
    return QMP_sum_double_array(x,length);
#else
    return sumT<double>(x,length);
#endif
  }

}; // End namespace
