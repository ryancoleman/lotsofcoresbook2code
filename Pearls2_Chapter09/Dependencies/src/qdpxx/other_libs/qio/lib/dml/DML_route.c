#include <qmp.h>
#include <qio_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#if ( defined(HAVE_QMP_ROUTE) && defined(QIO_USE_QMP_ROUTE) )
//#warning "Using native QMP_route since it is available and enabled"

/* Use native version of QMP_route since it is available */
QMP_status_t DML_grid_route(void* buffer, size_t count,
			    size_t src, size_t dest)
{
  return QMP_route(buffer, count, src, dest);
}
#else

#ifndef QIO_USE_FAST_ROUTE

//#warning "Using Balints slow DML GRID ROUTE"

/* Private implementation of route method */
QMP_status_t DML_grid_route(void* buffer, size_t count,
			    size_t src, size_t dest)
{
  int* l_src_coords;      /* Coordinates of the source */
  int* l_dst_coords;      /* Coordinates of the destination */
  int* l_disp_vec;           /* Displacement vector dst_coords - src_coords */
  char myname[] = "DML_grid_route";
  
  QMP_msgmem_t sendbufmm, recvbufmm; /* Message memory handles */

  QMP_mem_t* sendbuf_qmp_mem_t;      /* These are the opaque blocks returned */
  QMP_mem_t* recvbuf_qmp_mem_t;      /* by QMP memory allocation...          */

  /* Now I have to use QMP_get_pointer for these */
  void *sendbuf;                   /* These are the actual comms buffers */
  void *recvbuf;                 

  int i,j;                   /* Loop Counters. Use for directions too */
  int me;                    /* My node */

  int n_hops;                      /* Number of hops */
  int  direction_sign;       /* Direction of hops */

  QMP_msghandle_t send_handle;     /* A handle for sending */
  QMP_msghandle_t recv_handle;     /* A handle for receiving */

  QMP_status_t err;                /* Error status */
  
  /* The number of dimensions in our "grid" */
  int ndim;

  size_t bufsize;

  size_t alignment = 8;


  /* Check to see if the logical topology is declared or not */
  /*
  log_top_declP = QMP_logical_topology_is_declared();
  
  if ( log_top_declP == QMP_FALSE ) { 
    fprintf(stderr, "%s: QMP logical topology MUST be declared\n",myname);
    fprintf(stderr, "It appears not to be\n");
    return QMP_TOPOLOGY_EXISTS;
  }
  */
  /* Topology is declared */
  /* Get its details */
  ndim = QMP_get_logical_number_of_dimensions();

  /* Get my node number -- use it to see whether I am source or dest */
  me = QMP_get_node_number();

  /* Must free these later I think */
  /* Allocate space for the coordinates */
  l_src_coords =(int *)QMP_get_logical_coordinates_from(src);
  if( l_src_coords == (int *)NULL ) { 
    fprintf(stderr, "%s: QMP_get_logical_coordinates_from failed\n",myname);
    return QMP_NOMEM_ERR;
  }

  l_dst_coords = (int *)QMP_get_logical_coordinates_from(dest);
  if( l_dst_coords == (int *)NULL ) { 
    fprintf(stderr, "%s: QMP_get_logical_coordinates_from failed\n",myname);
    return QMP_NOMEM_ERR;
  }

  /* Will definitely have to free this */
  l_disp_vec = (int *)malloc(sizeof(int)*ndim);
  if( l_disp_vec == (int *)NULL ) {
    fprintf(stderr, "%s: Unable to allocate displacement array\n",myname);
    return QMP_NOMEM_ERR;
  }

  /* Compute the taxi driver displacement */
  for(i=0; i < ndim; i++) {
    l_disp_vec[i] = l_dst_coords[i] - l_src_coords[i];
  }

  /* Don't need these anymore */
  /* Intent (if not yet word) of standard is that I must free these */  
  free(l_src_coords);
  free(l_dst_coords);

  /* Pad the buffers so that their lengths are always divisible by 8 */
  /* This is a funky QCDOC-ism -- maybe */
  bufsize = count;
  if( count % 8 != 0 ) { 
    bufsize += (8 - (count % 8));
  }

  /* Will have to free these with QMP_free_memory */
  sendbuf_qmp_mem_t = (QMP_mem_t *)QMP_allocate_aligned_memory(bufsize,alignment,QMP_MEM_COMMS);
  if( sendbuf_qmp_mem_t == (QMP_mem_t *)NULL ) { 
    fprintf(stderr, "%s: Unable to allocate sendbuf in QMP_route\n",myname);
    return QMP_NOMEM_ERR;
  }

  recvbuf_qmp_mem_t =(QMP_mem_t *)QMP_allocate_aligned_memory(bufsize,alignment,QMP_MEM_COMMS);
  if( recvbuf_qmp_mem_t == (QMP_mem_t *)NULL ) { 
    fprintf(stderr ,"%s: Unable to allocate recvbuf in QMP_route\n",myname);
    return QMP_NOMEM_ERR;
  }

  // Now I need the aligned pointers from these...
  sendbuf = QMP_get_memory_pointer(sendbuf_qmp_mem_t);
  recvbuf = QMP_get_memory_pointer(recvbuf_qmp_mem_t);

  /* To start with -- the first thing I have to do, is to copy
     the message into my sendbuf if I am the sender. Otherwise 
     I really don't care what junk is in there. */

  if( me == src ) {
    memcpy( sendbuf, buffer, count);
  }
  else {
    /* I don't care what my buffer contains if I am not the source
       but it may be nice to set it to zero so I don't send complete 
       garbage */

    memset( sendbuf, 0, count);
  }
  /*   
       Now Roll around
  */

  /* Declare the message memories */
  sendbufmm = QMP_declare_msgmem(sendbuf, bufsize);
  recvbufmm = QMP_declare_msgmem(recvbuf, bufsize);

  /* For each dimension do */
  for(i=0; i < ndim; i++) { 
    
    /* If the displacement in this direction is nonzero */
    if( l_disp_vec[i] != 0 ) {    

      /* Get the number of hops */
      n_hops = abs(l_disp_vec[i]);

      /* Get the direction */
      direction_sign = ( l_disp_vec[i] > 0 ?  1 : -1 );

      /* Declare relative sends in that direction */
      /* Do N Hops , in the direction. */
      /* I can re-use the handles for this direction */

      /* Create a receive handle in -direction sign */
      recv_handle = QMP_declare_receive_relative(recvbufmm, 
						 i,
						 -direction_sign,
						 0);

      if( recv_handle == NULL) { 
	fprintf(stderr, "%s: QMP_declare_receive_relative returned NULL\n",myname);
	return QMP_BAD_MESSAGE;
      }

      /* Create a send handle in direction sign */
      send_handle = QMP_declare_send_relative(sendbufmm, 
					      i, 
					      direction_sign, 
					      0);

      if( send_handle == NULL ) { 
	fprintf(stderr, "%s: QMP_declare_send_relative returned NULL\n",myname);
	return QMP_BAD_MESSAGE;
      }
	
      /* Do the hops */
      for(j=0; j < n_hops; j++) { 
	/* Start receiving */
	err = QMP_start(recv_handle);
	if(err != QMP_SUCCESS ) { 
	  fprintf(stderr, "%s: QMP_start() failed on receive in DML_orute\n",myname); 
	  return QMP_ERROR;
	}

	/* Start sending */
	err = QMP_start(send_handle);
	if(err != QMP_SUCCESS ) { 
	  fprintf(stderr, "%s: QMP_start() failed on send in QMP_route\n",myname);
	  return QMP_ERROR;
	}
	
	/* Wait for send to complete */
	err = QMP_wait(send_handle);
	if( err != QMP_SUCCESS ) { 
	  fprintf(stderr, "%s: QMP_wait() failed on send in QMP_route\n",myname);
	  return QMP_ERROR;
	}

	/* Wait for receive to complete */
	err = QMP_wait(recv_handle);
	if( err != QMP_SUCCESS ) { 
	  fprintf(stderr, "%s: QMP_wait() recv on send in QMP_route\n",myname);
	  return QMP_ERROR;
	}

	/* Copy the contents of my recvbuf into my sendbuf, 
	   ready for the next hop  -- In theory I could 
	   pointer swap here, but this is clean if slow */

	memcpy(sendbuf, recvbuf, count);

      }  /* Data is now in sendbuf */

      /* We have now done n_hops shifts. We need to change 
	 direction, so I free the message handles */
      QMP_free_msghandle(send_handle);
      QMP_free_msghandle(recv_handle);

    }

    /* Next direction */
  }
  
  /* We have now rolled around all the dimensions */
  /* The data is in the send buffer */
  /* Take it out and put it in "buffer" on the destination node only */
  if( me == dest ) { 
    memcpy(buffer, sendbuf, count);
  }
  else {
    memset(buffer, 0, count);
  }


  /* We can now free a whole  bunch of stuff */
  QMP_free_msgmem(sendbufmm);
  QMP_free_msgmem(recvbufmm);
  QMP_free_memory(sendbuf_qmp_mem_t);
  QMP_free_memory(recvbuf_qmp_mem_t);

  /* Alloced with malloc */
  free(l_disp_vec);


  return(QMP_SUCCESS);
}
#else
//#warning "Using faster DML_route from James"

static int
get_path_dir(int src, int dest, int size)
{
  int len, dir=1;

  len = dest - src;
  if(len<0) dir *= -1;
  if(2*abs(len)>size) {
    dir *= -1;
  }

  return dir;
}

/* Private implementation of route method */
QMP_status_t DML_grid_route(void* buffer, size_t count,
			    size_t src, size_t dest)
{
  int *src_coords;      /* Coordinates of the source */
  int *dst_coords;      /* Coordinates of the destination */
  int *my_coords;       /* my coordinates */
  const int *machine_size;    /* size of machine */
  int ndim, me, i;
  int on_path, path_leg;
  QMP_mem_t *mem;
  QMP_msgmem_t msgmem;
  char myname[] = "DML_grid_route";

  /* Check to see if the logical topology is declared or not */
  if(QMP_logical_topology_is_declared() == QMP_FALSE) {
    QMP_fprintf(stderr, "%s: QMP logical topology not declared\n", myname);
    return QMP_TOPOLOGY_EXISTS;
  }

  /* Topology is declared */
  /* Get its details */
  /* I don't think I should free machine size since it's const */
  ndim = QMP_get_logical_number_of_dimensions();
  machine_size = QMP_get_logical_dimensions();

  /* Get my node number -- use it to see whether I am on the path */
  me = QMP_get_node_number();

  /* Allocate space for the coordinates */
  /* Must free these later */
  src_coords = QMP_get_logical_coordinates_from(src);
  if( src_coords == NULL ) { 
    QMP_fprintf(stderr, "%s: QMP_get_logical_coordinates_from failed\n", myname);
    return QMP_NOMEM_ERR;
  }

  dst_coords = QMP_get_logical_coordinates_from(dest);
  if( dst_coords == NULL ) { 
    QMP_fprintf(stderr, "%s: QMP_get_logical_coordinates_from failed\n", myname);
    return QMP_NOMEM_ERR;
  }

  my_coords = QMP_get_logical_coordinates_from(me);
  if( my_coords == NULL ) { 
    QMP_fprintf(stderr, "%s: QMP_get_logical_coordinates_from failed\n", myname);
    return QMP_NOMEM_ERR;
  }


  /* now see if we are on the path */
  on_path = 1;

  i = 0;
  while((i<ndim)&&(my_coords[i]==dst_coords[i])) i++;
  path_leg = i;  /* which leg of the path we are on */

  i = path_leg + 1;
  while((i<ndim)&&(my_coords[i]==src_coords[i])) i++;
  if(i<ndim) on_path = 0;

  if(path_leg<ndim) {
    int dir;
    dir = get_path_dir(src_coords[path_leg], dst_coords[path_leg],
		       machine_size[path_leg]);

    if(src_coords[path_leg] <= dst_coords[path_leg]) {
      if(dir==1) {
	if( (my_coords[path_leg]<src_coords[path_leg]) ||
	    (my_coords[path_leg]>dst_coords[path_leg]) ) on_path = 0;
      } else {
	if( (my_coords[path_leg]>src_coords[path_leg]) &&
	    (my_coords[path_leg]<dst_coords[path_leg]) ) on_path = 0;
      }
    } else {
      if(dir==1) {
	if( (my_coords[path_leg]<dst_coords[path_leg]) ||
	    (my_coords[path_leg]>src_coords[path_leg]) ) on_path = 0;
      } else {
	if( (my_coords[path_leg]>dst_coords[path_leg]) &&
	    (my_coords[path_leg]<src_coords[path_leg]) ) on_path = 0;
      }
    }
  }

  if(on_path) {
    int recv_axis, recv_dir=0;
    int send_axis, send_dir=0;

    /* figure out send and recv nodes */
    if(me==src) {
      recv_axis = -1;
    } else {
      i = path_leg;
      if(i==ndim) i--;
      while(my_coords[i] == src_coords[i]) i--;
      recv_dir = -get_path_dir(src_coords[i], dst_coords[i], machine_size[i]);
      recv_axis = i;
    }

    if(me==dest) {
      send_axis = -1;
    } else {
      i = path_leg;
      send_dir = get_path_dir(src_coords[i], dst_coords[i], machine_size[i]);
      send_axis = i;
    }

    if((recv_axis<0)||(send_axis<0)) {
      mem = NULL;
      msgmem = QMP_declare_msgmem(buffer, count);
    } else {
      mem = QMP_allocate_memory(count);
      msgmem = QMP_declare_msgmem(QMP_get_memory_pointer(mem), count);
    }

    /* do recv if necessary */
    if(recv_axis>=0) {
      QMP_msghandle_t mh;
      QMP_status_t err;

      mh = QMP_declare_receive_relative(msgmem, recv_axis, recv_dir, 0);
      if(mh == NULL) { 
	QMP_fprintf(stderr, "%s: QMP_declare_receive_relative returned NULL\n",myname);
	return QMP_BAD_MESSAGE;
      }

      err = QMP_start(mh);
      if(err != QMP_SUCCESS) { 
	QMP_fprintf(stderr, "%s: QMP_start() failed on receive in DML_route\n",myname); 
	return QMP_ERROR;
      }

      err = QMP_wait(mh);
      if( err != QMP_SUCCESS ) { 
	QMP_fprintf(stderr, "%s: QMP_wait() recv on send in DML_route\n",myname);
	return QMP_ERROR;
      }

      QMP_free_msghandle(mh);
    }

    /* do send if necessary */
    if(send_axis>=0) {
      QMP_msghandle_t mh;
      QMP_status_t err;

      mh = QMP_declare_send_relative(msgmem, send_axis, send_dir, 0);
      if(mh == NULL) { 
	QMP_fprintf(stderr, "%s: QMP_declare_receive_relative returned NULL\n",myname);
	return QMP_BAD_MESSAGE;
      }

      err = QMP_start(mh);
      if(err != QMP_SUCCESS) { 
	QMP_fprintf(stderr, "%s: QMP_start() failed on receive in DML_route\n",myname); 
	return QMP_ERROR;
      }

      err = QMP_wait(mh);
      if( err != QMP_SUCCESS ) { 
	QMP_fprintf(stderr, "%s: QMP_wait() recv on send in DML_route\n",myname);
	return QMP_ERROR;
      }

      QMP_free_msghandle(mh);
    }

    QMP_free_msgmem(msgmem);
    if(mem) QMP_free_memory(mem);

  }

  free(src_coords);
  free(dst_coords);
  free(my_coords);

  return(QMP_SUCCESS);
}
#endif

#endif
