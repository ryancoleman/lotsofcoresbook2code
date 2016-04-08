/* DML_parscalar.c */
/* QMP-dependent utilities for DML */

#include <lrl.h>
#include <qmp.h>
#include <dml.h>
#include <qio_stdint.h>

#ifndef HAVE_BGL
/* Hardwire this */
#define HAVE_BGL
#endif

/* Sum a uint64_t over all nodes (for 64 bit byte counts) */

void DML_peq_uint64_t(uint64_t *subtotal, uint64_t *addend)
{
  *subtotal += *addend;
}

void DML_sum_uint64_t(uint64_t *ipt)
{
  uint64_t work = *ipt;
  QMP_binary_reduction((void *)(&work), sizeof(work), 
		       (QMP_binary_func)DML_peq_uint64_t);
  *ipt = work;
}

/* Sum an int over all nodes (16 or 32 bit) */
void DML_sum_int(int *ipt)
{
  int work = *ipt;
  QMP_sum_int(&work);
  *ipt = work;
}

int DML_send_bytes(char *buf, size_t size, int tonode){
  QMP_msgmem_t mm;
  QMP_msghandle_t mh;

  mm = QMP_declare_msgmem(buf, size);
  mh = QMP_declare_send_to(mm, tonode, 0);
  QMP_start(mh);
  QMP_wait(mh);
  QMP_free_msghandle(mh);
  QMP_free_msgmem(mm);
  return 0;
}

int DML_get_bytes(char *buf, size_t size, int fromnode){
  QMP_msgmem_t mm;
  QMP_msghandle_t mh;

  mm = QMP_declare_msgmem(buf, size);
  mh = QMP_declare_receive_from(mm, fromnode, 0);
  QMP_start(mh);
  QMP_wait(mh);
  QMP_free_msghandle(mh);
  QMP_free_msgmem(mm);
  return 0;
}


void DML_broadcast_bytes(char *buf, size_t size, int this_node, int from_node)
{
  /* QMP broadcasts from node 0 only!  So move the data there first */
  if(from_node != 0){
    if(this_node == from_node)
      DML_send_bytes(buf, size, 0);
    if(this_node == 0)
      DML_get_bytes(buf, size, from_node);
  }
  
  QMP_broadcast(buf, size);
}

int DML_clear_to_send(char *scratch_buf, size_t size, 
		      int my_io_node, int new_node) 
{
#ifndef HAVE_BGL
  if (QMP_get_msg_passing_type() != QMP_GRID)
#endif
  {
    int this_node = QMP_get_node_number();

    if(this_node == my_io_node && new_node != my_io_node)
      DML_send_bytes(scratch_buf,size,new_node); /* junk message */

    if(this_node == new_node && new_node != my_io_node)
      DML_get_bytes(scratch_buf,size,my_io_node);
  }

  return 0;
}


int DML_route_bytes(char *buf, size_t size, int fromnode, int tonode) 
{
 /* Check if we are BGL native somehow... */
 /* EG: We can use the HAVE_BGL option. Currently this is specified
    in the QMP MPI headres (in the config.h.in  whic is included
    by the COMMON_P_H etc etc. However, one can always just make it a QIO
    define too, it is really no big deal.

    Non MPI QIOs probably dont have this defined - as always check before use
    caveat emptor etc etc. Biggest danger, have BGL defined and it meant something
    else, this can be circumvented by making this a QIO define eg: QIO_HAVE_BLUEGENE 
  */
#ifdef HAVE_BGL
 int grid_route_p = QMP_SWITCH;
#else
 int grid_route_p = QMP_get_msg_passing_type();
#endif

 if (grid_route_p == QMP_GRID)
 {
   DML_grid_route(buf, size, fromnode, tonode);
 }
 else
 {
   int this_node = QMP_get_node_number();

   if (this_node == tonode)
     DML_get_bytes(buf,size,fromnode);

   if (this_node == fromnode)
     DML_send_bytes(buf,size,tonode);
 }

 return 0;
} 

void DML_global_xor(uint32_t *x){
  unsigned long work = (unsigned long)*x;
  QMP_xor_ulong(&work);
  *x = (uint32_t)work;
}

void
DML_sync(void)
{
#if 1
  QMP_barrier();
#else
  // JCO -- not sure why this was here
  int my_int=5;
  QMP_sum_int(&my_int);
#endif
}

/* I/O layout */
/* Default choices.  Otherwise, set by user.  See QIO_Filesystem */
int DML_io_node(const int node){
  return QMP_io_node(node);
}

int DML_master_io_node(void){
  return QMP_master_io_node();
}
