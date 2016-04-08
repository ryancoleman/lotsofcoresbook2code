/* $Id: DML_scalar.c,v 1.10 2007-06-17 13:18:59 detar Exp $ */
/* Scalar versions of QMP-dependent utilities for DML */

#include <lrl.h>
#include <dml.h>
#include <qio_stdint.h>

/* Sum a uint64_t over all nodes (for 64 bit byte counts) */
void DML_peq_uint64_t(uint64_t *subtotal, uint64_t *addend)
{
  *subtotal += *addend;
}

void DML_sum_uint64_t(uint64_t *ipt) {}

/* Sum an int over all nodes (16 or 32 bit) */
void DML_sum_int(int *ipt){}

int DML_send_bytes(char *buf, size_t size, int tonode){
  printf("ERROR: called DML_send_bytes() in DML_vanilla.c\n");
  exit(1);
  return 1;
}

int DML_get_bytes(char *buf, size_t size, int fromnode){
  printf("ERROR: called DML_get_bytes() in DML_vanilla.c\n");
  exit(1);
  return 1;
}

void DML_broadcast_bytes(char *buf, size_t size, int this_node, int from_node) {}

int DML_clear_to_send(char *scratch_buf, size_t size, 
		      int my_io_node, int new_node)
{
  printf("ERROR: called DML_clear_to_send() in DML_vanilla.c\n");
  exit(1);
  return 1;
}

int DML_route_bytes(char *buf, size_t size, int fromnode, int tonode) {
  printf("ERROR: called DML_route_bytes() in DML_vanilla.c\n");
  exit(1);
  return 1;
}

void DML_global_xor(uint32_t *x){}

void DML_sync(void){}

/* I/O layout */
int DML_io_node(const int node){
  return 0;
}

int DML_master_io_node(void){
  return 0;
}


