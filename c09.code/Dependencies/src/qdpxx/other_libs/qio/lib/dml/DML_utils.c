/* DML_utils.c */
/* Utilities for DML */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <assert.h>
#include <sys/types.h>
#include <qio_stdint.h>
#include <sys/time.h>
//#include <qmp.h>

#undef DML_DEBUG

// duplicate hack to avoid XLC bug
double QIO_time2 (void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  //if(node==0) printf("QIO_time %g %g\n", (double)tv.tv_sec, (double)tv.tv_usec);
  return (double)tv.tv_sec + 1e-6*(double)tv.tv_usec;
}
#define timestart(t) { double dt = QIO_time2(); /*if(this_node==0) printf("start " #t " %g\n", dt);*/ t-=dt; }
#define timestop(t) { double dt = QIO_time2(); /*if(this_node==0) printf("stop " #t " %g\n", dt);*/ t+=dt; }

/* Iterators for lexicographic order */

/* Initialize to lower bound */
void DML_lex_init(int *dim, int coords[], int latdim, int lower[])
{
  int d;
  for(d = 0; d < latdim; d++)coords[d] = lower[d];
  *dim = 0;
}

/* Recursively update the coordinate counter */
/* Return 0 when finished */
int DML_lex_next(int *dim, int coords[], int latdim, int lower[], int upper[])
{
  if(++coords[*dim] < upper[*dim]){
    *dim = 0;
    return 1;
  }
  else{
    coords[*dim] = lower[*dim];
    if(++(*dim) < latdim)return DML_lex_next(dim, coords, latdim, 
					     lower, upper);
    else return 0;
  }
}

/*------------------------------------------------------------------*/
/* Convert linear lexicographic rank to lexicographic coordinate */

void DML_lex_coords(int coords[], const int latdim, const int latsize[], 
		    const DML_SiteRank rcv_coords)
{
  int dim;
  DML_SiteRank rank = rcv_coords;

  for(dim = 0; dim < latdim; dim++){
    coords[dim] = rank % latsize[dim];
    rank /= latsize[dim];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   DML_lex_coords) */

DML_SiteRank DML_lex_rank(const int coords[], int latdim, int latsize[])
{
  int dim;
  DML_SiteRank rank = coords[latdim-1];

  for(dim = latdim-2; dim >= 0; dim--){
    rank = rank * latsize[dim] + coords[dim];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Make temporary space for coords */

int *DML_allocate_coords(int latdim, const char *myname, int this_node){
  int *coords;

  coords = (int *)malloc(latdim*sizeof(int));
  if(!coords)printf("%s(%d) can't malloc coords\n",myname,this_node);
  return coords;
}


/*------------------------------------------------------------------*/
/* The message structure holds the site datum and site rank */

/* Accessor: Size of message */
size_t DML_msg_sizeof(size_t size){
  return size + sizeof(DML_SiteRank);
}

/*------------------------------------------------------------------*/
/* Constructor*/
char *DML_allocate_msg(size_t size, char *myname, int this_node){
  char *msg;
  size_t sizeof_msg = DML_msg_sizeof(size);

  msg = (char *)malloc(sizeof_msg);
  if(!msg)printf("%s(%d) can't malloc msg\n",myname,this_node);
  return msg;
}

/*------------------------------------------------------------------*/
/* Accessor: Pointer to datum member of msg */
char *DML_msg_datum(char *msg, size_t size){
  return msg;
}

/*------------------------------------------------------------------*/
/* Accessor: Pointer to rank member of msg */
DML_SiteRank *DML_msg_rank(char *msg, size_t size){
  return (DML_SiteRank *)(msg + size);
}

/*------------------------------------------------------------------*/
/* Count the sitelist for partitioned I/O format */
/* Return code 0 = success; 1 = failure */
int DML_count_partition_sitelist(DML_Layout *layout, DML_SiteList *sites){
  int *coords;
  int latdim = layout->latdim;
  int node;
  int this_node = layout->this_node;
  int number_of_nodes = layout->number_of_nodes;
  int my_io_node = layout->ionode(this_node);
  size_t number_of_io_sites;
  int number_of_my_ionodes;
  char myname[] = "DML_count_partition_sitelist";
#if 0
  DML_SiteRank rank;
  int *latsize = layout->latsize;
  size_t volume = layout->volume;
  int send_node;
#endif

  /* Space for a coordinate vector */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords)return 1;

  /* Iterate over all nodes, adding up the sites my partition
     writes and counting the number of nodes in my partition */
  number_of_my_ionodes = 0;
  number_of_io_sites = 0;
  for(node = 0; node < number_of_nodes; node++){
    if(layout->ionode(node) == my_io_node){
      /* ( If we are discovering the lattice dimensions, we won't know
	 the number of sites on any node ) */
      if(layout->latdim != 0)
	number_of_io_sites += layout->num_sites(node);
      number_of_my_ionodes++;
    }
  }

#if 0  
  /* Iterate over all sites, storing the lexicographic index
     only for the sites on my partition */
  number_of_io_sites = 0;
  /* Loop over all sites in lexicographic order */
  for(rank = 0; rank < volume; rank++)
    {
      /* Convert rank index to coordinates */
      DML_lex_coords(coords, latdim, latsize, rank);
      /* The node containing these coordinates */
      send_node = layout->node_number(coords);
      if(layout->ionode(send_node) == my_io_node)
	number_of_io_sites++;
    }
#endif


  sites->number_of_io_sites = number_of_io_sites;
  sites->number_of_my_ionodes = number_of_my_ionodes;
  free(coords);
  return 0;
}

/*------------------------------------------------------------------*/
/* Create sitelist structure (but not the sitelist, yet) */
DML_SiteList *DML_init_sitelist(int volfmt, int serpar, DML_Layout *layout){
  DML_SiteList *sites;
  int this_node = layout->this_node;
  char myname[] = "DML_init_sitelist";

  sites = (DML_SiteList *)malloc(sizeof(DML_SiteList));
  if(!sites){
    printf("%s: Can't allocate small sitelist structure\n",myname);
    return NULL;
  }
  sites->list               = NULL;
  sites->use_list           = 0;
  sites->first              = 0;
  sites->current_rank       = 0;
  sites->number_of_io_sites = 0;
  sites->number_of_my_ionodes = 0;
  sites->current_index      = 0;
  sites->use_subset         = 0;
  sites->subset_rank        = NULL;
  sites->subset_io_sites    = 0;

  /* Initialize number of I/O sites */

  if(volfmt == DML_SINGLEFILE && serpar == DML_SERIAL){
    /* Single files are always written in lexicographic order so
       single file format doesn't need a sitelist and all nodes count
       through the entire file */
      sites->number_of_io_sites = layout->volume;
  }

  else if(volfmt == DML_MULTIFILE){
    /* Multifile format requires a sitelist for each node */

    sites->use_list = 1;

    /* Each node reads/writes its own sites independently */
    sites->number_of_io_sites = layout->sites_on_node;
  }

  else if(volfmt == DML_PARTFILE || 
	  (volfmt == DML_SINGLEFILE && serpar == DML_PARALLEL)){
    /* Partitioned I/O requires a separate sitelist for each I/O
       partition.  Parallel I/O uses a sitelist to determine which
       sites to read or write. */
    /* Here we just count the number of sites in the list. */

    sites->use_list = 1;
    if(DML_count_partition_sitelist(layout,sites)){
      free(sites); return NULL;
    }
  }
  else{
    /* Bad volfmt */
    printf("%s(%d): bad volume format code = %d\n",
	   myname, this_node, volfmt);
    free(sites); return NULL;
  }

  return sites;

}

/*------------------------------------------------------------------*/
/* Compare site lists. Return 0 if all sites in the list agree and 1
   otherwise */

int DML_compare_sitelists(DML_SiteRank *lista, DML_SiteRank *listb, size_t n){
  size_t i;
  
  /* Scan site list for my node */
  for(i = 0; i < n; i++){
    if(lista[i] != listb[i])return 1;
  }
  
  return 0;
}

/*------------------------------------------------------------------*/
/* Free the sitelist structure */
void DML_free_sitelist(DML_SiteList *sites){
  if(sites == NULL)return;
  if(sites->list != NULL)free(sites->list);
  free(sites);
}

/*------------------------------------------------------------------*/
/* Fill the sitelist for multifile format */    
/* Return code 0 = success; 1 = failure */
int DML_fill_multifile_sitelist(DML_Layout *layout, DML_SiteList *sites){
  int *coords;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  int this_node = layout->this_node;
  size_t index;
  char myname[] = "DML_fill_multifile_sitelist";

  /* Each node dumps its own sites */

  /* Space for a coordinate vector */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords)return 1;
  /* Iterate over sites in storage order on this node */
  for(index = 0; index < layout->sites_on_node; index++){
    /* Convert storage order to coordinates */
    layout->get_coords(coords,this_node,index);
    /* Convert coordinate to lexicographic rank */
    sites->list[index] = DML_lex_rank(coords,latdim,latsize);
  }

  free(coords);
  return 0;
}

/*------------------------------------------------------------------*/
/* This version was slower for small partitions of large files */

/* Fill the sitelist for partitioned I/O format */
/* Return code 0 = success; 1 = failure */
int DML_fill_partition_sitelist_try(DML_Layout *layout, DML_SiteList *sites){
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  int this_node = layout->this_node;
  int my_io_node = layout->ionode(this_node);
  size_t number_of_io_sites = sites->number_of_io_sites;
  DML_SiteRank *list        = sites->list;
  size_t index;
  int *coords;
  DML_SiteRank rank;
  DML_SiteRank volume = (DML_SiteRank)layout->volume;
  char myname[] = "DML_fill_partition_sitelist_try";

  /* Space for a coordinate vector */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords)return 1;

  /* Scan all the sites in the lattice to find the sites in our I/O
     partition */
  /* The resulting list is automatically in order */
  for( index = 0, rank = 0; rank < volume; rank++ ){
    /* Map rank to coordinates */
    DML_lex_coords(coords, latdim, latsize, rank);
    /* If we have this site, add it to the list */
    /* (find the node that has the coords to node and then its I/O node) */
    if(layout->ionode(layout->node_number(coords)) == my_io_node){
      list[index] = rank;
      index++;
      if(index >= number_of_io_sites)break;  /* memory protection */
    }
  }

  if(index != number_of_io_sites){
    printf("%s(%d) Internal error. Can't count I/O sites\n",
	   myname,this_node);
    return 1;
  }

  free(coords);
  return 0;
}

/*------------------------------------------------------------------*/
/* Heap sort after Numerical Recipes */

/* Pull out the left element and work up the left side of a binary
   tree, stopping at an element that is smaller than or equal to the
   left one.  Push down any larger elements encountered on the way.
   Finally, put the left element in the remaining empty slot -- in
   effect producing a cyclic rotation. */

void DML_sift_down(DML_SiteRank *array, const long left, const long right)
{
  long j,jold;
  DML_SiteRank a;
  
  a = array[left];
  jold = left;
  j = 2*left+1;
  while (j <= right) {
    /* Use the larger of array[j] and array[j+1] */
    if (j < right && array[j] < array[j+1]) j++;
    if (a >= array[j]) break;
    array[jold] = array[j];
    jold = j;
    j = 2*j+1;
  }
  array[jold] = a;
}

void DML_hpsort(DML_SiteRank *array, long n)
{
  long i;
  DML_SiteRank tmp;

  for( i=n/2-1; i>=0; i--)
    DML_sift_down(array,i,n-1);
  for( i = n-1; i>0; i--) {
    tmp = array[0];
    array[0] = array[i];
    array[i] = tmp;
    DML_sift_down(array,0,i-1);
  }
}

/*------------------------------------------------------------------*/
/* Fill the sitelist for partitioned I/O format */
/* Return code 0 = success; 1 = failure */
int DML_fill_partition_sitelist(DML_Layout *layout, DML_SiteList *sites){
  size_t index, node_index;
  int *coords;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  int node;
  int node_sites;
  int this_node = layout->this_node;
  int my_io_node = layout->ionode(this_node);
  int number_of_nodes = layout->number_of_nodes;
  size_t number_of_io_sites = sites->number_of_io_sites;
  DML_SiteRank *list        = sites->list;
  char myname[] = "DML_fill_partition_sitelist";
#if 0
  DML_SiteRank rank;
  size_t volume = layout->volume;
  int send_node;
#endif

  /* Space for a coordinate vector */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords)return 1;

  /* Fill the list in storage order first */
  index = 0;
  for(node = 0; node < number_of_nodes; node++){
    /* Find the nodes on my partition */
    if(layout->ionode(node) == my_io_node){
      node_sites = layout->num_sites(node);
      for(node_index = 0; node_index < node_sites; node_index++){
	layout->get_coords(coords,node,node_index);
	list[index] = DML_lex_rank(coords,latdim,latsize);
	index++;
      }
    }
  }
  if(index != number_of_io_sites){
    printf("%s(%d) Internal error. Can't count I/O sites\n",
	   myname,this_node);
    return 1;
  }

  /* Put the site list in ascending lexicographic rank order */
  DML_hpsort(list, number_of_io_sites);
  free(coords);
  return 0;
}

/*------------------------------------------------------------------*/
/* Create and populate the sitelist for output */
/* Return code 0 = success; 1 = failure */
int DML_fill_sitelist(DML_SiteList *sites, int volfmt, int serpar,
		      DML_Layout *layout){
  int this_node = layout->this_node;
  char myname[] = "DML_fill_sitelist";

  if(sites->use_list == 0)return 0;

  /* Allocate the list */

  sites->list = 
    (DML_SiteRank *)malloc(sizeof(DML_SiteRank)*sites->number_of_io_sites);
  if(sites->list == NULL)return 1;

  /* Fill the list */

  if(volfmt == DML_MULTIFILE){
  /* Multifile format requires a sitelist */
    return DML_fill_multifile_sitelist(layout,sites);
  }
  else if(volfmt == DML_PARTFILE || 
	  (volfmt == DML_SINGLEFILE && serpar == DML_PARALLEL)){
    /* Partitioned I/O requires a sitelist on the I/O node */
    /* Singlefile parallel I/O requires the same sitelist as partfile
       to determine which sites to write/read */
    return DML_fill_partition_sitelist(layout,sites);
  }
  else {
    /* Bad volfmt */
    printf("%s(%d): bad volume format code = %d\n",
	   myname, this_node, volfmt);
    return 1;
  }
  return 0;
}

/*------------------------------------------------------------------*/
/* Read and check the sitelist for input */
/* return 0 for success and 1 for failure */
int DML_read_sitelist(DML_SiteList *sites, LRL_FileReader *lrl_file_in,
		      int volfmt, DML_Layout *layout,
		      LIME_type *lime_type){
  uint64_t check, announced_rec_size;
  int this_node = layout->this_node;
  LRL_RecordReader *lrl_record_in;
  DML_SiteRank *inputlist;
  int not_ok;
  int status;
  char myname[] = "DML_read_sitelist";

  if(sites->use_list == 0)return 0;

  /* Open sitelist record */
  lrl_record_in = LRL_open_read_record(lrl_file_in, &announced_rec_size, 
				       lime_type, &status);
  if(!lrl_record_in)return 1;

  /* Require that the record size matches expectations */
  check = sites->number_of_io_sites * sizeof(DML_SiteRank);
  /* Ignore a mismatch if we are trying to discover the lattice dimension */
  if(!layout->discover_dims_mode && announced_rec_size != check){
    printf("%s(%d): sitelist size mismatch: found %lu expected %lu lime type %s\n",
	   myname, this_node, (unsigned long)announced_rec_size,
	   (unsigned long)check, *lime_type);
    printf("%s(%d): latdim = %d\n",myname, this_node,layout->latdim);
    return 1;
  }

  /* Allocate check list according to record size */
  
  inputlist = (DML_SiteRank *)malloc(announced_rec_size);
  if(inputlist == NULL)return 1;
  
  /* Read the site list and close the record */
  check = LRL_read_bytes(lrl_record_in, (char *)inputlist, announced_rec_size);
  
  LRL_close_read_record(lrl_record_in);
    
#ifdef DML_DEBUG
  printf("%s(%d) sitelist record was read with %lu bytes\n",myname,
	 layout->this_node,(unsigned long)check);
#endif
 
  /* Check bytes read */
  if(check != announced_rec_size){
    printf("%s(%d): bytes read %lu != expected rec_size %lu\n",
	   myname, this_node, (unsigned long)check, 
	   (unsigned long)announced_rec_size);
    free(inputlist); return 1;
  }
 
  /* Byte reordering for entire sitelist */
  if (! DML_big_endian())
    DML_byterevn((char *)inputlist, announced_rec_size, sizeof(DML_SiteRank));
  
  /* All input sitelists must agree exactly with what we expect */
  /* Unless we are reading in discovery mode */
  if(!layout->discover_dims_mode){
    not_ok = DML_compare_sitelists(sites->list, inputlist, 
				   sites->number_of_io_sites);
    
    if(not_ok)
      printf("%s(%d): sitelist does not conform to I/O layout.\n",
	     myname,this_node);

    /* Return 1 if not OK and 0 if OK */
    free(inputlist); 
    return not_ok;
  }
  else
    return 0;
}


/*------------------------------------------------------------------*/
/* First site for I/O processing */

DML_SiteRank DML_init_site_loop(DML_SiteList *sites){
  sites->current_index = 0;
  if(sites->use_list){
    return sites->list[sites->current_index];
  }
  else{
    sites->current_rank = sites->first;
    return sites->current_rank;
  }
}

/*------------------------------------------------------------------*/
/* Iterator for sites processed by I/O */
/* Returns 0 when iteration is complete. 1 when not and updates rank. */
int DML_next_site(DML_SiteRank *rank, DML_SiteList *sites)
{
  sites->current_index++;
  if(sites->current_index >= sites->number_of_io_sites)return 0;
  if(sites->use_list){
    *rank = sites->list[sites->current_index];
  }
  else{
    sites->current_rank++;
    *rank = sites->current_rank;
  }
  return 1;
}
 
/*------------------------------------------------------------------*/
/* Copy subset data into DML layout structure                       */
/* return 1 for failure (bad hypercube bounds) and 0 for success */
int DML_insert_subset_data(DML_Layout *layout, int recordtype,
			   int *lower, int *upper, int n)
{
  int i;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  int *hyperupper = layout->hyperupper;
  int *hyperlower = layout->hyperlower;
  int this_node = layout->this_node;
  size_t subsetvolume;


  /* Insert record type */
  layout->recordtype = recordtype;

  if(recordtype != DML_HYPER){

    /* Not a hypercube record, so the subset includes everything */
    layout->subsetvolume = layout->volume;

  } else {

    for(i = 0; i < latdim; i++){
      if(i < n){
	hyperlower[i] = lower[i];
	hyperupper[i] = upper[i];
      }
      else{
	/* This shouldn't happen, but do something sensible if it does */
	hyperlower[i] = 0;
	hyperupper[i] = 0;
      }
    }
    
    /* Compute the subset volume */
    subsetvolume = 1;
    for(i = 0; i < latdim; i++){
      if(hyperlower[i] >= 0 && 
	 hyperupper[i] <= latsize[i] && 
	 hyperlower[i] <= hyperupper[i])
	subsetvolume *= (upper[i] - lower[i] + 1);
      else {
	printf("DML_insert_subset_data(%d): Bad hypercube bounds\n", this_node);
	return 1;
      }
    }
    layout->subsetvolume = subsetvolume;
  }    
  return 0;
}

/*------------------------------------------------------------------*/
/* Check whether site specified by coordinate is within the
   hypercube subset  */
/* Returns 0 if not and 1 if so */
int DML_coord_inside_subset(int coords[], DML_Layout *layout){
  int *upper   = layout->hyperupper;
  int *lower   = layout->hyperlower;
  int latdim   = layout->latdim;
  int i, status;

  /* For QIO_FIELD records, we take all sites (no subset) */
  if(layout->recordtype == DML_FIELD)return 1;

  status = 1;
  for(i = 0; i < latdim; i++)
    if(lower[i] > coords[i] || upper[i] < coords[i]){
      status = 0;
      break;
    }

  return status;
}


/*------------------------------------------------------------------*/
/* Check whether site specified by lexicographic rank is within the
   hypercube subset  */
/* Returns 0 if not and 1 if so */
int DML_rank_inside_subset(DML_SiteRank rank, DML_Layout *layout){
  int latdim   = layout->latdim;
  int *latsize = layout->latsize;
  int *coords;
  char myname[] = "DML_rank_inside_subset";
  int this_node = layout->this_node;
  int status;

  /* For QIO_FIELD records, we take all sites (no subset) */
  if(layout->recordtype == DML_FIELD)return 1;

  /* Allocate lattice coordinate */
  coords = DML_allocate_coords(latdim, myname, this_node);
  /* Convert lexicographic rank to coordinates */
  DML_lex_coords(coords, latdim, latsize, rank);

  status = DML_coord_inside_subset(coords, layout);

  free(coords); 
  return status;
}


/*------------------------------------------------------------------*/
/* Table lookup for sorted table.  Return index if found and -1 if not
   found. Binary search for exact match. */

int DML_table_lookup(DML_SiteRank list[], size_t n, DML_SiteRank r)
{
  int ju,jm,jl;
  
  if ( n == 0) return -1;
  if ( r < list[0] || r > list[n-1] ) return -1;
  if ( r == list[0] ) return 0;
  if ( n == 1 ) return -1;
  if ( r == list[n-1] ) return n-1;
  if ( n == 2 ) return -1;

  jl=0;
  ju=n-1;
  while (ju-jl > 1) {
    jm=(ju+jl)/2;
    if ( r == list[jm] ) return jm;
    if ( r > list[jm]  ) jl=jm;
    else                 ju=jm;
  }
  return -1;
}

/*------------------------------------------------------------------*/
/* Find the physical location in the record of the site with
   lexicographic index "rank".  Returns 1 on error.  0 for success. */

int DML_lookup_subset_rank(DML_SiteRank *seek, DML_SiteRank rank, 
			   DML_SiteList *sites)
{
  int current_index;
  int status;

  if(sites->use_list){
    current_index = DML_table_lookup(sites->list, sites->number_of_io_sites, 
				     rank);
    if(current_index < 0)
      return 1;
  }
  else{
    current_index = rank;
  }

  if(sites->use_subset){
    status = sites->subset_rank[current_index];
    if(status < 0)return 1;
    *seek = status;
  } else {
    *seek = current_index;
  }
  return 0;
}

/*------------------------------------------------------------------*/
/* Iterator for I/O sites in subset */
/* Returns 0 when iteration is complete. 1 when not and updates rank. */
int DML_next_subset_site(DML_SiteRank *rank, DML_SiteList *sites)
{
  int status;

  status = DML_next_site(rank, sites);
  if(sites->use_subset)
    while(sites->subset_rank[sites->current_index] < 0){
      if( (status = DML_next_site(rank, sites)) == 0) break;
    }

  return status;
}
 
/*------------------------------------------------------------------*/
/* Return position of site in record - same as subset rank          */
/* Note: this routine is used only in conjunction with the
   DML_next_subset_site iterator.  */
int
DML_subset_rank(DML_SiteRank rank, DML_SiteList *sites)
{
  //printf("node %i rank %i\n", QMP_get_node_number(), rank);
  if(sites->use_subset)
    return sites->subset_rank[sites->current_index];
  else
    return rank;
}

/*------------------------------------------------------------------*/
/* Initialize site iterator for I/O processing */
/* Returns 0 when there are no sites to process. Otherwise 1. */
int DML_init_subset_site_loop(DML_SiteRank *rank, DML_SiteList *sites){

  /* Our first site */
  *rank = DML_init_site_loop(sites);

  /* If we are doing a subset and our first site is not in that
     subset, scan forward to find the first site in the subset */
  if(sites->use_subset)
    if(sites->subset_rank[sites->current_index] < 0)
      return DML_next_subset_site(rank, sites);
  return 1;
}

/*------------------------------------------------------------------*/
/* Search an ordered table of site ranks for a given rank.  Return the
   index of the match, if found, or -1 if not found */

int DML_lookup_site_rank(DML_SiteRank *t, int n, DML_SiteRank r){
  int s1 = 0;
  int s2 = n-1;
  int snew;

  if(t[s1] > r || t[s2] < r)return -1;
  if(t[s1] == r)return s1;
  if(t[s2] == r)return s2;
  /* Binary search */
  while(s2 > s1 + 1){
    snew = (s1 + s2)/2;
    if(t[snew] == r)return snew;
    if(t[snew] > r)s2 = snew;
    else s1 = snew;
  }
  return -1;
}
/*------------------------------------------------------------------*/
/* Create the subset rank list for single-file parallel I/O.        */
/* See DML_create_subset_rank below for a definition of this list   */

/* Because we are dealing with a single file, the data that this I/O
   partition handles (if any) may be scattered throughout the record.
   They are arranged in lexicographic coordinate order, but with
   omissions, because this is a subset record.  So our data locations
   do not follow a regular pattern of offsets.  Therefore, we have to
   count through all the sites in the subset and not just those
   belonging to our I/O partition to determine the data locations.
   Even though we are dealing with a subset, we create a subset_rank
   list with an entry for all the sites in our partition.  We simply
   mark the entries that are not in the subset. */

/* Return value 0 for success and 1 for malloc failure */
int DML_create_subset_rank_parallel(DML_SiteList *sites, DML_Layout *layout){

  int latdim   = layout->latdim;
  int *latsize = layout->latsize;
  int this_node = layout->this_node;
  int *upper = layout->hyperupper;
  int *lower = layout->hyperlower;
  DML_SiteRank r, *ranklist;
  int *coords,*ubound;
  int d,s,t, dim;
  char myname[] = "DML_create_subset_rank_parallel";

  sites->use_subset = 1;

  sites->subset_rank = 
    (int *)malloc(sizeof(int)*sites->number_of_io_sites); /* Could be less */
  if(sites->subset_rank == NULL)return 1;

  ranklist = 
    (DML_SiteRank *)malloc(sizeof(DML_SiteRank)*sites->number_of_io_sites);
  if(ranklist == NULL)return 1;

  /* Allocate lattice coordinate */
  coords = DML_allocate_coords(latdim, myname, this_node);

  /* Tabulate all sites in our I/O partition in lexicographic
     (processing) order.  */
  r = DML_init_site_loop(sites);
  s = 0;
  do {
    /* At first flag them as outside the subset. */
    sites->subset_rank[s] = -1;
    /* The rank of the site in the whole lattice */
    ranklist[s] = r;
    s++;
  } while(DML_next_site(&r, sites));

  sites->subset_io_sites = s;

  /* Upper limits for DML_lex_next iteration */
  ubound  = DML_allocate_coords(latdim, myname, this_node);
  for(d = 0; d < latdim; d++)
    ubound[d] = upper[d]+1;

  /* Iterate over all sites in the subset in lexicographic order, i.e. the
     sequence in the file record */

  t = 0;   /* Rank of site in the record */
  DML_lex_init(&dim, coords, latdim, lower);
  do {
    /* The rank of this site in the whole lattice */
    r = DML_lex_rank(coords, latdim, latsize);
    /* Is this site in our I/O partition? If so, get its rank in our
       sitelist */
    s = DML_lookup_site_rank(ranklist, sites->subset_io_sites, r);
    if(s >= 0)
      sites->subset_rank[s] = t;
    t++;
  } while(DML_lex_next(&dim, coords, latdim, lower, ubound));
	  
  free(coords);
  free(ubound);
  free(ranklist);

  return 0;
}

/*------------------------------------------------------------------*/
/* Create the subset rank list for serial reading (single/part/multifile). */
/* See DML_create_subset_rank below for a definition of this list   */

/* Return value 0 for success and 1 for malloc failure */
int DML_create_subset_rank_serial(DML_SiteList *sites, DML_Layout *layout){

  DML_SiteRank r;
  int s;

  sites->use_subset = 1;
  sites->subset_rank = 
    (int *)malloc(sizeof(int)*sites->number_of_io_sites);  /* Could be less */
  if(sites->subset_rank == NULL)return 1;
  r = DML_init_site_loop(sites);
  s = 0;
  do {
    if(DML_rank_inside_subset(r, layout))
      sites->subset_rank[sites->current_index] = s++;
    else
      sites->subset_rank[sites->current_index] = -1;
  } while(DML_next_site(&r, sites));

  sites->subset_io_sites = s;

  return 0;
}

/*------------------------------------------------------------------*/
/* Create the subset rank list                                      */
/* This list maps the sites handled by this I/O partition to the
   positions of their data in the file record, relative to the
   beginning of the record.
   There is one entry in this list for each site in the I/O partition
   to which this node belongs.  The order of entries is lexicographic
   by site coordinate.   */
/* Return value 0 for success and 1 for malloc failure */
int DML_create_subset_rank(DML_SiteList *sites, DML_Layout *layout,
			   int volfmt, int serpar){
  sites->use_subset = 0;
  sites->subset_io_sites = sites->number_of_io_sites;
  if(layout->recordtype == DML_FIELD)return 0;

  if(volfmt == DML_SINGLEFILE && serpar == DML_PARALLEL)
    return DML_create_subset_rank_parallel(sites, layout);
  else
    return DML_create_subset_rank_serial(sites, layout);
}

/*------------------------------------------------------------------*/
void DML_destroy_subset_rank(DML_SiteList *sites){
  if(sites->use_subset == 1)
    if(sites->subset_rank != NULL)
      free(sites->subset_rank);
}

/*------------------------------------------------------------------*/
/* Checksum "class" */
/* We do a crc32 sum on the site data -- then do two lexicographic-
   rank-based bit rotations and XORs on the resulting crc32
   checksum */

/* Initialize checksums */
void DML_checksum_init(DML_Checksum *checksum){
  checksum->suma = 0;
  checksum->sumb = 0;
}

/* Accumulate checksums */
void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank, 
			char *buf, size_t size){

  DML_SiteRank rank29 = rank;
  DML_SiteRank rank31 = rank;
  uint32_t work = DML_crc32(0, (unsigned char*)buf, size);

  rank29 %= 29; rank31 %= 31;

  checksum->suma ^= work<<rank29 | work>>(32-rank29);
  checksum->sumb ^= work<<rank31 | work>>(32-rank31);
}

/* Combine checksums over all nodes */
void DML_checksum_combine(DML_Checksum *checksum){
  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);
}


/* Add single checksum set to the total */
void DML_checksum_peq(DML_Checksum *total, DML_Checksum *checksum){
  total->suma ^= checksum->suma;
  total->sumb ^= checksum->sumb;
}

/*------------------------------------------------------------------*/
/* Is this a big endian architecture? Return 1 or 0. */
int DML_big_endian(void)
{
  union {
    int  l;
    char c[sizeof(int)];
  } u;
  u.l = 1;

  return (u.c[sizeof(int)-1] == 1 ? 1 : 0); 
}


/* Do byte reversal on n contiguous 32-bit words */
void DML_byterevn32(uint32_t w[], size_t n)
{
  uint32_t old,newv;
  size_t j;

  assert(sizeof(uint32_t) == 4);
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      newv = old >> 24 & 0x000000ff;
      newv |= old >> 8 & 0x0000ff00;
      newv |= old << 8 & 0x00ff0000;
      newv |= old << 24 & 0xff000000;
      w[j] = newv;
    }
}


/* Do byte reversal on n contiguous 64-bit words */
void DML_byterevn64(uint32_t w[], size_t n)
{
  uint32_t tmp;
  size_t j;

  assert(sizeof(uint32_t) == 4);
  
  /* First swap pairs of 32-bit words */
  for(j=0; j<n; j++){
    tmp = w[2*j];
    w[2*j] = w[2*j+1];
    w[2*j+1] = tmp;
  }

  /* Then swap bytes in 32-bit words */
  DML_byterevn32(w, 2*n);
}

/* Do byte reversal on size bytes of contiguous words,
   each word consisting of word_size bytes 
   word_size = 1, 4 or 8 are the only choices. */

void DML_byterevn(char *buf, size_t size, int word_size)
{
  if(word_size == 1) {
    /* NOP */
  }
  else if(word_size == 4) {
    DML_byterevn32((uint32_t *)buf, size/word_size);
  }
  else if(word_size == 8) {
    DML_byterevn64((uint32_t *)buf, size/word_size);
  }
  else{
    printf("DML_byterevn: illegal word_size %d\n",word_size);
  }
}

/*------------------------------------------------------------------*/
/* Read and write buffer management */

/* Compute number of sites worth of data that fit in allowed space */
/* The number is supposed to be a multiple of "factor" */
size_t DML_max_buf_sites(size_t size, int factor){
  return ((size_t)DML_BUF_BYTES/(size*factor))*factor;
}

/*------------------------------------------------------------------*/
char *
DML_allocate_buf(size_t size, size_t *max_buf_sites)
{
  char *lbuf = NULL;
  while(*max_buf_sites>0) {
    lbuf = (char*) malloc(*max_buf_sites*size);
    if(lbuf!=NULL) break;
    *max_buf_sites /= 2;
  }
  return lbuf;
}

/*------------------------------------------------------------------*/
/* Write buffer to the current position in the file */

int DML_write_buf_current(LRL_RecordWriter *lrl_record_out, 
			  char *lbuf, size_t buf_sites, size_t size, 
			  uint64_t *nbytes, char *myname, 
			  int this_node){

  if(LRL_write_bytes(lrl_record_out, lbuf, buf_sites*size)
     != buf_sites*size){
    printf("%s(%d) write error\n",myname,this_node);
    return 1;
  }
  *nbytes += buf_sites*size;

  return 0;
}

/*------------------------------------------------------------------*/
/* Write buffer to a seek position in the file */

/* If the outbuf contains data for more than one site, we assume the
   sites are in the correct order for writing to the requested point
   in the file */

int DML_write_buf_seek(LRL_RecordWriter *lrl_record_out, 
		       DML_SiteRank seeksite, 
		       char *lbuf, size_t buf_sites, size_t size,
		       uint64_t *nbytes, char *myname, int this_node){

  /* Seek to the appropriate position */
  if(LRL_seek_write_record(lrl_record_out,(off_t)size*seeksite)
     != LRL_SUCCESS){
    printf("%s(%d) seek error\n",myname,this_node);
    return 1;
  }

  /* Then write */
  return DML_write_buf_current(lrl_record_out, lbuf, buf_sites, size,
			       nbytes, myname, this_node);
}

/*------------------------------------------------------------------*/
/* Seek if requested and read buffer. */
int
DML_read_buf(LRL_RecordReader *lrl_record_in, char *buf,
	     DML_SiteRank firstrank, size_t size, int num, int doseek)
{
  //printf("node %i firstrank %i size %li num %i doseek %i\n", QMP_get_node_number(), firstrank, size, num, doseek);
  if(doseek) {
    if(LRL_seek_read_record(lrl_record_in,(off_t)size*firstrank)
       != LRL_SUCCESS) {
      return -1;
    }
  }
  size *= num;
  if(LRL_read_bytes(lrl_record_in, buf, size) != size) {
    return -1;
  }
  return 0;
}

/*------------------------------------------------------------------*/
/* Seek and read buffer.  We work with only one site at a time for now */


size_t DML_read_buf_seek(LRL_RecordReader *lrl_record_in, 
			 DML_SiteRank seeksite, size_t size,
			 char *lbuf, size_t *buf_extract, size_t buf_sites, 
			 size_t max_buf_sites, size_t isite, 
			 size_t max_send_sites, 
			 uint64_t *nbytes, char *myname, int this_node,
			 int *err){
  /* Number of available sites in read buffer */
  size_t new_buf_sites = buf_sites;   

  /* Seeking makes buffering more complicated.  We take the easy way
     out until we are forced to find an intelligent way to do this */

  max_buf_sites = 1;  /* Force a max one-site buffer. Should be
			 changed in the future. */
  *err = 0;

  if(*buf_extract == buf_sites){  
    /* new buffer length  = remaining sites, but never bigger 
       than the buffer capacity */
    new_buf_sites = max_send_sites - isite;
    if(new_buf_sites > max_buf_sites) new_buf_sites = max_buf_sites; 

    /* Seek to the appropriate position in the record */
    if(LRL_seek_read_record(lrl_record_in,(off_t)size*seeksite)
       != LRL_SUCCESS){
      *err = -1;
      return 0;
    }

    /* Fill the buffer */
    if( LRL_read_bytes(lrl_record_in, lbuf, new_buf_sites*size) 
	!= new_buf_sites*size){
      printf("%s(%d) read error\n", myname,this_node); 
      *err = -1;
      return 0;
    }
    *nbytes += new_buf_sites*size;
    *buf_extract = 0;  /* reset counter */
  }  /* end of the buffer read */

  return new_buf_sites;
}

/*------------------------------------------------------------------*/
/* Get new buffer if no data remains to be processed */

size_t DML_read_buf_next(LRL_RecordReader *lrl_record_in, size_t size,
			 char *lbuf, size_t *buf_extract, size_t buf_sites, 
			 size_t max_buf_sites, size_t isite, 
			 size_t max_send_sites, 
			 uint64_t *nbytes, char *myname, int this_node,
			 int *err){
  /* Number of available sites in read buffer */
  size_t new_buf_sites = buf_sites;   

  *err = 0;

  if(*buf_extract == buf_sites){  
    /* new buffer length  = remaining sites, but never bigger 
       than the buffer capacity */
    new_buf_sites = max_send_sites - isite;
    if(new_buf_sites > max_buf_sites) new_buf_sites = max_buf_sites; 
    /* Fill the buffer */
    if( LRL_read_bytes(lrl_record_in, lbuf, new_buf_sites*size) 
	!= new_buf_sites*size){
      printf("%s(%d) read error\n", myname,this_node); 
      *err = -1;
      return 0;
    }
    *nbytes += new_buf_sites*size;
    *buf_extract = 0;  /* reset counter */
  }  /* end of the buffer read */

  return new_buf_sites;
}

/*------------------------------------------------------------------*/
/* Determine the node that does my I/O */
int DML_my_ionode(int volfmt, int serpar, DML_Layout *layout){

  if(volfmt == DML_SINGLEFILE){
    if(serpar == DML_SERIAL)
      return layout->master_io_node;
    else
      return layout->ionode(layout->this_node);
  }
  else if(volfmt == DML_MULTIFILE){
    return layout->this_node;
  }
  else if(volfmt == DML_PARTFILE){
    return layout->ionode(layout->this_node);
  }
  else {
    printf("DML_my_ionode: Bad volfmt code %d\n",volfmt);
    return 0;
  }
}

/*------------------------------------------------------------------*/
/* Synchronize the writers */

int DML_synchronize_out(LRL_RecordWriter *lrl_record_out, DML_Layout *layout){
  void *state_ptr;
  size_t state_size;
  int master_io_node = layout->master_io_node;

  /* DML isn't supposed to know the inner workings of LRL or LIME,
     so the state is captured as a string of bytes that only LRL
     understands.  All nodes create their state structures. */
  LRL_get_writer_state(lrl_record_out, &state_ptr, &state_size);

  /* The broadcast assumes all nodes are in the synchronization group 
     which is what we want for singlefile parallel I/O.
     If we decide later to do partfile parallel I/O we will need to
     change it */
  DML_broadcast_bytes((char *)state_ptr, state_size, layout->this_node, 
		      master_io_node);

  /* All nodes but the master node set their states */
  if(layout->this_node != master_io_node)
    LRL_set_writer_state(lrl_record_out, state_ptr);

  LRL_destroy_writer_state_copy(state_ptr);

  return 0;
}

/*------------------------------------------------------------------*/

/* The following four procedures duplicate (for the most part) the
   functionality of DML_partition_out.  They were broken out to allow
   finer high-level control of record writing for the file format
   conversion utilities.  Those utilities are currently run on a
   single processor.  Note that DML_partition_out now buffers
   messages, but these fragments have not yet been updated for that
   purpose. */

/* See DML_partition_out below for a description */
/* This routine opens a record and prepares to write a field */

DML_RecordWriter *DML_partition_open_out(
	   LRL_RecordWriter *lrl_record_out, size_t size, 
	   size_t set_buf_sites, DML_Layout *layout, DML_SiteList *sites,
	   int volfmt, int serpar, DML_Checksum *checksum)
{
  char *outbuf;
  int *coords;
  int this_node = layout->this_node;
  int my_io_node;
  int latdim = layout->latdim;
  size_t max_buf_sites;
  DML_RecordWriter *dml_record_out;
  char myname[] = "DML_partition_open_out";

  dml_record_out = (DML_RecordWriter *)malloc(sizeof(DML_RecordWriter));
  if(!dml_record_out){
    printf("%s(%d): No space for DML_RecordWriter\n",myname,this_node);
    return NULL;
  }

  /* Get my I/O node */
  my_io_node = DML_my_ionode(volfmt, serpar, layout);

  /* Allocate buffer for writing or sending data */
  /* I/O node needs a large buffer.  Others only enough for one site */
  /* But if set_buf_sites > 0 it determines the size of the buffer */
  if(this_node == my_io_node){
    max_buf_sites = DML_max_buf_sites(size,1);
    if(set_buf_sites > 0)max_buf_sites = set_buf_sites;
  }
  else{
    max_buf_sites = 1;
  }

  outbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!outbuf){
    printf("%s(%d) can't malloc outbuf\n",myname,this_node);
    return 0;
  }

  /* Allocate lattice coordinate */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords){free(outbuf);return NULL;}
  
  /* Initialize checksum */
  DML_checksum_init(checksum);
  
  /* Save/set the initial state */

  dml_record_out->lrl_rw          = lrl_record_out;
  dml_record_out->outbuf          = outbuf;
  dml_record_out->buf             = outbuf;
  dml_record_out->coords          = coords;
  dml_record_out->checksum        = checksum;
  dml_record_out->current_node    = my_io_node;
  dml_record_out->my_io_node      = my_io_node;
  dml_record_out->nbytes          = 0;
  dml_record_out->isite           = 0;
  dml_record_out->buf_sites       = 0;
  dml_record_out->max_buf_sites   = max_buf_sites;
  dml_record_out->max_dest_sites  = sites->subset_io_sites;

  return dml_record_out;
}

/*------------------------------------------------------------------*/
/* See DML_partition_out below for a description */
/* This procedure writes one site's worth of data with lexicographic
   rank snd_coords to an open record at a position specified by the
   site rank "seeksite" relative to the beginning of the binary
   data for that record */

int DML_partition_subset_sitedata_out(DML_RecordWriter *dml_record_out,
	   void (*get)(char *buf, size_t index, int count, void *arg),
           DML_SiteRank subset_rank, DML_SiteRank snd_coords, int count, 
	   size_t size, int word_size, void *arg, DML_Layout *layout)
{

  LRL_RecordWriter *lrl_record_out = dml_record_out->lrl_rw;
  char *outbuf                     = dml_record_out->outbuf;
  char *buf                        = dml_record_out->buf;
  int *coords                      = dml_record_out->coords;
  DML_Checksum *checksum           = dml_record_out->checksum;
  int current_node	           = dml_record_out->current_node;
  int my_io_node                   = dml_record_out->my_io_node;
  uint64_t nbytes                  = dml_record_out->nbytes;
  size_t isite                     = dml_record_out->isite;
  size_t buf_sites                 = dml_record_out->buf_sites;
  size_t max_buf_sites             = dml_record_out->max_buf_sites;
  size_t max_dest_sites            = dml_record_out->max_dest_sites;

  int this_node    = layout->this_node;
  int latdim       = layout->latdim;
  int *latsize     = layout->latsize;

  int new_node;
  int status;
  char scratch_buf[4];
  char myname[] = "DML_partition_subset_sitedata_out";

  scratch_buf[0] = scratch_buf[1] = scratch_buf[2] = scratch_buf[3] = '\0';

  /* Convert lexicographic rank to coordinates */
  DML_lex_coords(coords, latdim, latsize, snd_coords);
  
  /* Node that has this data and sends it to my_io_node */
  new_node = layout->node_number(coords);
  
  /* Send nodes must wait for a ready signal from the I/O node
     to prevent message pileups on the I/O node */
  
  /* CTS only if changing data source node */
  if(new_node != current_node){
    DML_clear_to_send(scratch_buf,4,my_io_node,new_node);
    current_node = new_node;
  }
  
  /* Fetch site data and copy to the write buffer */
  if(this_node == current_node){
    /* Fetch directly to the buffer */
    buf = outbuf + size*buf_sites;
    get(buf,layout->node_index(coords),count,arg);
    buf_sites++;
  }
  
  /* Send result to my I/O node. Avoid I/O node sending to itself. */
  if (current_node != my_io_node) 
    {
#if 1
      /* Data from any other node is received in the I/O node write buffer */
      if(this_node == my_io_node){
	buf = outbuf + size*buf_sites;
      }
      DML_route_bytes(buf,size,current_node,my_io_node);
      if(this_node == my_io_node)buf_sites++;
      if(this_node == current_node)buf_sites--;
#else
      /* Data from any other node is received in the I/O node write buffer */
      if(this_node == my_io_node){
	buf = outbuf + size*buf_sites;
	DML_get_bytes(buf,size,current_node);
	buf_sites++;
      }
      
      /* All other nodes send the data */
      if(this_node == current_node){
	DML_send_bytes(buf,size,my_io_node);
	buf_sites--;
      }
#endif
    }
  
  /* my_io_node writes the data */
  if(this_node == my_io_node)
    {
      /* Do byte reordering before checksum */
      if (! DML_big_endian())
	DML_byterevn(buf, size, word_size);
      
      /* Update checksum */
      DML_checksum_accum(checksum, snd_coords, buf, size);
      
      /* Write the buffer when full */

      if( (buf_sites >= max_buf_sites) || (isite >= max_dest_sites-1) )
	{
	  status = DML_write_buf_seek(lrl_record_out, subset_rank,
				      outbuf, buf_sites, size, &nbytes,
				      myname, this_node);
	  buf_sites = 0;
	  if(status !=  0) {return 1;}
	}
    }
  isite++;

  /* Save changes to state */

  dml_record_out->current_node    = current_node;
  dml_record_out->nbytes          = nbytes;
  dml_record_out->isite           = isite;
  dml_record_out->buf_sites       = buf_sites;

  return 0;
}
  
/*------------------------------------------------------------------*/
/* Write one site's worth of data with lexicographic rank snd_coords */

int DML_partition_sitedata_out(DML_RecordWriter *dml_record_out,
	   void (*get)(char *buf, size_t index, int count, void *arg),
           DML_SiteRank snd_coords, int count, size_t size, int word_size, 
           void *arg, DML_Layout *layout, DML_SiteList *sites){
  DML_SiteRank subset_rank;
  int this_node = layout->this_node;
  char myname[] = "DML_partition_sitedata_out";

  /* Convert lexicographic index for this partition to the physical
     rank in the record */
  if(DML_lookup_subset_rank(&subset_rank, snd_coords, sites)!= 0){
    printf("%s(%d) Request to write a site %d not planned for the record.\n",
	   myname, this_node, snd_coords);
    return 1;
  }
  if(DML_partition_subset_sitedata_out(dml_record_out, get, subset_rank, 
            snd_coords, count, size, word_size, arg, layout) != 0)
    return 1;

  return 0;
}

/*------------------------------------------------------------------*/
/* See DML_partition_out below for a description */
/* This routine closes an open record and cleans up */

/* See DML_partition_out below for a description */
/* This writes one site's worth of data to an open record */

int DML_partition_allsitedata_out(DML_RecordWriter *dml_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
           DML_Layout *layout, DML_SiteList *sites)
{

  DML_SiteRank snd_coords;

  /* Iterate over all sites processed by this I/O partition */

  if(DML_init_subset_site_loop(&snd_coords, sites)==0)
    return 0;
  do {
    if(DML_partition_sitedata_out(dml_record_out, get, snd_coords, 
		  count, size, word_size, arg, layout, sites) != 0)
      return 1;
    
  } while(DML_next_subset_site(&snd_coords, sites));
  
  return 0;
}

  
/*------------------------------------------------------------------*/
/* See DML_partition_out below for a description */
/* This routine closes an open record and cleans up */

uint64_t DML_partition_close_out(DML_RecordWriter *dml_record_out)
{

  uint64_t nbytes = dml_record_out->nbytes;

  if(dml_record_out == NULL)return 0;
  if(dml_record_out->coords != NULL)
    free(dml_record_out->coords);
  if(dml_record_out->outbuf != NULL)
    free(dml_record_out->outbuf);
  free(dml_record_out);

  /* Number of bytes written by this node only */
  return nbytes;
}

/*------------------------------------------------------------------*/
/* Flush the outputbuffer to the file */
static int DML_flush_outbuf(LRL_RecordWriter *lrl_record_out, int serpar,
			    DML_SiteRank snd_coords, 
			    char *outbuf, size_t buf_sites, size_t size,
			    uint64_t *nbytes,  int this_node)
{
  char myname[] = "DML_flush_outbuf";
  int status;

  if(buf_sites == 0)return 0;

  if(serpar == DML_SERIAL)
    status = DML_write_buf_current(lrl_record_out,
				   outbuf, buf_sites, size, nbytes,  
				   myname, this_node);
  else
    status = DML_write_buf_seek(lrl_record_out, snd_coords, 
				outbuf, buf_sites, size, nbytes,
				myname, this_node);
  return status;
}

#if defined(QIO_USE_DML_OUT_BUFFERING)
/*------------------------------------------------------------------*/
/* Flush message buffer to IO buffer.  Do byte reordering if needed.
   Accumulate checksums. (tbuf_sites is not reset here)
 */

static void DML_flush_tbuf_to_outbuf(size_t size, 
				    char *outbuf, size_t buf_sites, 
				    char *tbuf, size_t tbuf_sites )
{
  if(tbuf_sites == 0)return;

  /* Copy tbuf to outbuf */
  memcpy((void *)(outbuf + size*buf_sites), 
	 (void *)tbuf, size*tbuf_sites);
}

/*------------------------------------------------------------------*/
/* Each I/O node (or the master node) receives data from all of its
   nodes and writes it to its file.
   Returns the checksum and number of bytes written by this node only */

/* In order to be nondeadlocking, this algorithm requires that the set
   of nodes containing sites belonging to any single I/O node are
   disjoint from the corresponding set for any other I/O node.  This
   algorithm is intended for SINGLEFILE/SERIAL, MULTIFILE, and
   PARTFILE modes. */

uint64_t DML_partition_out(LRL_RecordWriter *lrl_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
	   DML_Layout *layout, DML_SiteList *sites, int volfmt, 
	   int serpar, DML_Checksum *checksum)
{
  double dtall=0, dtwrite=0, dtsend=0;
  char *buf,*outbuf = NULL,*tbuf = NULL, *scratch_buf;
  int current_node, new_node;
  int *coords;
  int this_node = layout->this_node;
  int my_io_node;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  size_t isite,buf_sites,tbuf_sites,max_buf_sites=0,max_tbuf_sites;
  int status;
  DML_SiteRank subset_rank;
  DML_SiteRank snd_coords,prev_coords,outbuf_coords;
  uint64_t nbytes = 0;
  char myname[] = "DML_partition_out";

  timestart(dtall);

  /* Get my I/O node */
  my_io_node = DML_my_ionode(volfmt, serpar, layout);

  /* Allocate buffers for writing and sending data */

  /* All nodes need a temporary message buffer for holding some number
     of lexicographically contiguous sites.  */

  max_tbuf_sites = DML_TBUF_BYTES/size;
  if(max_tbuf_sites<1) max_tbuf_sites = 1;

  /* I/O node needs a large output buffer. */
  if(this_node == my_io_node){
    max_buf_sites = DML_max_buf_sites(size,max_tbuf_sites);
    /* This buffer can't be smaller than the message buffer */
    if(max_tbuf_sites > max_buf_sites)max_buf_sites = max_tbuf_sites;
  }

#if 0
  /* For parallel I/O we don't try to buffer for messaging.  When each
     node can do I/O the data being written is local.  If we start
     doing PARTFILE parallel I/O we may want to buffer. */
  if(serpar == DML_PARALLEL){
    max_buf_sites = 1;
    max_tbuf_sites = 1;
  }
#endif

  /* Only the I/O node has an output buffer */
  if(this_node == my_io_node){
    outbuf = DML_allocate_buf(size, &max_buf_sites);
    if(!outbuf){
      printf("%s(%d) can't malloc outbuf\n",myname,this_node);
      return 0;
    }
  }

  tbuf = DML_allocate_buf(size, &max_tbuf_sites);
  if(!tbuf){
    printf("%s(%d) can't malloc tbuf\n",myname,this_node);
    free(outbuf);
    return 0;
  }

  /* Scratch for clear-to-send signal */
  { size_t one=1; scratch_buf = DML_allocate_buf(4, &one); }
  if(!scratch_buf){
    printf("%s(%d) can't malloc scratch_buf\n",myname,this_node);
    free(outbuf); free(tbuf);
    return 0;
  }
  memset(scratch_buf,0,4);

  /* Allocate lattice coordinate */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords){free(outbuf);free(tbuf);free(scratch_buf);return 0;}
  
  /* Initialize checksum */
  DML_checksum_init(checksum);
  
#ifdef DML_DEBUG
  if (! DML_big_endian())
    printf("%s(%d): byte reversing %d\n",myname,this_node,word_size);
#endif
  
  /* Maximum number of sites to be processed */
  //max_dest_sites = sites->subset_io_sites;
  isite = 0;  /* Running count of all sites */
  
  current_node = my_io_node;
  
  /* Loop over the sending coordinates */
  buf = outbuf;
  buf_sites = 0;   /* Count of sites in the output buffer */
  tbuf_sites = 0;  /* Count of sites in the message tbuf */

  if(DML_init_subset_site_loop(&snd_coords, sites) == 0){
    free(outbuf); free(tbuf); free(scratch_buf); free(coords); 
    return 0;
  }

  /* To track lexicographic order */
  prev_coords = snd_coords-1;
  /* Remember coordinate of first site in outbuf so we know where to
     seek */
  outbuf_coords = snd_coords;
  /* The subset_rank locates the datum for outbuf_coords in the
     record our I/O partition is writing */
  subset_rank = (DML_SiteRank)DML_subset_rank(outbuf_coords, sites);

  do {
    /* Convert lexicographic rank to coordinates */
    DML_lex_coords(coords, latdim, latsize, snd_coords);

    /* Node that sends data */
    new_node = layout->node_number(coords);

    /* A node sends its message buffer to the io_node when changing
       nodes or when its message buffer is full or when the
       lexicographic order is broken */
    if(new_node != current_node || 
       tbuf_sites >= max_tbuf_sites ||
       snd_coords != prev_coords + 1){
      if(tbuf_sites > 0){
	/* Node with data sends its message buffer to the I/O node's tbuf */
	if(current_node != my_io_node) {
	  timestart(dtsend);
	  DML_route_bytes(tbuf,size*tbuf_sites,current_node,my_io_node);
	  timestop(dtsend);
	}
	/* The I/O node flushes its tbuf and accumulates the checksum */
	if(this_node == my_io_node){
	  DML_flush_tbuf_to_outbuf(size, outbuf, buf_sites, tbuf, tbuf_sites);
	  buf_sites += tbuf_sites;
	  /* The I/O node writes the I/O buffer when full or when the
	     lexicographic order is broken */
	  if(buf_sites > max_buf_sites - max_tbuf_sites ||
	     snd_coords != prev_coords + 1){
	    timestart(dtwrite);
	    status = DML_flush_outbuf(lrl_record_out, serpar, subset_rank, 
			     outbuf, buf_sites, size, &nbytes, this_node);
	    timestop(dtwrite);
	    buf_sites = 0;
	    outbuf_coords = snd_coords;
	    subset_rank = (DML_SiteRank)DML_subset_rank(outbuf_coords, sites);
	    if(status != 0) {
	      free(outbuf); free(tbuf); free(scratch_buf); free(coords); 
	      return 0;
	    }
	  }
	}
	tbuf_sites = 0;
      }
      
      /* Send nodes must wait for a ready signal from the I/O node
	 to prevent message pileups on the I/O node */
      
      /* CTS only if changing data source node */
      if(new_node != current_node){
	timestart(dtsend);
	DML_clear_to_send(scratch_buf,4,my_io_node,new_node);
	timestop(dtsend);
	current_node = new_node;
      }
    } /* current_node != newnode || tbuf_sites >= max_tbuf_sites */

    /* The node with the data just appends it to its message buffer */
    if(this_node == current_node){
      /* Fetch to the message buffer */
      buf = tbuf + size*tbuf_sites;
      get(buf,layout->node_index(coords),count,arg);

      /* Do byte reordering and update checksum */
      if (! DML_big_endian()) DML_byterevn(buf, size, word_size);
      DML_checksum_accum(checksum, snd_coords, buf, size);
    }

    /* The I/O node and current node count sites together */
    if(this_node == current_node || this_node == my_io_node)tbuf_sites++;

    isite++;
    prev_coords = snd_coords;
  } while(DML_next_subset_site(&snd_coords, sites));

  /* Purge any remaining data */

  if(tbuf_sites > 0){
    if(current_node != my_io_node) {
      timestart(dtsend);
      DML_route_bytes(tbuf,size*tbuf_sites,current_node,my_io_node);
      timestop(dtsend);
    }
  }

  if(this_node == my_io_node){
    DML_flush_tbuf_to_outbuf(size, outbuf, buf_sites, tbuf, tbuf_sites);
    buf_sites += tbuf_sites;
    tbuf_sites = 0;
    timestart(dtwrite);
    status = DML_flush_outbuf(lrl_record_out, serpar, subset_rank, 
			      outbuf, buf_sites, size, &nbytes, this_node);
    timestop(dtwrite);
    buf_sites = 0;
    if(status !=  0) nbytes = 0;
  }

  free(coords);
  free(scratch_buf);
  free(outbuf);
  free(tbuf);

  timestop(dtall);
  if(QIO_verbosity()>=QIO_VERB_LOW && this_node==layout->master_io_node) {
    printf("%s times: write %.2f  send %.2f  total %.2f\n",
	   myname, dtwrite, dtsend, dtall);
  }
  /* Number of bytes written by this node only */
  return nbytes;
}

/*------------------------------------------------------------------*/
#else  /* not defined(QIO_USE_DML_OUT_BUFFERING) */

/* Each I/O node (or the master node) receives data from all of its
   nodes and writes it to its file.
   Returns the checksum and number of bytes written by this node only */

/* In order to be nondeadlocking, this algorithm requires that the set
   of nodes containing sites belonging to any single I/O node are
   disjoint from the corresponding set for any other I/O node.  This
   algorithm is intended for SINGLEFILE/SERIAL, MULTIFILE, and
   PARTFILE modes. */

/* This is the old algorithm that sent only one site's worth at a time */

uint64_t DML_partition_out(LRL_RecordWriter *lrl_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
	   DML_Layout *layout, DML_SiteList *sites, int volfmt, 
	   int serpar, DML_Checksum *checksum)
{
  char *buf,*outbuf,*scratch_buf;
  int current_node, new_node;
  int *coords;
  int this_node = layout->this_node;
  int my_io_node;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  size_t isite,buf_sites,max_buf_sites,max_dest_sites;
  int status;
  DML_SiteRank snd_coords, subset_rank;
  uint64_t nbytes = 0;
  char myname[] = "DML_partition_out";

  /* Get my I/O node */
  my_io_node = DML_my_ionode(volfmt, serpar, layout);

  /* Allocate buffer for writing or sending data */
  /* I/O node needs a large buffer.  Others only enough for one site */
  if(this_node == my_io_node)
    max_buf_sites = DML_max_buf_sites(size,1);
  else
    max_buf_sites = 1;

  if(serpar == DML_PARALLEL)
    max_buf_sites = 1;

  outbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!outbuf){
    printf("%s(%d) can't malloc outbuf\n",myname,this_node);
    return 0;
  }

  { size_t one=1; scratch_buf = DML_allocate_buf(4, &one); }
  if(!scratch_buf){
    printf("%s(%d) can't malloc scratch_buf\n",myname,this_node);
    free(outbuf);
    return 0;
  }
  memset(scratch_buf,0,4);

  /* Allocate lattice coordinate */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords){free(outbuf);return 0;}
  
  /* Initialize checksum */
  DML_checksum_init(checksum);
  
#ifdef DML_DEBUG
  if (! DML_big_endian())
    printf("%s(%d): byte reversing %d\n",myname,this_node,word_size);
#endif
  
  /* Maximum number of sites to be processed */
  max_dest_sites = sites->subset_io_sites;
  isite = 0;  /* Running count of all sites */
  
  current_node = my_io_node;
  
  /* Loop over the sending coordinates */
  buf = outbuf;
  buf_sites = 0;   /* Count of sites in the output buffer */
  if(DML_init_subset_site_loop(&snd_coords, sites) == 0){
    free(outbuf); free(coords); free(scratch_buf);
    return 0;
  }

  do {
    /* Convert lexicographic rank to coordinates */
    DML_lex_coords(coords, latdim, latsize, snd_coords);

    /* Node that sends data */
    new_node = layout->node_number(coords);

    /* Send nodes must wait for a ready signal from the I/O node
       to prevent message pileups on the I/O node */

    /* CTS only if changing data source node */
    if(new_node != current_node){
      DML_clear_to_send(scratch_buf,4,my_io_node,new_node);
      current_node = new_node;
    }

    /* Copy to the write buffer */
    if(this_node == current_node){
      /* Fetch directly to the buffer */
      buf = outbuf + size*buf_sites;
      get(buf,layout->node_index(coords),count,arg);
      buf_sites++;
    }

    /* Send result to my I/O node. Avoid I/O node sending to itself. */
    if (current_node != my_io_node) 
    {
#if 1
      /* Data from any other node is received in the I/O node write buffer */
      if(this_node == my_io_node){
	buf = outbuf + size*buf_sites;
      }
      DML_route_bytes(buf,size,current_node,my_io_node);
      if(this_node == my_io_node)buf_sites++;
      if(this_node == current_node)buf_sites--;
#else
      /* Data from any other node is received in the I/O node write buffer */
      if(this_node == my_io_node){
	buf = outbuf + size*buf_sites;
	DML_get_bytes(buf,size,current_node);
	buf_sites++;
      }
    
      /* All other nodes send the data */
      if(this_node == current_node){
	DML_send_bytes(buf,size,my_io_node);
	buf_sites--;
      }
#endif
    }

    /* Now write data */
    if(this_node == my_io_node)
    {
      /* Do byte reordering before checksum */
      if (! DML_big_endian())
	DML_byterevn(buf, size, word_size);
      
      /* Update checksum */
      DML_checksum_accum(checksum, snd_coords, buf, size);
      
      /* Write the buffer when full */

      if( (buf_sites >= max_buf_sites) || (isite == max_dest_sites-1) )
	{
	  /* The subset_rank locates the datum for snd_coords in the
	     record that our I/O partition is writing */
	  subset_rank = (DML_SiteRank)DML_subset_rank(snd_coords, sites);
	  status = DML_flush_outbuf(lrl_record_out, serpar, subset_rank,
				   outbuf, buf_sites, size, &nbytes, 
				   this_node);
	  buf_sites = 0;
	  if(status != 0) {free(outbuf); free(coords); return 0;}
	}
    }
    isite++;
  } while(DML_next_subset_site(&snd_coords, sites));
  
  free(coords);
  free(scratch_buf);
  free(outbuf);

  /* Number of bytes written by this node only */
  return nbytes;
}
#endif  /* if defined(QIO_USE_DML_OUT_BUFFERING) */


/*------------------------------------------------------------------*/
/* The master node fetches the global data in one call to "get"  and writes */
/* Returns the number of bytes written */

size_t DML_global_out(LRL_RecordWriter *lrl_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
           DML_Layout *layout, int volfmt, 
	   DML_Checksum *checksum)
{
  char *buf;
  int this_node = layout->this_node;
  size_t nbytes = 0;
  char myname[] = "DML_global_out";
  
  /* Allocate buffer for datum */
  buf = (char *)malloc(size);
  if(!buf){
    printf("%s(%d) can't malloc buf\n",myname,this_node);
    return 0;
  }
  
  /* Initialize checksum */
  DML_checksum_init(checksum);
  
#ifdef DML_DEBUG
  if (! DML_big_endian())
    printf("%s(%d): byte reversing %d\n",myname,this_node,word_size);
#endif
  
  /* Master node writes all the data */
  if(this_node == layout->master_io_node){
    /* Get all the data.  0 for the unused site index */
    get(buf,0,count,arg);
    
    /* Do byte reordering before checksum */
    if (! DML_big_endian())
      DML_byterevn(buf, size, word_size);
    
    /* Do checksum.  Straight crc32. */
    DML_checksum_accum(checksum, 0, buf, size);
    
    /* Write all the data */
    nbytes = LRL_write_bytes(lrl_record_out,(char *)buf,size);
    if( nbytes != size){
      free(buf); return 0;}
  }
  
  free(buf);
  return nbytes;
}

/*--------------------------------------------------------------------*/
/* THIS PROCEDURE IS OBSOLETE */
/* Each node writes its data to its own private file.  The order of
   sites is sequential according to the layout storage order, which
   is not likely to be lexicographic. */

/* Returns the number of bytes written by this node alone */

uint64_t DML_multifile_out(LRL_RecordWriter *lrl_record_out, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      int count, size_t size, int word_size, void *arg, 
	      DML_Layout *layout, DML_Checksum *checksum)
{
  
  size_t max_buf_sites, buf_sites;
  size_t isite, max_dest_sites;
  uint64_t nbytes = 0;
  DML_SiteRank rank;
  int this_node = layout->this_node;
  char myname[] = "DML_multifile_out";
  char *lbuf, *buf=NULL;
  int *coords;
  int status;

  /* Allocate buffer for writing */
  max_buf_sites = DML_max_buf_sites(size,1);
  lbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!lbuf){
    printf("%s(%d): Can't malloc lbuf\n",myname,this_node);
    return 0;
  }

  /* Allocate coordinate */
  coords = DML_allocate_coords(layout->latdim,myname,this_node);
  if(!coords){free(buf); return 0;}

  /* Initialize checksum */
  DML_checksum_init(checksum);

  buf_sites = 0;

  max_dest_sites = layout->sites_on_node;

  /** TODO: VECTORIZE THE TRANSFER - CD **/
  /* Loop over the storage order index of sites on the local node */
  for(isite = 0; isite < max_dest_sites; isite++){

    /* The coordinates of this site */
    layout->get_coords(coords, this_node, isite);

    /* The lexicographic rank of this site */
    rank = DML_lex_rank(coords, layout->latdim, layout->latsize);

    /* Fetch directly to the buffer */
    buf = lbuf + size*buf_sites;
    get(buf, isite, count, arg);
    buf_sites++;

    /* Accumulate checksums as the values are inserted into the buffer */
    
    /* Do byte reversal if needed */
    if (! DML_big_endian())
      DML_byterevn(buf, size, word_size);
    
    DML_checksum_accum(checksum, rank, buf, size);
    
    /* Write buffer when full or last site processed */
    if( (buf_sites >= max_buf_sites) || (isite == max_dest_sites - 1))
      {
	status = DML_write_buf_current(lrl_record_out, 
				       lbuf, buf_sites, size, &nbytes,
				       myname, this_node);
	buf_sites = 0;
	if(status != 0) {free(lbuf); free(coords); return 0;}
      }
  } /* isite */

  free(lbuf);   free(coords);
  
  /* Return the number of bytes written by this node only */
  return nbytes;

} /* DML_multifile_out */

/*------------------------------------------------------------------*/
/* Each node reads its data from its own private file.  The order of
   sites is assumed to be sequential according to the layout storage
   order. */

/* Returns the number of bytes read by this node alone */

uint64_t DML_multifile_in(LRL_RecordReader *lrl_record_in, 
	     DML_SiteRank sitelist[],
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     int count, size_t size, int word_size, void *arg, 
	     DML_Layout *layout, DML_Checksum *checksum)
{
  size_t buf_sites, buf_extract, max_buf_sites;
  size_t isite, max_send_sites;
  uint64_t nbytes = 0;
  DML_SiteRank rank;
  int this_node = layout->this_node;
  char myname[] = "DML_multifile_in";
  char *lbuf, *buf;
  int *coords;
  int err;

  /* Allocate buffer for reading */
  max_buf_sites = DML_max_buf_sites(size,1);
  lbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!lbuf)return 0;

  /* Allocate coordinate */
  coords = DML_allocate_coords(layout->latdim, myname, this_node);
  if(!coords){free(lbuf);return 0;}

  /* Initialize checksum */
  DML_checksum_init(checksum);

  buf_sites = 0;      /* Length of current read buffer */
  buf_extract = 0;    /* Counter for current site in read buffer */
  max_send_sites = layout->sites_on_node;

  /** TODO: VECTORIZE THE TRANSFER - CD **/
  /* Loop over the storage order site index for this node */
  for(isite = 0; isite < max_send_sites; isite++){

    /* The coordinates of this site */
    layout->get_coords(coords, this_node, isite);

    /* The lexicographic rank of this site */
    rank = DML_lex_rank(coords, layout->latdim, layout->latsize);

    /* Refill buffer if necessary */
    buf_sites = DML_read_buf_next(lrl_record_in, size, lbuf, 
				  &buf_extract, buf_sites, max_buf_sites, 
				  isite, max_send_sites, &nbytes, 
				  myname, this_node, &err);
    if(err < 0){free(lbuf);free(coords);return 0;}
    
    /* Copy data directly from the buffer */
    buf = lbuf + size*buf_extract;

    /* Accumulate checksums as the values are inserted into the buffer */
    DML_checksum_accum(checksum, rank, buf, size);
    
    /* Do byte reversal after checksum if needed */
    if (! DML_big_endian())
      DML_byterevn(buf, size, word_size);

    put(buf, isite, count, arg);
    
    buf_extract++;
  } /* isite */

  free(lbuf);   free(coords);
  
  /* Return the number of bytes read by this node only */
  return nbytes;
}

/*------------------------------------------------------------------*/
/* Synchronize the readers */

int DML_synchronize_in(LRL_RecordReader *lrl_record_in, DML_Layout *layout){
  void *state_ptr;
  size_t state_size;
  int master_io_node = layout->master_io_node;

  /* DML isn't supposed to know the inner workings of LRL or LIME,
     so the state is captured as a string of bytes that only LRL
     understands */
  LRL_get_reader_state(lrl_record_in, &state_ptr, &state_size);

  /* The broadcast assumes all nodes are in the synchronization group 
     which is what we want for singlefile parallel I/O.
     If we decide later to do partfile parallel I/O we will need to
     change it */
  DML_broadcast_bytes((char *)state_ptr, state_size, layout->this_node, 
		      master_io_node);
  
  /* All nodes but the master node set their states */
  if(layout->this_node != master_io_node)
    LRL_set_reader_state(lrl_record_in, state_ptr);

  LRL_destroy_reader_state_copy(state_ptr);

  return 0;
}

/*------------------------------------------------------------------*/
/* The following four procedures duplicate the functionality
   of DML_partition_in.  They were broken out to allow
   finer high-level control of record reading.  After they
   have been tested, they will replace DML_partition_in */

/* See DML_partition_in below for a description */
/* This routine opens a record and prepares to read a field */

DML_RecordReader *DML_partition_open_in(LRL_RecordReader *lrl_record_in, 
	  size_t size, size_t set_buf_sites, DML_Layout *layout, 
 	  DML_SiteList *sites, int volfmt, int serpar, DML_Checksum *checksum)
{
  char *inbuf;
  int my_io_node;
  int *coords;
  int this_node = layout->this_node;
  int latdim = layout->latdim;
  size_t max_buf_sites;
  DML_RecordReader *dml_record_in;
  char myname[] = "DML_partition_open_in";

  dml_record_in = (DML_RecordReader *)malloc(sizeof(DML_RecordReader));
  if(!dml_record_in){
    printf("%s(%d): No space for DML_RecordReader\n",myname,this_node);
    return NULL;
  }

  /* Get my I/O node */
  my_io_node = DML_my_ionode(volfmt, serpar, layout);

  /* Allocate buffer for reading or receiving data */
  /* I/O node needs a large buffer.  Others only enough for one site */
  if(this_node == my_io_node){
    max_buf_sites = DML_max_buf_sites(size,1);
    if(set_buf_sites > 0)max_buf_sites = set_buf_sites;
  }
  else
    max_buf_sites = 1;

 
  inbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!inbuf){
    printf("%s(%d) can't malloc inbuf\n",myname,this_node);
    return 0;
  }
 

  /* Allocate coordinate counter */
  coords = DML_allocate_coords(latdim, myname, this_node);
  if(!coords){free(inbuf); return 0;}
  
  /* Initialize checksum */
  DML_checksum_init(checksum);

#ifdef DML_DEBUG
  if (! DML_big_endian())
    printf("%s(%d): byte reversing %d\n",myname,this_node,word_size);
#endif

  /* Save/set the initial state */

  dml_record_in->lrl_rr          = lrl_record_in;
  dml_record_in->inbuf           = inbuf;
  dml_record_in->coords          = coords;
  dml_record_in->checksum        = checksum;
  dml_record_in->current_node    = my_io_node;
  dml_record_in->my_io_node      = my_io_node;
  dml_record_in->nbytes          = 0;
  dml_record_in->isite           = 0;
  dml_record_in->buf_sites       = 0;
  dml_record_in->buf_extract     = 0;
  dml_record_in->max_buf_sites   = max_buf_sites;
  dml_record_in->max_send_sites  = sites->subset_io_sites;

  return dml_record_in;
}

/*------------------------------------------------------------------*/
/* See DML_partition_in below for a description */
/* Seek within the record to a given linear index subset_rank 
   and read one site's worth of data */
int DML_partition_subset_sitedata_in(DML_RecordReader *dml_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
          DML_SiteRank subset_rank, DML_SiteRank rcv_coords, int count, 
          size_t size, int word_size,  void *arg, DML_Layout *layout)
{

  LRL_RecordReader *lrl_record_in  = dml_record_in->lrl_rr;
  char *inbuf                      = dml_record_in->inbuf;
  int *coords                      = dml_record_in->coords;
  DML_Checksum *checksum           = dml_record_in->checksum;
  int current_node	           = dml_record_in->current_node;
  int my_io_node                   = dml_record_in->my_io_node;
  uint64_t nbytes                  = dml_record_in->nbytes;
  size_t isite                     = dml_record_in->isite;
  size_t buf_sites                 = dml_record_in->buf_sites;
  size_t buf_extract               = dml_record_in->buf_extract;
  size_t max_buf_sites             = dml_record_in->max_buf_sites;
  size_t max_send_sites            = dml_record_in->max_send_sites;
  char *buf=NULL;

  int this_node = layout->this_node;
  int latdim    = layout->latdim;
  int *latsize  = layout->latsize;

  int dest_node;
  int err;
  char myname[] = "DML_partition_subset_sitedata_in";

  /* Convert lexicographic rank to coordinates */
  DML_lex_coords(coords, latdim, latsize, rcv_coords);
  
  /* The node that gets the next datum */
  dest_node = layout->node_number(coords);
  
  if(this_node == my_io_node){
    /* I/O node reads the next value */
    buf_sites = DML_read_buf_seek(lrl_record_in, subset_rank, size,
				  inbuf, &buf_extract, buf_sites,
				  max_buf_sites, isite, 
				  max_send_sites, &nbytes,
				  myname, this_node, &err);
    
    if(err < 0){
      printf("%s(%d) DML_read_buf_seek returns error\n",
	     myname,this_node);
      return 1;
    }

    /* Location of new datum on I/O node */
    buf = inbuf + size*buf_extract;
  }

  /* Send result to destination node. Avoid I/O node sending to itself. */
  if (dest_node != my_io_node) {
#if 1
    DML_route_bytes(buf,size,my_io_node,dest_node);
#else
    /* If destination elsewhere, send it */
    if(this_node == my_io_node){
      DML_send_bytes(buf, size, dest_node);
    }
    
    /* Other nodes receive from the master node */
    if(this_node == dest_node){
      DML_get_bytes(buf, size, my_io_node);
    }
#endif
  }
  
  /* Process data before inserting */
  if(this_node == dest_node){
    
    /* Accumulate checksum */
    DML_checksum_accum(checksum, rcv_coords, buf, size);
    
    /* Do byte reversal if necessary */
    if (! DML_big_endian())
      DML_byterevn(buf, size, word_size);
    
    /* Store the data */
    put(buf,layout->node_index(coords),count,arg);
  }
  
  buf_extract++;
  isite++;

  /* Save changes to state */

  dml_record_in->current_node    = current_node;
  dml_record_in->nbytes          = nbytes;
  dml_record_in->isite           = isite;
  dml_record_in->buf_sites       = buf_sites;
  dml_record_in->buf_extract     = buf_extract;

  return 0;
}

/*------------------------------------------------------------------*/
int DML_partition_sitedata_in(DML_RecordReader *dml_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  DML_SiteRank rcv_coords, int count, size_t size, int word_size, 
          void *arg, DML_Layout *layout, DML_SiteList *sites)
{
  DML_SiteRank subset_rank;
  int this_node = layout->this_node;
  char myname[] = "DML_partition_sitedata_in";

  /* Convert lexicographic index for this partition to the physical
     rank in the record */
  if(DML_lookup_subset_rank(&subset_rank, rcv_coords, sites)!= 0){
    printf("%s(%d) Request for a site %d not found in the record.\n",
	   myname, this_node, rcv_coords);
  return 1;
  }
  if(DML_partition_subset_sitedata_in(dml_record_in, put, subset_rank,
	      rcv_coords, count, size, word_size, arg, layout) != 0)
    return 1;

  return 0;
}

/*------------------------------------------------------------------*/
/* See DML_partition_in below for a description */
/* This routine reads all sites in the record */

int DML_partition_allsitedata_in(DML_RecordReader *dml_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  int count, size_t size, int word_size, void *arg, 
	  DML_Layout *layout, DML_SiteList *sites, int volfmt,
	  DML_Checksum *checksum)
{

  DML_SiteRank rcv_coords;

  if(DML_init_subset_site_loop(&rcv_coords, sites) == 0)
    return 0;
  do {
    if(DML_partition_sitedata_in(dml_record_in, put, rcv_coords,
		count, size, word_size, arg, layout, sites) != 0)
      return 1;

  }  while(DML_next_subset_site(&rcv_coords, sites));

  return 0;
}

/*------------------------------------------------------------------*/
/* See DML_partition_in below for a description */
/* This routine closes a record */

uint64_t DML_partition_close_in(DML_RecordReader *dml_record_in)
{
  uint64_t nbytes = dml_record_in->nbytes;

  if(dml_record_in == NULL)return 0;
  if(dml_record_in->coords != NULL)
    free(dml_record_in->coords);
  if(dml_record_in->inbuf != NULL)
    free(dml_record_in->inbuf);
  free(dml_record_in);

  /* return the number of bytes read by this node only */
  return nbytes;
}

/*------------------------------------------------------------------*/
/* Each I/O node (or the master I/O node) reads data from its file and
   distributes it to its nodes.
   Returns the number of bytes read by this node only */

/* In order to be nondeadlocking, this algorithm requires that the set
   of nodes containing sites belonging to any single I/O node are
   disjoint from the corresponding set for any other I/O node.  This
   algorithm is intended for SINGLEFILE/SERIAL, MULTIFILE, and
   PARTFILE modes. */
uint64_t
DML_partition_in(LRL_RecordReader *lrl_record_in, 
		 void (*put)(char *buf, size_t index, int count, void *arg),
		 int count, size_t size, int word_size, void *arg, 
		 DML_Layout *layout, DML_SiteList *sites, int volfmt,
		 int serpar, DML_Checksum *checksum)
{
  double dtall=0, dtread=0, dtsend=0;
  char *buf, *inbuf;
  int my_io_node;
  int *coords;
  int this_node = layout->this_node;
  int latdim = layout->latdim;
  int *latsize = layout->latsize;
  size_t nbytes=0, max_buf_sites=1;

  timestart(dtall);

  /* Get my I/O node */
  my_io_node = DML_my_ionode(volfmt, serpar, layout);

  /* Allocate buffer for reading or receiving data */
  /* I/O node needs a large buffer.  Others only enough for one site */
  if(this_node == my_io_node) max_buf_sites = DML_max_buf_sites(size,1);
  if(max_buf_sites<1) max_buf_sites = 1;

  inbuf = DML_allocate_buf(size, &max_buf_sites);
  if(!inbuf){
    printf("%s(%d) can't malloc inbuf\n",__func__,this_node);
    return 0;
  }

  /* Allocate coordinate counter */
  coords = DML_allocate_coords(latdim, __func__, this_node);
  if(!coords) { free(inbuf); return 0; }

  /* Initialize checksum */
  DML_checksum_init(checksum);

#ifdef DML_DEBUG
  if (! DML_big_endian())
    printf("%s(%d): byte reversing %d\n",__func__,this_node,word_size);
#endif

  /* Loop over the receiving sites */
  DML_SiteRank rcv_coords;
  if(DML_init_subset_site_loop(&rcv_coords, sites) == 0) {
    free(inbuf); free(coords);
    return 0;
  }

  DML_SiteRank rcoords[max_buf_sites];
  DML_SiteRank firstrank=0, nextrank=0;
  int dest_node[max_buf_sites];
  int node_index[max_buf_sites];
  int notdone = 1;
  while(notdone) {
    int k = 0;
    do { // get list of file contiguous sites
      /* The subset_rank locates the datum for rcv_coords in the
	 record our I/O partition is reading */
      DML_SiteRank subset_rank = nextrank;
      if(serpar == DML_PARALLEL)
	subset_rank = (DML_SiteRank) DML_subset_rank(rcv_coords, sites);
      if(k==0) firstrank = subset_rank;
      else if(subset_rank!=firstrank+k) break;
      /* Convert lexicographic rank to coordinates */
      DML_lex_coords(coords, latdim, latsize, rcv_coords);
      rcoords[k] = rcv_coords;
      /* The node that gets the next datum */
      dest_node[k] = layout->node_number(coords);
      node_index[k] = layout->node_index(coords);
      k++;
      notdone = DML_next_subset_site(&rcv_coords, sites);
    } while(k<max_buf_sites && notdone);

    /* I/O node reads the next value */
    if(this_node == my_io_node) {
      int doseek = (nextrank != firstrank);
      timestart(dtread);
      int err = DML_read_buf(lrl_record_in, inbuf, firstrank, size, k, doseek);
      timestop(dtread);
      nbytes += k*size;

      if(err < 0) {
        printf("%s(%d) DML_read_buf returns error\n", __func__, this_node);
        free(inbuf); free(coords);
        return 0;
      }
    }
    nextrank = firstrank + k;

    for(int i=0; i<k; i++) {
      buf = inbuf + i*size;
      /* Send result to destination node. Avoid I/O node sending to itself. */
      if (dest_node[i] != my_io_node) {
	timestart(dtsend);
	DML_route_bytes(buf, size, my_io_node, dest_node[i]);
	timestop(dtsend);
      }
      /* Process data before inserting */
      if(this_node == dest_node[i]) {
	/* Accumulate checksum */
	DML_checksum_accum(checksum, rcoords[i], buf, size);
	/* Do byte reversal if necessary */
	if (! DML_big_endian()) DML_byterevn(buf, size, word_size);
	/* Store the data */
	put(buf, node_index[i], count, arg);
      }
    }
  }
  free(inbuf); free(coords);

  timestop(dtall);
  if(this_node==0) printf("%s times: read %.2f  send %.2f  total %.2f\n", __func__, dtread, dtsend, dtall);
  /* return the number of bytes read by this node only */
  return nbytes;
}

/*------------------------------------------------------------------*/
/* The master node reads all the global data at once,
   broadcasts to all nodes and calls "put" */

/* Returns the number of bytes read */

size_t DML_global_in(LRL_RecordReader *lrl_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  int count, size_t size, int word_size, void *arg, 
          DML_Layout* layout, int volfmt, int broadcast_globaldata,
	  DML_Checksum *checksum)
{
  char *buf;
  int this_node = layout->this_node;
  size_t nbytes = 0;
  char myname[] = "DML_global_in";

  /* Allocate buffer for datum */
  buf = (char *)malloc(size);
  if(!buf){
    printf("%s(%d) can't malloc buf\n",myname,this_node);
    return 0;
  }

  /* Initialize checksum */
  DML_checksum_init(checksum);

  if(this_node == layout->master_io_node){
    /* Read all the data */
    nbytes = LRL_read_bytes(lrl_record_in, (char *)buf, size);
    if(nbytes != size){
      free(buf); return 0;
    }
    
    /* Do checksum.  Straight crc32. */
    DML_checksum_accum(checksum, 0, buf, size);

    /* Do byte reordering if needed */
    if (! DML_big_endian())
      DML_byterevn(buf, size, word_size);
    
  }

  /* We turn off broadcasting, for example, if we are doing
     single-processor file conversion */
  if(broadcast_globaldata){
    /* Broadcast the result to node bufs */
    DML_broadcast_bytes(buf, size, this_node, layout->master_io_node);
    /* All nodes store their data. Unused site index is 0. */
    put(buf,0,count,arg);
  }
  else{
    /* Only the master I/O node stores its data */
    if(this_node == layout->master_io_node)
      put(buf,0,count,arg);
  }

  free(buf);
  return nbytes;
}

