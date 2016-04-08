/* Utility for converting gauge configuration files from NERSC 
   to USQCD SciDAC (or optionally ILDG) format */

/* This is single processor code */

/* Usage ...

   qio-convert-nersc [--ildg] nersc_file scidac_file
   [LFN string]

   The LFN string is used as the logical file name for ILDG usage.

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <qio.h>

#define LATDIM 4

static int nx, ny, nz, nt;
static size_t sites_on_node;
static size_t volume;
static int this_node;
static int number_of_nodes;

/*----------------------------------------------------------------------*/
/* NERSC reader definitions */
/*----------------------------------------------------------------------*/
#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512
#define GAUGE_VERSION_NUMBER_ARCHIVE 0x42454749  /* 1111836489 decimal */
#define ARCHIVE_3x2   0
#define ARCHIVE_3x3   1

/* Gauge Connection Header */

typedef struct {  /* Structure to hold archive header tokens */
  int ntoken;
  char **token;
  char **value;
} QCDheader ;

/* Generic header */

typedef struct {
  int            dims[LATDIM];  /* Full lattice dimensions */
  QCDheader      *hdr;          /* Archive header */
  int            dataformat;    /* 3x3 or 3x2 matrix (archive files only) */
  int            precision;     /* precision 1 or 2 (archive files for now) */
  uint32_t       nersc_chksum;  /* Checksum as announced in header */
} gauge_header;

/* Gauge file pointer */

typedef struct {
  FILE *         fp;            /* File pointer */
  gauge_header*  gh;            /* Pointer to header for file */
  char *         filename;      /* Pointer to file name string */
  int            byterevflag;   /* Byte reverse flag - used only for reading */
} gauge_file;

/* NERSC reader state */

typedef struct {
  gauge_file *gf;
  int siterank;
  double * dbuf;
  float * fbuf;
  uint32_t crc;
  uint32_t nersc_chksum;        /* Checksum computed as the file is read */
} r_serial_site_reader;

/*--------------------------------------------------------------------*/
/* Do byte reversal on n 32-bit words */

void byterevn(uint32_t w[], int n)
{
  register uint32_t old,newv;
  int j;

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
} /* byterevn */

/*--------------------------------------------------------------------*/
/* Do byte reversal on n contiguous 64-bit words */

void byterevn64(uint32_t w[], int n)
{
  uint32_t tmp;
  int j;

  assert(sizeof(uint32_t) == 4);
  
  /* First swap pairs of 32-bit words */
  for(j=0; j<n; j++){
    tmp = w[2*j];
    w[2*j] = w[2*j+1];
    w[2*j+1] = tmp;
  }

  /* Then swap bytes in 32-bit words */
  byterevn(w, 2*n);
}

/*--------------------------------------------------------------------*/
QCDheader * qcdhdr_get_hdr(FILE *in)
{
  char line[MAX_LINE_LENGTH];
  int n,len;
  QCDheader *hdr;
  char **tokens, **values;
  char *p, *q;

  /* Begin reading, and check for "BEGIN_HEADER" token */
  fgets(line,MAX_LINE_LENGTH,in);
  /*
  if (strcmp(line,"BEGIN_HEADER\n")!=0)
    error_exit("qcdhdr_get_hdr: Missing \"BEGIN_HEADER\"; punting \n");
  */
  /* Allocate space for QCDheader and its pointers */
  tokens = (char **) malloc(MAX_TOKENS*sizeof(char *));
  values = (char **) malloc(MAX_TOKENS*sizeof(char *));
  hdr = (QCDheader *) malloc(sizeof(QCDheader));
  if(tokens == NULL || values == NULL || hdr == NULL){
    fprintf(stderr,"Can't malloc space for archive header\n");
    return NULL;
  }
  hdr->token = tokens;
  hdr->value = values;

  /* Begin loop on tokens */
  n = 0;
  printf("Archive header:\n");
  while (1) {
    fgets(line,MAX_LINE_LENGTH,in);
    printf("%s", line);

    if (strcmp(line,"END_HEADER\n")==0) break;

    /* Tokens are terminated by a space */
    q = strchr(line, (int)' ');

    /* Overwrite space with a terminating null */
    *q = '\0';
    len = strlen(line);

    /* allocate space and copy the token in to it */
    p = (char *)malloc(len+1);
    hdr->token[n] = p;
    strcpy(p,line);

    q = strchr(++q, (int)'='); q++;
    len = strlen(q);
    q[len-1] = 0;
    p = (char *)malloc(len);
    hdr->value[n] = p;
    strcpy(p,q);
    n++;
  }
  hdr->ntoken = n;
  return hdr;
}

/*--------------------------------------------------------------------*/
/* Destroy header */
void qcdhdr_destroy_hdr(QCDheader *hdr){
  int i;
  
  if(hdr == NULL)return;

  for(i = 0; i < hdr->ntoken; i++){
    free(hdr->value[i]);
    free(hdr->token[i]);
  }

  free(hdr->token);
  free(hdr->value);
  free(hdr);
}

/*--------------------------------------------------------------------*/
gauge_file *create_gauge_file(char *filename){

  gauge_file *gf;

  /* Allocate gauge file structure */
  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL){
    fprintf(stderr, "Can't malloc file structure.\n");
    return NULL;
  }

  /* Allocate gauge header structure */
  gf->gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gf->gh == NULL){
    fprintf(stderr, "Can't malloc header structure.\n");
    return NULL;
  }

  /* File name */
  gf->filename = (char *)malloc(strlen(filename)+1);
  if(gf->filename == NULL){
    fprintf(stderr, "Can't malloc filename.\n");
    return NULL;
  }
  strcpy(gf->filename,filename);

  return gf;
}

/*--------------------------------------------------------------------*/
void destroy_gauge_file(gauge_file *gf){
  if(gf == NULL)return;
  if(gf->gh != NULL){
    if(gf->gh->hdr != NULL)qcdhdr_destroy_hdr(gf->gh->hdr);
    free(gf->gh);
  }
  free(gf);
}

/*--------------------------------------------------------------------*/
gauge_file *open_file_read_magic_no(char *filename){

  gauge_file *gf;
  uint32_t magic_number, tmp;

  /* Allocate gauge file pointer */
  gf = create_gauge_file(filename);
  if(gf == NULL){
    return NULL;
  }

  /* Open the file */
  gf->fp = fopen(filename, "rb");
  if(gf->fp == NULL){
    fprintf(stderr, "Can't open %s\n", filename);
    return NULL;
  }

  /* Read the magic number */
  if(fread(&magic_number, sizeof(magic_number), 1, gf->fp) != 1){
    fprintf(stderr, "Can't read magic number in %s.\n",filename);
    return NULL;
  }

  /* Verify magic number.  Use it to determine whether byte reversal
     will be required. */
  if(magic_number != GAUGE_VERSION_NUMBER_ARCHIVE){
    tmp = magic_number;
    byterevn(&tmp, 1);
    if(tmp != GAUGE_VERSION_NUMBER_ARCHIVE){
      fprintf(stderr, "%s does not look like a NERSC archive file.\n",
	      filename);
      fprintf(stderr, "Quitting.\n");
      return NULL;
    }
    gf->byterevflag = 1;
  } else {
    gf->byterevflag = 0;
  }
  return gf;
}

/*--------------------------------------------------------------------*/
/* NERSC reader utilities */
/*--------------------------------------------------------------------*/
int qcdhdr_get_str(char *s, QCDheader *hdr, char **q) {     
  /* find a token and return the value */
  int i;
  for (i=0; i<(char)(*hdr).ntoken; i++) {
    if (strcmp(s,(char *)(*hdr).token[i])==0) {
      *q = (*hdr).value[i];
      return SUCCESS;
    }
  }
  *q = NULL;
  return FAILURE;
}
  
int qcdhdr_get_int(char *s,QCDheader *hdr,int *q) {
  char *p;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return FAILURE;
  sscanf(p,"%d",q);
  return SUCCESS;
}
int qcdhdr_get_int32x(char *s,QCDheader *hdr,uint32_t *q) {
  char *p;
  int r;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return FAILURE;
  sscanf(p,"%x",&r);
  *q = r;
  return SUCCESS;
}

void error_exit(char *s) { fprintf(stderr,"%s\n",s); exit(1);}

void complete_U(float *u) {
  u[12] = u[ 2]*u[10] - u[ 4]*u[ 8] - u[ 3]*u[11] + u[ 5]*u[ 9];
  u[13] = u[ 4]*u[ 9] - u[ 2]*u[11] + u[ 5]*u[ 8] - u[ 3]*u[10];
  u[14] = u[ 4]*u[ 6] - u[ 0]*u[10] - u[ 5]*u[ 7] + u[ 1]*u[11];
  u[15] = u[ 0]*u[11] - u[ 4]*u[ 7] + u[ 1]*u[10] - u[ 5]*u[ 6];
  u[16] = u[ 0]*u[ 8] - u[ 2]*u[ 6] - u[ 1]*u[ 9] + u[ 3]*u[ 7];
  u[17] = u[ 2]*u[ 7] - u[ 0]*u[ 9] + u[ 3]*u[ 6] - u[ 1]*u[ 8];
}


void complete_Ud(double *u) {
  u[12] = u[ 2]*u[10] - u[ 4]*u[ 8] - u[ 3]*u[11] + u[ 5]*u[ 9];
  u[13] = u[ 4]*u[ 9] - u[ 2]*u[11] + u[ 5]*u[ 8] - u[ 3]*u[10];
  u[14] = u[ 4]*u[ 6] - u[ 0]*u[10] - u[ 5]*u[ 7] + u[ 1]*u[11];
  u[15] = u[ 0]*u[11] - u[ 4]*u[ 7] + u[ 1]*u[10] - u[ 5]*u[ 6];
  u[16] = u[ 0]*u[ 8] - u[ 2]*u[ 6] - u[ 1]*u[ 9] + u[ 3]*u[ 7];
  u[17] = u[ 2]*u[ 7] - u[ 0]*u[ 9] + u[ 3]*u[ 6] - u[ 1]*u[ 8];
}

/*--------------------------------------------------------------------*/
/* Read the entire archive file header */

gauge_file *read_archive_lat_hdr(char *filename){
  FILE* fp;
  gauge_file *gf;
  gauge_header *gh;
  QCDheader *hdr;
  char *datatype;
  char *floatpt;

  gf = open_file_read_magic_no(filename);
  if(gf == NULL){
    return NULL;
  }

  fp = gf->fp;
  gh = gf->gh;

  /* Read the entire header of the archive file */
  hdr = qcdhdr_get_hdr(fp);
  gh->hdr = hdr;
  
  /* Get dimensions */
  if (qcdhdr_get_int("DIMENSION_1",hdr,gh->dims+0)==FAILURE)
    error_exit("DIMENSION_1 not present");
  if (qcdhdr_get_int("DIMENSION_2",hdr,gh->dims+1)==FAILURE)
    error_exit("DIMENSION_2 not present");
  if (qcdhdr_get_int("DIMENSION_3",hdr,gh->dims+2)==FAILURE)
    error_exit("DIMENSION_3 not present");
  if (qcdhdr_get_int("DIMENSION_4",hdr,gh->dims+3)==FAILURE)
    error_exit("DIMENSION_4 not present");
  
  /* Get archive checksum */
  if (qcdhdr_get_int32x("CHECKSUM",hdr,&gh->nersc_chksum)==FAILURE)
    error_exit("CHECKSUM not present");
  
  /* Get archive datatype */
  if (qcdhdr_get_str("DATATYPE",hdr,&datatype)==FAILURE)
    error_exit("DATATYPE not present");
  /* Two choices currently */
  gh->dataformat = ARCHIVE_3x2;
  if(strcmp(" 4D_SU3_GAUGE_3x3",datatype) == 0)
    gh->dataformat = ARCHIVE_3x3;
  
  /* Get archive floating point format */
  gh->precision = 1;
  if (qcdhdr_get_str("FLOATING_POINT",hdr,&floatpt)==FAILURE)
    fprintf(stderr,"FLOATING_POINT tag not present.  Assuming IEEE32BIG.\n");
  else if(strcmp(" IEEE64BIG",floatpt) == 0)
    gh->precision = 2;
  
  return gf;
}


/*----------------------------------------------------------------------*/
/* Echo the archive header in the user record XML 
   and locate specific metadata required for the USQCD record XML */

static QIO_String *xml_record;

QIO_String *create_recxml(QCDheader *hdr){
  QIO_String *buf;
  QIO_USQCDLatticeInfo *record_info;
  char missing[] = "missing";
  char *plaqstring;
  char *linktrstring;
  int i;

  buf = QIO_string_create();
  if(buf == NULL){
    fprintf(stderr,"No room for record XML\n");
  }

  /* Get specific metadata */

  if (qcdhdr_get_str("PLAQUETTE",hdr,&plaqstring)==FAILURE){
    fprintf(stderr,"Warning: PLAQUETTE tag not present");
    plaqstring = missing;
  }

  if (qcdhdr_get_str("LINK_TRACE",hdr,&linktrstring)==FAILURE){
    fprintf(stderr,"Warning: LINK_TRACE tag not present");
    linktrstring = missing;
  }

  /* Copy tokens and values to buf */
  for(i = 0; i < hdr->ntoken; i++){
    QIO_string_append(buf, hdr->token[i]);
    QIO_string_append(buf, " =");  /* Values seem to have a leading blank */
    QIO_string_append(buf, hdr->value[i]);
    QIO_string_append(buf, "\n");
  }

  record_info = QIO_create_usqcd_lattice_info(plaqstring, linktrstring, 
					      QIO_string_ptr(buf));
  xml_record = QIO_string_create();
  QIO_encode_usqcd_lattice_info(xml_record, record_info);
  QIO_destroy_usqcd_lattice_info(record_info);
  
  return xml_record;
}

/*------------------------------------------------------------------*/
/* Layout utilities */
/*------------------------------------------------------------------*/

/* This is a scalar application, so we simply use lexicographic order */

static int squaresize[LATDIM];	   /* dimensions of hypercubes =
				      lattice dimensions */

/*------------------------------------------------------------------*/
/* Convert rank to coordinates */

static void lex_coords(int coords[], const int dim, const int size[], 
	   const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinates to rank */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Layout initialization */

void setup_layout(){

  squaresize[0] = nx; squaresize[1] = ny;
  squaresize[2] = nz; squaresize[3] = nt;

  this_node = 0;
  number_of_nodes = 1;

  /* Number of sites on node */
  sites_on_node =
    squaresize[0]*squaresize[1]*squaresize[2]*squaresize[3];
}

/*------------------------------------------------------------------*/
/* Layout utilities */
/*------------------------------------------------------------------*/
int node_number(const int coords[]) {
  return 0;
}

/*------------------------------------------------------------------*/
int node_index(const int coords[]) {
  return lex_rank(coords,4,squaresize);
}

/*------------------------------------------------------------------*/
int num_sites(int node) {
  return  sites_on_node;
}

/*------------------------------------------------------------------*/
/* Map node number and index to coordinates  */
/* (The inverse of node_number and node_index) */

void get_coords(int coords[], int node, const int index){
  assert(node == 0);
  lex_coords(coords, LATDIM, squaresize, index);
}

/*------------------------------------------------------------------*/
int io_node(int node){
  return node;
}

/*----------------------------------------------------------------------*/
/* Initialize archive reader */

void r_serial_reader_start(gauge_file *gf, r_serial_site_reader *state)
{
  /* gf  = gauge configuration file structure */
  /* state of the writer for a single site */

  char myname[] = "r_serial_reader_start";
  gauge_header *gh = gf->gh;
  int dataformat = gh->dataformat;
  double *dbuf;
  float *fbuf;
  int realspersite;

  if(dataformat == ARCHIVE_3x2)realspersite = 48;
  else realspersite = 72;

  /* Allocate read buffer according to precision.  We will read only one site
     at a time */
  fbuf = NULL; dbuf = NULL;
  if(gh->precision == 1)
    fbuf = (float *)malloc(realspersite*sizeof(float));
  else
    dbuf = (double *)malloc(realspersite*sizeof(double));
  if((fbuf == NULL) && (dbuf == NULL))
    {
      printf("%s: Can't malloc read buffer\n",myname);
      fflush(stdout);
      exit(1);
    }
  
  /* File should already be positioned at the start of the data */
  
  /* Save read state */
  state->gf              = gf;
  state->siterank        = 0;
  state->dbuf            = dbuf;
  state->fbuf            = fbuf;
  state->crc             = 0;
  state->nersc_chksum    = 0;
}
/*----------------------------------------------------------------------*/
/* Factory function called by the SciDAC reader.  Reads the data from
   the archive file in the order this function is called, regardless of
   "index".  Count should be 4 for four su3_matrices.  The
   pass-through arg is a pointer to the reader state */

void r_serial_reader(char *buf, size_t index, int count, void *arg)
{
  r_serial_site_reader *state = (r_serial_site_reader *)arg;
  gauge_file *gf    = state->gf;
  double *dbuf = state->dbuf;
  double *qd;
  double Ud[18];
  float  *fbuf = state->fbuf;
  float  *q;
  float  U[18];
  char myname[]     = "r_serial_reader";

  FILE *fp = gf->fp;
  gauge_header *gh = gf->gh;
  int precision = gh->precision;
  int dataformat  = gh->dataformat;
  int realspersite;
  int mu, p;

  if(dataformat == ARCHIVE_3x2)realspersite = 48;
  else realspersite = 72;

  if(state->siterank != index){
    printf("%s: expected index %d but got index %lu\n",
	   myname,state->siterank,index);
  }

  assert(count == 4);

  if(precision == 1){
    if( fread(fbuf,realspersite*sizeof(float),1,fp) != 1)
      {
	printf("%s: gauge configuration read error %d file %s\n",
	       myname,errno,gf->filename); 
	fflush(stdout); exit(1);
      }
    /* Byte reverse if necessary */
    if(gf->byterevflag)
      byterevn((uint32_t *)fbuf, realspersite);

    /* Compute NERSC checksum, do third row reconstruction if needed,
       do crc32 checksum and copy to buf for QIO */
    q = fbuf;
    for (mu=0;mu<4;mu++) {
      for (p=0;p<realspersite/4;p++) {
	state->nersc_chksum += *(uint32_t *) q;
	U[p] = (float) *(q++);
      }

      if(dataformat == ARCHIVE_3x2) complete_U(U);

      state->crc = 
	DML_crc32(state->crc, (char *)U, 18*sizeof(float));

      memcpy(buf + mu*18*sizeof(float), U, 18*sizeof(float));

    }

  } else { /* precision == 2 */
    if( fread(dbuf,realspersite*sizeof(double),1,fp) != 1)
      {
	printf("%s: gauge configuration read error %d file %s\n",
	       myname,errno,gf->filename); 
	fflush(stdout); exit(1);
      }
    
    if(gf->byterevflag)
      byterevn64((uint32_t *)dbuf, realspersite);
    
    qd = dbuf;
    for (mu=0;mu<4;mu++) {
      for (p=0;p<realspersite/4;p++) {
	state->nersc_chksum += *(uint32_t *) qd;
	state->nersc_chksum += *((uint32_t *) qd + 1);
	Ud[p] = (double) *(qd++);
      }
      if(dataformat == ARCHIVE_3x2) complete_Ud(Ud);
      
      state->crc = 
	DML_crc32(state->crc, (char *)Ud, 18*sizeof(double));

      memcpy(buf + mu*18*sizeof(double), Ud, 18*sizeof(double));
    }
  }

  state->siterank++;

} /* r_serial_reader */

/*----------------------------------------------------------------------*/
void r_serial_reader_close(r_serial_site_reader *state)
{
  double *dbuf = state->dbuf;
  float  *fbuf = state->fbuf;

  if(dbuf != NULL)free(dbuf);
  if(fbuf != NULL)free(fbuf);
}

/*----------------------------------------------------------------------*/
void build_qio_layout(QIO_Layout *layout){
  static int lattice_size[LATDIM];

  lattice_size[0] = nx;
  lattice_size[1] = ny;
  lattice_size[2] = nz;
  lattice_size[3] = nt;

  layout->node_number     = node_number;
  layout->node_index      = node_index;
  layout->get_coords      = get_coords;
  layout->num_sites       = num_sites;
  layout->latsize         = lattice_size;
  layout->latdim          = LATDIM;
  layout->volume          = volume;
  layout->sites_on_node   = sites_on_node;
  layout->this_node       = this_node;
  layout->number_of_nodes = number_of_nodes;
}

/*----------------------------------------------------------------------*/
void build_qio_filesystem(QIO_Filesystem *fs){
  fs->number_io_nodes = 0;
  fs->type = QIO_SINGLE_PATH;
  fs->my_io_node = io_node;   /* Partfile I/O uses io_node from layout*.c */
  fs->master_io_node = NULL;  /* Serial I/O uses default: node 0 */
  fs->io_node = NULL;
  fs->node_path = NULL;
}

/*----------------------------------------------------------------------*/
QIO_Writer *open_scidac_output(char *filename, int volfmt, 
			       int serpar, int ildgstyle, 
			       char *stringLFN, QIO_Layout *layout,
			       QIO_Filesystem *fs,
			       QIO_String *xml_write_file){
  QIO_Writer *outfile;
  QIO_Oflag oflag;

  /* Create the output flag structure */
  oflag.serpar = serpar;
  oflag.ildgstyle = ildgstyle;
  if(stringLFN != NULL){
    oflag.ildgLFN = QIO_string_create();
    QIO_string_set(oflag.ildgLFN, stringLFN);
  }
  else
    oflag.ildgLFN = NULL;
  oflag.mode = QIO_TRUNC;

  /* Open the file for writing */
#ifdef QIO_TRELEASE
  QIO_set_trelease(0,QIO_TRELEASE);
#endif
  outfile = QIO_open_write(xml_write_file, filename, volfmt, layout, 
			   fs, &oflag);
  if(outfile == NULL){
    printf("open_scidac_output(%d): QIO_open_write returned NULL\n",this_node);
    return NULL;
  }
  return outfile;
}

#define MAX_ILDGLFN 513

/*----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  gauge_file *gf;
  char *filename_nersc,*filename_scidac;
  QIO_Layout layout;
  QIO_Writer *outfile;
  QIO_Filesystem fs;
  QIO_RecordInfo *rec_info;
  int status;
  int ildgstyle;
  int count = 4;
  int word_size, datum_size;
  int length;
  int *dims;
  r_serial_site_reader state;
  QIO_String *xml_record_out;
  char ildg_lfn[MAX_ILDGLFN];
  QIO_String *filexml;
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Converted NERSC gauge configuration</title>";


  /* Process command line arguments */

  if(argc < 3 || 
     (argc < 4 && argv[1][0] == '-') || 
     (argc >= 4 && (strcmp(argv[1],"--ildg") != 0)))
    {
      fprintf(stderr,"Usage %s [--ildg] <NERSC file> <SciDAC file>\n",argv[0]);
      return 1;
    }

  if(argc < 4){
    filename_nersc   = argv[1];
    filename_scidac  = argv[2];
    ildgstyle = QIO_ILDGNO;
  }
  else{
    filename_nersc   = argv[2];
    filename_scidac = argv[3];
    ildgstyle = QIO_ILDGLAT;
    if(fgets(ildg_lfn, MAX_ILDGLFN, stdin) == NULL){
      fprintf(stderr,"Couldn't read the LFN\n");
      return 1;
    }
    else{
      /* Chop end-of-line character */
      length = strlen(ildg_lfn);
      if(ildg_lfn[length-1] == '\n'){
	ildg_lfn[length-1] = '\0';
      }
    }
  }

  /* Announcement */
  printf("Converting NERSC file %s to SciDAC file %s\n",
	 filename_nersc, filename_scidac);
  if(ildgstyle == QIO_ILDGLAT)
    printf("in ILDG compatible format with LFN\n%s\n",ildg_lfn);

  /* Open the NERSC file and read its header */
  gf = read_archive_lat_hdr(filename_nersc);
  if(gf == NULL)return 1;

  /* Create the user record XML */
  xml_record_out = create_recxml(gf->gh->hdr);
  if(xml_record_out == NULL)return 1;

  dims = gf->gh->dims;
  nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
  volume = nx*ny*nz*nt;

  /* Define the layout */
  setup_layout();

  /* Build the QIO layout structure */
  build_qio_layout(&layout);

  /* Open the SciDAC file for writing */
  build_qio_filesystem(&fs);
  filexml = QIO_string_create();
  QIO_string_set(filexml, default_file_xml);
  outfile = open_scidac_output(filename_scidac, QIO_SINGLEFILE,
			       QIO_SERIAL, ildgstyle, ildg_lfn, &layout,
			       &fs, filexml);
  if(outfile == NULL)exit(1);
  QIO_string_destroy(filexml);

  /* Initialize reading the NERSC lattice data */
  r_serial_reader_start(gf, &state);

  /* Data for private record XML depends on precision */

  if(gf->gh->precision == 1){
    word_size = sizeof(float);
    datum_size = 18*word_size;
    rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				      "USQCD_F3_ColorMatrix", "F", 
				      3, 0, datum_size, 4);
  } else {
    word_size = sizeof(double);
    datum_size = 18*word_size;
    rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				      "USQCD_D3_ColorMatrix", "D", 
				      3, 0, datum_size, 4);
  }

  /* Write the SciDAC record. The factory function "r_serial_reader"
     reads the site links from the NERSC file */

  status = QIO_write(outfile, rec_info, xml_record_out,
		     r_serial_reader, datum_size*count, word_size, 
		     (void *)&state);
  if(status != QIO_SUCCESS)exit(1);

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  /* Verify NERSC checksums */
  if(gf->gh->nersc_chksum != state.nersc_chksum)
    printf("NERSC checksum mismatch.  Wanted %x.  Got %x\n",
	    gf->gh->nersc_chksum, state.nersc_chksum);
  else
    printf("NERSC checksum %x OK\n",state.nersc_chksum);

  printf("SciDAC checksums %x %x\n",
	 QIO_get_writer_last_checksuma(outfile),
	 QIO_get_writer_last_checksumb(outfile));
  printf("ILDG crc32 checksum %lu\n",state.crc);

  /* Close the SciDAC file */
  QIO_close_write(outfile);

  /* Clean up */
  r_serial_reader_close(&state);
  destroy_gauge_file(gf);

  return 0;
}
