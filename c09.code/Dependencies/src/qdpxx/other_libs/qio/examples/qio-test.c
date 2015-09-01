/* Test QIO with specified volume format and real (float) fields */
/* Writes a field of values for an 8x4^3 lattice, then a global array,
   closes the file, then reads the field data, peeks at the record
   info, then reads it and compares, then reads the global array and
   compares. */

/* C. DeTar */
/* October 18, 2004 */
/* June 12, 2005 modified */

#include <qio.h>
#include <stdio.h>
#define MAIN
#include "qio-test.h"

#define NREAL 2
#define NARRAY 3
#define NMATRIX 4 

/* Open a QIO file for writing */

QIO_Writer *open_test_output(char *filename, int volfmt, 
			     int serpar, int ildgstyle, char *myname){
  QIO_String *xml_file_out;
  char xml_write_file[] = "Dummy user file XML";
  QIO_Writer *outfile;
  QIO_Oflag oflag;

  oflag.serpar = serpar;
  oflag.ildgstyle = ildgstyle;
  oflag.ildgLFN = QIO_string_create();
  QIO_string_set(oflag.ildgLFN,"TestLFN");
  oflag.mode = QIO_TRUNC;

  QIO_verbose(QIO_VERB_DEBUG);

  /* Create the file XML */
  xml_file_out = QIO_string_create();
  QIO_string_set(xml_file_out,xml_write_file);

  /* Open the file for writing */
  outfile = QIO_open_write(xml_file_out, filename, volfmt,  
			   &layout, NULL, &oflag);
  if(outfile == NULL){
    printf("%s(%d): QIO_open_write returned NULL\n",myname,this_node);
    fflush(stdout);
    return NULL;
  }
  QIO_string_destroy(xml_file_out);
  return outfile;
}

int write_real_field(QIO_Writer *outfile, int count, 
		     float *field_out[], char *myname){
  QIO_String *xml_record_out;
  char xml_write_field[] = "Dummy user record XML for real field";
  int status;
  QIO_RecordInfo *rec_info;

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, 
				    "QDP_F_Real", "F", 0,
				    0, sizeof(float), count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_field);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_R, 
		     count*sizeof(float), sizeof(float), field_out);
  printf("%s(%d): QIO_write returns status %d\n",myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_real_field_subset(QIO_Writer *outfile, int count, 
			    float *field_out[], int lower[],
			    int upper[], int dim, char *myname){
  QIO_String *xml_record_out;
  char xml_write_field[] = "Dummy user record XML for subset of real field";
  int status;
  QIO_RecordInfo *rec_info;

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_HYPER, lower, upper, dim, 
				    "QDP_F_Real", "F", 0,
				    0, sizeof(float), count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_field);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_R, 
		     count*sizeof(float), sizeof(float), field_out);
  printf("%s(%d): QIO_write returns status %d\n",myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_su3_field(QIO_Writer *outfile, int count, 
		     suN_matrix *field_out[], char *myname){
  QIO_String *xml_record_out;
  char xml_write_field[] = "Dummy user record XML for su3 field";
  int status;
  QIO_RecordInfo *rec_info;

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				    "QDP_F3_ColorMatrix", "F", 3,
				    0, sizeof(suN_matrix), count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_field);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_M, 
		     count*sizeof(suN_matrix), sizeof(float), field_out);
  printf("%s(%d): QIO_write returns status %d\n",myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_real_global(QIO_Writer *outfile, int count, float array_out[], 
		      char *myname){
  QIO_String *xml_record_out;
  char xml_write_global[] = "Dummy user record XML for global";
  int status;
  QIO_RecordInfo *rec_info;

  /* Create the record info for the global array */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_global);
  rec_info = QIO_create_record_info(QIO_GLOBAL, NULL, NULL, 0,
				    "QLA_F_Real", "F", 0,
				    0, sizeof(float), count);
  /* Write the array to a file */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_r, 
		     count*sizeof(float), sizeof(float), array_out);
  printf("%s(%d): QIO_write returns status %d\n",myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;
  
  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

QIO_Reader *open_test_input(char *filename, int volfmt, int serpar,
			    char *myname){
  QIO_String *xml_file_in;
  QIO_Reader *infile;
  QIO_Iflag iflag;

  iflag.serpar = serpar;
  iflag.volfmt = volfmt;

  /* Create the file XML */
  xml_file_in = QIO_string_create();

  /* Open the file for reading */
  infile = QIO_open_read(xml_file_in, filename, &layout, NULL, &iflag);
  if(infile == NULL){
    printf("%s(%d): QIO_open_read returns NULL.\n",myname,this_node);
    return NULL;
  }

  printf("%s(%d): QIO_open_read done.\n",myname,this_node);
  printf("%s(%d): User file info is \"%s\"\n",myname,this_node,
	 QIO_string_ptr(xml_file_in));

  QIO_string_destroy(xml_file_in);
  return infile;
}

int peek_record_info(QIO_Reader *infile, char *myname){
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;

  printf("%s(%d) peeking\n",myname,this_node);
  /* Create the record XML */
  xml_record_in = QIO_string_create();
  /* Create placeholder for record_info */

  status = QIO_read_record_info(infile, &rec_info, xml_record_in);
  if(status != QIO_SUCCESS)return 1;
  
  printf("%s(%d): User record info is \"%s\"\n",myname,this_node,
	 QIO_string_ptr(xml_record_in));

  QIO_string_destroy(xml_record_in);
  return 0;
}


int read_real_field(QIO_Reader *infile, int count, 
		    float *field_in[], char *myname)
{
  int status;
  
  /* Read the field record */
  status = QIO_read_record_data(infile, vput_R, sizeof(float)*count, 
				sizeof(float), field_in);
  printf("%s(%d): QIO_read_record_data returns status %d\n",
	 myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;
  return 0;
}
 
int read_real_field_subset(QIO_Reader *infile, int count, 
			   float *field_in[], char *myname)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  
  /* Create the record XML */
  xml_record_in = QIO_string_create();

  /* Read the field record */
  status = QIO_read(infile, &rec_info, xml_record_in,
		    vput_R, sizeof(float)*count, 
		    sizeof(float), field_in);
  printf("%s(%d): QIO_read_record_data returns status %d\n",
	 myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;
  return 0;
}
 
int read_su3_field(QIO_Reader *infile, int count, 
		    suN_matrix *field_in[], char *myname)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  
  /* Create the record XML */
  xml_record_in = QIO_string_create();

  /* Read the field record */
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_M, sizeof(suN_matrix)*count, sizeof(float), field_in);
  printf("%s(%d): QIO_read_record_data returns status %d\n",
	 myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;
  return 0;
}

int read_real_global(QIO_Reader *infile, int count, 
		    float array_in[], char *myname)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;

  /* Create the record XML */
  xml_record_in = QIO_string_create();

  /* Read the global array record */
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_r, sizeof(float)*count, sizeof(float), array_in);
  printf("%s(%d): QIO_read returns status %d\n",
	 myname,this_node,status);
  if(status != QIO_SUCCESS)return 1;
  printf("%s(%d): User record info is \"%s\"\n",myname,this_node,
	 QIO_string_ptr(xml_record_in));

  QIO_string_destroy(xml_record_in);
  return 0;
}


int qio_test(int output_volfmt, int output_serpar, int ildgstyle, 
	     int input_volfmt, int input_serpar, int argc, char *argv[]){

  float array_in[NARRAY], array_out[NARRAY];
  float *field_in[NREAL], *subset_in[NREAL], 
    *field_out[NREAL], *subset_out[NREAL];
  suN_matrix *field_su3_out[NMATRIX], *field_su3_in[NMATRIX];
  QIO_Writer *outfile;
  QIO_Reader *infile;
  float diff_field = 0, diff_array = 0, diff_su3 = 0, diff_subset = 0;
  QMP_thread_level_t provided;
  int status;
  int sites_on_node = 0;
  int i,volume;
  char filename[] = "binary_test";
  int dim = 4;
  int lower[4] = {1, 0, 0, 2};
  int upper[4] = {2, 3, 3, 2};
  char myname[] = "qio_test";
  
  /* Start message passing */
  QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &provided);

  this_node = mynode();
  printf("%s(%d) QMP_init_msg_passing done\n",myname,this_node);

  /* Lattice dimensions */
  lattice_dim = 4;
  lattice_size[0] = 8;
  lattice_size[1] = 4;
  lattice_size[2] = 4;
  lattice_size[3] = 4;

  volume = 1;
  for(i = 0; i < lattice_dim; i++){
    volume *= lattice_size[i];
  }

  /* Set the mapping of coordinates to nodes */
  if(setup_layout(lattice_size, 4, QMP_get_number_of_nodes())!=0)
    return 1;
  printf("%s(%d) layout set for %d nodes\n",myname,this_node,
	 QMP_get_number_of_nodes());
  sites_on_node = num_sites(this_node);

  /* Build the layout structure */
  layout.node_number     = node_number;
  layout.node_index      = node_index;
  layout.get_coords      = get_coords;
  layout.num_sites       = num_sites;
  layout.latsize         = lattice_size;
  layout.latdim          = lattice_dim;
  layout.volume          = volume;
  layout.sites_on_node   = sites_on_node;
  layout.this_node       = this_node;
  layout.number_of_nodes = QMP_get_number_of_nodes();

  /* Open the test output file */
  outfile = open_test_output(filename, output_volfmt, output_serpar, 
			     ildgstyle, myname);
  if(outfile == NULL)return 1;

  /* If this is not the ILDG file test */
  if(ildgstyle == QIO_ILDGNO){
    /* Create the test output field */
    status = vcreate_R(field_out, NREAL);
    if(status)return status;
    
    /* Set some values for the field */
    vset_R(field_out, NREAL);
    
    /* Write the real test field */
    status = write_real_field(outfile, NREAL, field_out, myname);
    if(status)return status;
    
    /* Write a subset of the real test field */
    status = write_real_field_subset(outfile, NREAL, field_out, 
				     lower, upper, dim, myname);
    if(status)return status;
    
    /* Set some values for the global array */
    for(i = 0; i < NARRAY; i++)
      array_out[i] = i;
    
    /* Write the real global array */
    status = write_real_global(outfile, NARRAY, array_out, myname);
    if(status)return status;
  }

  /* Create the test output su3 field */
  status = vcreate_M(field_su3_out, NMATRIX);
  if(status)return status;

  /* Set some values for the su3 field */
  vset_M(field_su3_out, NMATRIX);

  /* Write the su3 test field */
  status = write_su3_field(outfile, NMATRIX, field_su3_out, myname);
  if(status)return status;

  /* Close the file */
  QIO_close_write(outfile);
  printf("%s(%d): Closed file for writing\n",myname,this_node);

  /* Set up a dummy input field */
  status = vcreate_R(field_in, NREAL);
  if(status)return status;
    
  /* Set up a dummy input field for subset */
  status = vcreate_R(subset_in, NREAL);
  if(status)return status;
    
  /* Set up a dummy input SU(N) field */
  status = vcreate_M(field_su3_in, NMATRIX);
  if(status)return status;

  /* Open the test file for reading */
  infile = open_test_input(filename, input_volfmt, input_serpar, myname);
  if(infile == NULL)return 1;

  if(ildgstyle == QIO_ILDGNO){
    /* Peek at the field record */
    status = peek_record_info(infile, myname);
    if(status != QIO_SUCCESS)return status;
    /* Skip the record */

#if(0)
    
    /* Skip the field */
    status = QIO_next_record(infile);
    if(status != QIO_SUCCESS)return status;
    
#else
    
    /* Read the field record */
    printf("%s(%d) reading real field\n",myname,this_node); fflush(stdout);
    status = read_real_field(infile, NREAL, field_in, myname);
    if(status)return status;
    
#endif

    /* Read the subset of the field */
    printf("%s(%d) reading subset of real field\n",
	   myname,this_node); fflush(stdout);
    status = read_real_field_subset(infile, NREAL, subset_in, myname);
    if(status)return status;
    
    /* Read the global array record */
    printf("%s(%d) reading global field\n",myname,this_node); fflush(stdout);
    status = read_real_global(infile, NARRAY, array_in, myname);
    if(status)return status;

  }    

  /* Read the su3 field record */
  printf("%s(%d) reading su3 field\n",myname,this_node); fflush(stdout);
  status = read_su3_field(infile, NMATRIX, field_su3_in, myname);
  if(status)return status;

  /* Close the file */
  QIO_close_read(infile);
  printf("%s(%d): Closed file for reading\n",myname,this_node);

  if(ildgstyle == QIO_ILDGNO){

    /* Compare the input and output fields */
    diff_field = vcompare_R(field_out, field_in, NREAL);
    if(this_node == 0){
      printf("%s(%d): Comparison of in and out real fields |in - out|^2 = %e\n",
	     myname,this_node,diff_field);
    }
    
    /* Create the subset output field */
    status = vcreate_R(subset_out, NREAL);
    if(status)return status;

    /* Copy the subset */
    vsubset_R(subset_out, field_out, lower, upper, NREAL);
    
    /* Compare the input and output subsets */
    diff_subset = vcompare_R(subset_out, subset_in, NREAL);
    if(this_node == 0){
      printf("%s(%d): Comparison of subsets of in and out real fields |in - out|^2 = %e\n",
	     myname,this_node,diff_subset);
    }
    
    /* Compare the input and output global arrays */
    diff_array = vcompare_r(array_out, array_in, NREAL);
    if(this_node == 0){
      printf("%s(%d): Comparison of in and out real global arrays |in - out|^2 = %e\n",
	     myname, this_node, diff_array);
    }
  }

  /* Compare the input and output suN fields */
  diff_su3 = vcompare_M(field_su3_out, field_su3_in, NMATRIX);
  if(this_node == 0){
    printf("%s(%d): Comparison of in and out suN fields |in - out|^2 = %e\n",
	   myname, this_node, diff_field);
  }

  /* Clean up */
  if(ildgstyle == QIO_ILDGNO){
    vdestroy_R(field_out, NREAL);
    vdestroy_R(field_in, NREAL);
    vdestroy_R(subset_in, NREAL);
    vdestroy_R(subset_out, NREAL);
  }
  vdestroy_M(field_su3_in, NMATRIX);
  vdestroy_M(field_su3_out, NMATRIX);

  /* Shut down QMP */
  QMP_finalize_msg_passing();

  /* Report result */
  if(diff_field + diff_subset + diff_su3 + diff_array > 0){
    printf("%s(%d): Test failed\n",myname,this_node);
    return 1;
  }
  printf("%s(%d): Test passed\n",myname,this_node);

  return 0;
}
