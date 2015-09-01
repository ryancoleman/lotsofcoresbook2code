/* Read and write private file, record, and checksum info for SciDAC
   binary file format */

#include <qio_config.h>
#include <stdio.h>
#include <string.h>
#include <qio.h>
#include <qioxml.h>
#include <qio_string.h>
#include <qio_stdint.h>
#include <sys/types.h>
#include <time.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif


int QIO_decode_record_info(QIO_RecordInfo *record_info, 
			QIO_String *record_string){
  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_RecordInfoWrapper wrapper = QIO_RECORD_INFO_WRAPPER;
  QIO_RecordInfo templ = QIO_RECORD_INFO_TEMPLATE;
  char *left_angle;

  /* Compatibility */
  QIO_RecordInfo_v1p0 record_info_v1p0 = QIO_RECORD_INFO_TEMPLATE_v1p0;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_RecordInfo));

  /* Start parsing record_string */
  /* Check leading tag, which is probably the info phrase "<?xml ...?>" */
  /* We ignore it if it is there */
  tmp_pt = QIO_next_tag(parse_pt, tag, &left_angle);
  if(strcmp(tag,QIO_QUESTXML)==0){
    /* Found ?xml, so resume parsing after the closing ">", ignoring
       the field. Otherwise, leave the parse_pt at its initial value */
    parse_pt = tmp_pt;
  }

  /* Open top-level tag (wrapper) and extract string containing tags */
  parse_pt = QIO_get_tag_value(parse_pt, tag, tags_string);
  QIO_decode_as_string (tag, tags_string, &wrapper.recordinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.recordinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_record_info_tag_string(&wrapper);

  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);
    
    QIO_decode_as_string (tag,value_string,&record_info->version);
    QIO_decode_as_string (tag,value_string,&record_info->date);
    QIO_decode_as_int    (tag,value_string,&record_info->recordtype);
    QIO_decode_as_int    (tag,value_string,&record_info->spacetime);
    QIO_decode_as_intlist(tag,value_string,&record_info->hyperlower); 
    QIO_decode_as_intlist(tag,value_string,&record_info->hyperupper); 
    QIO_decode_as_string (tag,value_string,&record_info->datatype);
    QIO_decode_as_string (tag,value_string,&record_info->precision);
    QIO_decode_as_int    (tag,value_string,&record_info->colors);
    QIO_decode_as_int    (tag,value_string,&record_info->spins);
    QIO_decode_as_int    (tag,value_string,&record_info->typesize);
    QIO_decode_as_int    (tag,value_string,&record_info->datacount);
    /* compatibility */	 
    QIO_decode_as_int    (tag,value_string,&record_info_v1p0.globaldata);
  }

  /* Backward compatibility */

  /* Convert version 1.0 record_info structure to version 1.1 */

  if(strcmp("1.0",QIO_get_record_info_version(record_info)) == 0){

    /* Version 1.1 added a new record type: hypercube subset and
       requires a list of lower and upper coordinate bounds to specify
       the hypercube.  These parameters are ignored with the
       QIO_GLOBAL and QIO_FIELD record types, the only ones used in
       version 1, so we don't need to set default values. */

    /* An earlier defective version (also labeled 1.0) was missing the
       "globaldata" parameter altogether. */
    /* If the old globaldata tag is missing, insert a default value */
    
    if(QIO_check_int_occur(&record_info_v1p0.globaldata) != 0){
      record_info_v1p0.globaldata.occur = 1;
      /* Default is "field" record type */
      record_info_v1p0.globaldata.value = QIO_FIELD;
    }

    /* Also the "globaldata" member was renamed "recordtype".  So just
       copy the old parameter value. */

    record_info->recordtype.occur = 1;
    record_info->recordtype.value = QIO_get_globaldata(&record_info_v1p0);

  }

  /* Check for completeness */
  
  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_int_occur   (&record_info->recordtype);
  errors += QIO_check_string_occur(&record_info->datatype);
  errors += QIO_check_string_occur(&record_info->precision);
  errors += QIO_check_int_occur   (&record_info->typesize);
  errors += QIO_check_int_occur   (&record_info->datacount);

  /* Requirements for hypercube record type */

  if(QIO_get_recordtype(record_info) == QIO_HYPER){
    errors += QIO_check_int_occur (&record_info->spacetime);
    errors += QIO_check_intarray_occur (&record_info->hyperlower);
    errors += QIO_check_intarray_occur (&record_info->hyperupper);
  }

  return errors;
}

void QIO_encode_record_info(QIO_String *record_string, 
			  QIO_RecordInfo *record_info){
  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  QIO_RecordInfoWrapper wrapper = QIO_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->date, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->recordtype, &remainder);
  if(QIO_get_recordtype(record_info) == QIO_HYPER){
    buf = QIO_encode_as_int    (buf,&record_info->spacetime, &remainder);
    n = record_info->spacetime.value;
    buf = QIO_encode_as_intlist(buf,&record_info->hyperlower, n, &remainder);
    buf = QIO_encode_as_intlist(buf,&record_info->hyperupper, n, &remainder);
  }
  buf = QIO_encode_as_string(buf,&record_info->datatype, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->precision, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->colors, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->spins, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->typesize, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->datacount, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_record_tag_string(&wrapper, recordinfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(record_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(record_string);
  remainder = QIO_string_length(record_string);
  
  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("QIO_encode_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.recordinfo_tags, &remainder);
  }
}

/* Decode private SciDAC file info string */

int QIO_decode_file_info(QIO_FileInfo *file_info, 
			  QIO_String *file_string){
  char *parse_pt = QIO_string_ptr(file_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_FileInfoWrapper wrapper = QIO_FILE_INFO_WRAPPER;
  QIO_FileInfo templ = QIO_FILE_INFO_TEMPLATE;
  char *left_angle;

  /* Compatibility */
  QIO_FileInfo_v1p0 file_info_v1p0  = QIO_FILE_INFO_TEMPLATE_v1p0;
  
  /* Initialize file info structure from a template */
  memcpy(file_info, &templ, sizeof(QIO_FileInfo));

  /* Start parsing file_string */
  /* Check leading tag, which is probably the info phrase "<?xml ...?>" */
  /* We ignore it if it is there */
  tmp_pt = QIO_next_tag(parse_pt, tag, &left_angle);
  if(strcmp(tag,QIO_QUESTXML)==0){
    /* Found ?xml, so resume parsing after the closing ">", ignoring
       the field. Otherwise, leave the parse_pt at its initial value. */
    parse_pt = tmp_pt;
  }

  /* Open top-level tag (wrapper) and extract string containing tags */
  parse_pt = QIO_get_tag_value(parse_pt, tag, tags_string);
  QIO_decode_as_string (tag, tags_string, &wrapper.fileinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.fileinfo_tags))return QIO_BAD_XML;
  /* Otherwise start parsing the enclosed string of tags */
  parse_pt = QIO_get_file_info_tag_string(&wrapper);

  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);
    
    QIO_decode_as_string (tag,value_string,&file_info->version);
    QIO_decode_as_int    (tag,value_string,&file_info->spacetime);
    QIO_decode_as_intlist(tag,value_string,&file_info->dims);
    QIO_decode_as_int    (tag,value_string,&file_info->volfmt);
    /* compatibility */
    QIO_decode_as_int    (tag,value_string,&file_info_v1p0.multifile);
  }

  /* Check for completeness */
  
  errors += QIO_check_string_occur  (&file_info->version);
  errors += QIO_check_int_occur     (&file_info->spacetime);
  errors += QIO_check_intarray_occur(&file_info->dims);
  errors += QIO_check_int_occur     (&file_info->volfmt);

  /* Did we get all the spacetime dimensions */
  if(file_info->spacetime.value != file_info->dims.n){
    printf("QIO_decode_file_info: mismatch in spacetime dimensions\n");
    errors++;
  }

  /* Backward compatibility */

  /* Convert version 1.0 file_info structure to version 1.1 */

  if(strcmp("1.0",QIO_get_file_version(file_info)) == 0){

    /* Version 1.0 had a multifile parameter where the volfmt paramter
       now appears.  The multifile parameter gave the file count. A 1
       implied singlefile and greater than 1 implied multifile.  There
       was no partfile format in 1.0. In version 1.1 the multifile
       flag was changed to specify the volume format: SINGLEFILE,
       MULTIFILE, PARTFILE */
    if(QIO_get_multifile(&file_info_v1p0) == 1)
      QIO_insert_volfmt(file_info,QIO_SINGLEFILE);
    else
      QIO_insert_volfmt(file_info,QIO_MULTIFILE);
  }

  return errors;
}

/* Encode private SciDAC file info string */

void QIO_encode_file_info(QIO_String *file_string, 
			  QIO_FileInfo *file_info){
  char *buf;
  int remainder,n;
  char fileinfo_tags[QIO_MAXVALUESTRING];
  QIO_FileInfoWrapper wrapper = QIO_FILE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = fileinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string (buf,&file_info->version, &remainder);
  buf = QIO_encode_as_int    (buf,&file_info->spacetime, &remainder);
  buf = QIO_encode_as_intlist(buf,&file_info->dims, 
			      file_info->spacetime.value, &remainder);
  buf = QIO_encode_as_int    (buf,&file_info->volfmt, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_file_tag_string(&wrapper, fileinfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(file_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(file_string);
  remainder = QIO_string_length(file_string);

  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("QIO_encode_file_info: file_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.fileinfo_tags, &remainder);
  }
}

int QIO_decode_checksum_info(QIO_ChecksumInfo *checksum, 
			     QIO_String *file_string){
  char *parse_pt = QIO_string_ptr(file_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_ChecksumInfoWrapper wrapper = QIO_CHECKSUM_INFO_WRAPPER;
  QIO_ChecksumInfo templ = QIO_CHECKSUM_INFO_TEMPLATE;
  char *left_angle;
  
  /* Initialize from template */
  memcpy(checksum, &templ, sizeof(QIO_ChecksumInfo));

  /* Start parsing checksum_string */
  /* Check leading tag, which is probably the info phrase "<?xml ...?>" */
  /* We ignore it if it is there */
  tmp_pt = QIO_next_tag(parse_pt, tag, &left_angle);
  if(strcmp(tag,QIO_QUESTXML)==0){
    /* Found ?xml, so resume parsing after the closing ">", ignoring
       the field. Otherwise, leave the parse_pt at its initial value. */
    parse_pt = tmp_pt;
  }

  /* Open top-level tag (wrapper) and extract string containing tags */
  parse_pt = QIO_get_tag_value(parse_pt, tag, tags_string);
  QIO_decode_as_string (tag, tags_string, &wrapper.checksuminfo_tags);

  /* If outer wrapper has bad tag exit with error status */
  if(QIO_check_string_occur(&wrapper.checksuminfo_tags))return QIO_BAD_XML;
  /* Otherwise start parsing the enclosed string of tags */
  parse_pt = QIO_get_checksum_info_tag_string(&wrapper);

  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);
    
    QIO_decode_as_string (tag,value_string,&checksum->version);
    QIO_decode_as_hex32  (tag,value_string,&checksum->suma);
    QIO_decode_as_hex32  (tag,value_string,&checksum->sumb);
  }

  /* Check for completeness */
  
  errors += QIO_check_string_occur  (&checksum->version);
  errors += QIO_check_hex32_occur   (&checksum->suma);
  errors += QIO_check_hex32_occur   (&checksum->sumb);

  return errors;
}

void QIO_encode_checksum_info(QIO_String *checksum_string, 
			      QIO_ChecksumInfo *checksum){
  char *buf;
  int remainder,n;
  char checksuminfo_tags[QIO_MAXVALUESTRING];
  QIO_ChecksumInfoWrapper wrapper = QIO_CHECKSUM_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = checksuminfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string (buf,&checksum->version, &remainder);
  buf = QIO_encode_as_hex32  (buf,&checksum->suma, &remainder);
  buf = QIO_encode_as_hex32  (buf,&checksum->sumb, &remainder);

  /* Insert inner tag string into checksum wrapper structure */
  QIO_insert_checksum_tag_string(&wrapper, checksuminfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(checksum_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(checksum_string);
  remainder = QIO_string_length(checksum_string);

  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("QIO_encode_checksum_info: checksum_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.checksuminfo_tags, &remainder);
  }
}

/* Utilities for loading file_info values */

int QIO_insert_file_tag_string(QIO_FileInfoWrapper *wrapper, 
			       char *fileinfo_tags){
  wrapper->fileinfo_tags.occur = 0;
  if(!fileinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->fileinfo_tags.value, fileinfo_tags, QIO_MAXVALUESTRING-1);
  wrapper->fileinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->fileinfo_tags.occur = 1;
  if(strlen(fileinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_file_version(QIO_FileInfo *file_info, char *version){
  file_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(file_info->version.value, version, QIO_MAXVALUESTRING-1);
  file_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_spacetime_dims(QIO_FileInfo *file_info, 
			      int spacetime, int *dims){
  int i;

  file_info->spacetime.occur = 0;
  file_info->dims.occur = 0;
  if(!spacetime)return QIO_BAD_ARG;
  if(!dims)return QIO_BAD_ARG;
  file_info->spacetime.value =  spacetime;
  if(spacetime > QIO_MAXINTARRAY){
    printf("QIO_insert_spacetime_dims: spacetime %d exceeds max %d\n",
	   spacetime, QIO_MAXINTARRAY);
    return QIO_BAD_ARG;
  }
  for(i = 0; i < spacetime; i++)
    file_info->dims.value[i] = dims[i];
  file_info->spacetime.occur = 1;
  file_info->dims.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_volfmt(QIO_FileInfo *file_info, int volfmt){
  file_info->volfmt.occur = 0;
  if(volfmt!=QIO_SINGLEFILE && 
     volfmt!=QIO_MULTIFILE &&
     volfmt!=QIO_PARTFILE){
    printf("QIO_insert_volfmt: Bad volfmt parameter %d\n",volfmt);
    return QIO_BAD_ARG;
  }
  file_info->volfmt.value = volfmt;
  file_info->volfmt.occur = 1;
  return QIO_SUCCESS;
}

/* Utilities for loading record_info values */

int QIO_insert_record_tag_string(QIO_RecordInfoWrapper *wrapper, 
				 char *recordinfo_tags){
  wrapper->recordinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->recordinfo_tags.value, recordinfo_tags, 
	  QIO_MAXVALUESTRING-1);
  wrapper->recordinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->recordinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_record_version(QIO_RecordInfo *record_info, char *version){
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_record_date(QIO_RecordInfo *record_info, char *date_string){
  int n;

  record_info->date.occur = 0;
  if(!date_string)return QIO_BAD_ARG;
  strncpy(record_info->date.value, date_string, QIO_MAXVALUESTRING-1);
  /* Edit date: replace trailing end-of-line by blank and add UTC */
  n = strlen(record_info->date.value);
  if(record_info->date.value[n-1] == '\n')record_info->date.value[n-1] = ' ';
  strncpy(record_info->date.value + n,"UTC",QIO_MAXVALUESTRING - n);
  record_info->date.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->date.occur = 1;
  if(strlen(date_string) + 3 >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_recordtype(QIO_RecordInfo *record_info, int recordtype){
  record_info->recordtype.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->recordtype.value = recordtype;
  record_info->recordtype.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_hypercube_bounds(QIO_RecordInfo *record_info, 
				int *lower, int *upper, int n){
  int i;
  record_info->spacetime.occur = 0;
  record_info->hyperlower.occur = 0;
  record_info->hyperupper.occur = 0;
  if(n > QIO_MAXINTARRAY){
    printf("QIO_insert_hypercube_bounds: dimension %d exceeds max %d\n",
	   n, QIO_MAXINTARRAY);
    return QIO_BAD_ARG;
  }
  if(lower == NULL || upper == NULL)n = 0;
  for(i = 0; i < n; i++){
    record_info->hyperlower.value[i] = lower[i];
    record_info->hyperupper.value[i] = upper[i];
  }
  record_info->spacetime.value = n;
  record_info->spacetime.occur = 1;
  record_info->hyperlower.occur = 1;
  record_info->hyperupper.occur = 1;
  return QIO_SUCCESS;
}



int QIO_insert_datatype(QIO_RecordInfo *record_info, char* datatype){
  record_info->datatype.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->datatype.value, datatype, QIO_MAXVALUESTRING-1);
  record_info->datatype.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->datatype.occur = 1;
  if(strlen(datatype) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_precision(QIO_RecordInfo *record_info, char* precision){
  record_info->precision.occur = 0;
  if(!precision)return QIO_BAD_ARG;
  strncpy(record_info->precision.value, precision, QIO_MAXVALUESTRING-1);
  record_info->precision.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->precision.occur = 1;
  if(strlen(precision) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_colors(QIO_RecordInfo *record_info, int colors){
  record_info->colors.occur = 0;
  if(!colors)return QIO_BAD_ARG;
  record_info->colors.value = colors;
  record_info->colors.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_spins(QIO_RecordInfo *record_info, int spins){
  record_info->spins.occur = 0;
  if(!spins)return QIO_BAD_ARG;
  record_info->spins.value = spins;
  record_info->spins.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_typesize(QIO_RecordInfo *record_info, int typesize){
  record_info->typesize.occur = 0;
  if(!typesize)return QIO_BAD_ARG;
  record_info->typesize.value = typesize;
  record_info->typesize.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_datacount(QIO_RecordInfo *record_info, int datacount){
  record_info->datacount.occur = 0;
  if(!datacount)return QIO_BAD_ARG;
  record_info->datacount.value = datacount;
  record_info->datacount.occur = 1;
  return QIO_SUCCESS;
}

/* Utility for loading checksum values */

int QIO_insert_checksum_tag_string(QIO_ChecksumInfoWrapper *wrapper, 
				   char *checksuminfo_tags){
  wrapper->checksuminfo_tags.occur = 0;
  if(!checksuminfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->checksuminfo_tags.value, checksuminfo_tags, 
	  QIO_MAXVALUESTRING-1);
  wrapper->checksuminfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->checksuminfo_tags.occur = 1;
  if(strlen(checksuminfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_checksum_version(QIO_ChecksumInfo *checksum_info, 
				char *version){
  checksum_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(checksum_info->version.value, version, QIO_MAXVALUESTRING-1);
  checksum_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  checksum_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_suma_sumb(QIO_ChecksumInfo *checksum_info, 
			 uint32_t suma, uint32_t sumb){
  checksum_info->suma.occur = 0;
  checksum_info->sumb.occur = 0;
  if(!suma || !sumb)return QIO_BAD_ARG;
  checksum_info->suma.value = suma;
  checksum_info->sumb.value = sumb;
  checksum_info->suma.occur = 1;
  checksum_info->sumb.occur = 1;
  return QIO_SUCCESS;
}

/* Accessors for file info */

char *QIO_get_file_info_tag_string(QIO_FileInfoWrapper *wrapper){
  return wrapper->fileinfo_tags.value;
}

char *QIO_get_file_version(QIO_FileInfo *file_info){
  return file_info->version.value;
}

int QIO_get_spacetime(QIO_FileInfo *file_info){
  return file_info->spacetime.value;
}

int *QIO_get_dims(QIO_FileInfo *file_info){
  return file_info->dims.value;
}

int QIO_get_volfmt(QIO_FileInfo *file_info){
  return file_info->volfmt.value;
}

int QIO_get_multifile(QIO_FileInfo_v1p0 *file_info){
  return file_info->multifile.value;
}

int QIO_defined_spacetime(QIO_FileInfo *file_info){
  return file_info->spacetime.occur;
}

int QIO_defined_dims(QIO_FileInfo *file_info){
  return file_info->dims.occur;
}

int QIO_defined_volfmt(QIO_FileInfo *file_info){
  return file_info->volfmt.occur;
}


/* Accessors for record info */

char *QIO_get_record_info_tag_string(QIO_RecordInfoWrapper *wrapper){
  return wrapper->recordinfo_tags.value;
}

char *QIO_get_record_info_version(QIO_RecordInfo *record_info){
  return record_info->version.value;
}

int QIO_get_recordtype(QIO_RecordInfo *record_info){
  return record_info->recordtype.value;
}

int QIO_get_hyper_spacetime(QIO_RecordInfo *record_info){
  return record_info->spacetime.value;
}

int *QIO_get_hyperupper(QIO_RecordInfo *record_info){
  return record_info->hyperupper.value;
}

int *QIO_get_hyperlower(QIO_RecordInfo *record_info){
  return record_info->hyperlower.value;
}

int QIO_get_globaldata(QIO_RecordInfo_v1p0 *record_info){
  return record_info->globaldata.value;
}

char *QIO_get_datatype(QIO_RecordInfo *record_info){
  return record_info->datatype.value;
}

char *QIO_get_precision(QIO_RecordInfo *record_info){
  return record_info->precision.value;
}

char *QIO_get_record_date(QIO_RecordInfo *record_info){
  return record_info->date.value;
}

int QIO_get_colors(QIO_RecordInfo *record_info){
  return record_info->colors.value;
}

int QIO_get_spins(QIO_RecordInfo *record_info){
  return record_info->spins.value;
}

int QIO_get_typesize(QIO_RecordInfo *record_info){
  return record_info->typesize.value;
}

int QIO_get_datacount(QIO_RecordInfo *record_info){
  return record_info->datacount.value;
}

int QIO_defined_recordtype(QIO_RecordInfo *record_info){
  return record_info->recordtype.occur;
}

void QIO_set_recordtype(QIO_RecordInfo *record_info, int recordtype){
  record_info->recordtype.value = recordtype;
  record_info->recordtype.occur = 1;
}

void QIO_set_datatype(QIO_RecordInfo *record_info, char *datatype){
  strncpy(record_info->datatype.value,datatype,QIO_MAXVALUESTRING-1);
  record_info->datatype.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->datatype.occur = 1;
}

void QIO_set_precision(QIO_RecordInfo *record_info, char *precision){
  strncpy(record_info->precision.value,precision,QIO_MAXVALUESTRING-1);
  record_info->precision.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->precision.occur = 1;
}

void QIO_set_record_date(QIO_RecordInfo *record_info, char *date){
  strncpy(record_info->date.value,date,QIO_MAXVALUESTRING-1);
  record_info->date.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->date.occur = 1;
}

void QIO_set_colors(QIO_RecordInfo *record_info, int colors){
  record_info->colors.value = colors;
  record_info->colors.occur = 1;
}

void QIO_set_spins(QIO_RecordInfo *record_info, int spins){
  record_info->spins.value = spins;
  record_info->spins.occur = 1;
}

void QIO_set_typesize(QIO_RecordInfo *record_info, int typesize){
  record_info->typesize.value = typesize;
  record_info->typesize.occur = 1;
}

void QIO_set_datacount(QIO_RecordInfo *record_info, int datacount){
  record_info->datacount.value = datacount;
  record_info->datacount.occur = 1;
}

int QIO_defined_datatype(QIO_RecordInfo *record_info){
  return record_info->datatype.occur;
}

int QIO_defined_precision(QIO_RecordInfo *record_info){
  return record_info->precision.occur;
}

int QIO_defined_colors(QIO_RecordInfo *record_info){
  return record_info->colors.occur;
}

int QIO_defined_spins(QIO_RecordInfo *record_info){
  return record_info->spins.occur;
}

int QIO_defined_typesize(QIO_RecordInfo *record_info){
  return record_info->typesize.occur;
}

int QIO_defined_datacount(QIO_RecordInfo *record_info){
  return record_info->datacount.occur;
}

/* Accessors for checksum info */

char *QIO_get_checksum_info_tag_string(QIO_ChecksumInfoWrapper *wrapper){
  return wrapper->checksuminfo_tags.value;
}

uint32_t QIO_get_suma(QIO_ChecksumInfo *checksum_info){
  return checksum_info->suma.value;
}

uint32_t QIO_get_sumb(QIO_ChecksumInfo *checksum_info){
  return checksum_info->sumb.value;
}

int QIO_defined_suma(QIO_ChecksumInfo *checksum_info){
  return checksum_info->suma.occur;
}

int QIO_defined_sumb(QIO_ChecksumInfo *checksum_info){
  return checksum_info->sumb.occur;
}

/* Utilities for creating structures from templates */

QIO_FileInfo *QIO_create_file_info(int spacetime, int *dims, int volfmt){
  QIO_FileInfo templ = QIO_FILE_INFO_TEMPLATE;
  QIO_FileInfo *file_info;
  
  file_info = (QIO_FileInfo *)malloc(sizeof(QIO_FileInfo));
  if(!file_info)return NULL;
  
  memcpy(file_info, &templ, sizeof(QIO_FileInfo));
  QIO_insert_file_version(file_info,QIO_FILEFORMATVERSION);
  QIO_insert_spacetime_dims(file_info,spacetime,dims);
  QIO_insert_volfmt(file_info,volfmt);
  return file_info;
}

void QIO_destroy_file_info(QIO_FileInfo *file_info){
  free(file_info);
}

/* Compare only fields that occur in the expected structure */
int QIO_compare_file_info(QIO_FileInfo *found, QIO_FileInfo *expect,
			  char *myname, int this_node){
  int i, n, ok;
  int *dims_expect, *dims_found;

  if(QIO_defined_spacetime(expect))
    if(!QIO_defined_spacetime(found) &&
       QIO_get_spacetime(found) != QIO_get_spacetime(expect))
      {
	printf("%s(%d):Spacetime dimension mismatch expected %d found %d \n",
	       myname, this_node,
	       QIO_get_spacetime(expect), QIO_get_spacetime(found));
	return QIO_ERR_FILE_INFO;
      }
  
  if(QIO_defined_dims(expect)){
    if(!QIO_defined_dims(found))
      {
	printf("%s(%d):Dimensions missing\n",myname,this_node);
	return QIO_ERR_FILE_INFO;
      }
    
    dims_expect = QIO_get_dims(expect);
    dims_found = QIO_get_dims(found);
    n = QIO_get_spacetime(expect);
    ok = 1;
    
    for(i = 0; i < n; i++)
      if(dims_expect[i] != dims_found[i])ok = 0;
    
    if(!ok){
      printf("%s(%d): lattice dimensions do not match\n",myname,this_node);
      printf("Expected ");
      for(i = 0; i < n; i++)printf(" %d", dims_expect[i]);
      printf("\nFound   ");
      for(i = 0; i < n; i++)printf(" %d", dims_found[i]);
      printf("\n");
      return QIO_ERR_FILE_INFO;
    }
  }
  
  if(QIO_defined_volfmt(expect))
    if(!QIO_defined_volfmt(found) &&
       QIO_get_volfmt(found) != QIO_get_volfmt(expect))
      {
	printf("%s(%d):Volfmt parameter mismatch: expected %d found %d \n",
	       myname,this_node,
	       QIO_get_volfmt(expect),QIO_get_volfmt(found));
	return QIO_ERR_FILE_INFO;
      }
  
  return QIO_SUCCESS;
}

QIO_RecordInfo *QIO_create_record_info(int recordtype, int lower[],
				       int upper[], int n,
				       char *datatype, char *precision, 
				       int colors, int spins, int typesize, 
				       int datacount){
  QIO_RecordInfo templ = QIO_RECORD_INFO_TEMPLATE;
  QIO_RecordInfo *record_info;
  time_t cu_time;
  
  record_info = (QIO_RecordInfo *)malloc(sizeof(QIO_RecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  memcpy(record_info, &templ, sizeof(QIO_RecordInfo));
  QIO_insert_record_version(record_info,QIO_RECORDFORMATVERSION);
  QIO_insert_record_date(record_info,asctime(gmtime(&cu_time)));
  QIO_insert_recordtype(record_info,recordtype);
  if(recordtype == QIO_HYPER)
    QIO_insert_hypercube_bounds(record_info, lower, upper, n);
  QIO_insert_datatype(record_info,datatype);
  QIO_insert_precision(record_info,precision);
  QIO_insert_colors(record_info,colors);
  QIO_insert_spins(record_info,spins);
  QIO_insert_typesize(record_info,typesize);
  QIO_insert_datacount(record_info,datacount);
  return record_info;
}

void QIO_destroy_record_info(QIO_RecordInfo *record_info){
  free(record_info);
}

/* Compare only fields that occur in the expected record info */
int QIO_compare_record_info(QIO_RecordInfo *found, QIO_RecordInfo *expect){
  char myname[] = "QIO_compare_record_info";

  if(QIO_defined_recordtype(expect))
    if(!QIO_defined_recordtype(found) &&
       QIO_get_recordtype(found)  != QIO_get_recordtype(expect))
      {
	printf("%s:Recordtype flag mismatch expected %d found %d \n",myname,
	       QIO_get_recordtype(expect),QIO_get_recordtype(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_datatype(expect))
    if(!QIO_defined_datatype(found) && 
       strncmp(QIO_get_datatype(found),QIO_get_datatype(expect),
	       QIO_MAXVALUESTRING))
      {
	printf("%s:Datatype mismatch expected %s found %s \n",myname,
	       QIO_get_datatype(expect),QIO_get_datatype(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_precision(expect))
    if(!QIO_defined_precision(found) &&
       strncmp(QIO_get_precision(found),QIO_get_precision(expect),
	       QIO_MAXVALUESTRING))
      {
	printf("%s:Precision mismatch expected %s found %s \n",myname,
	       QIO_get_precision(expect),QIO_get_precision(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_colors(expect))
    if(!QIO_defined_colors(found) &&
       QIO_get_colors(found) != QIO_get_colors(expect))
      {
	printf("%s:Colors mismatch expected %d found %d \n",myname,
	       QIO_get_colors(expect),QIO_get_colors(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_spins(expect))
    if(!QIO_defined_spins(found) &&
       QIO_get_spins(found)  != QIO_get_spins(expect))
      {
	printf("%s:Spins mismatch expected %d found %d \n",myname,
	       QIO_get_spins(expect),QIO_get_spins(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_typesize(expect))
    if(!QIO_defined_typesize(found) &&
       QIO_get_typesize(found) != QIO_get_typesize(expect))
      {
	printf("%s:Typesize mismatch expected %d found %d \n",myname,
	       QIO_get_typesize(expect),QIO_get_typesize(found));
	return QIO_ERR_REC_INFO;
      }

  if(QIO_defined_datacount(expect))
    if(!QIO_defined_datacount(found) &&
       QIO_get_datacount(found) != QIO_get_datacount(expect))
      {
	printf("%s:Datacount mismatch expected %d found %d \n",myname,
	       QIO_get_datacount(expect),QIO_get_datacount(found));
	return QIO_ERR_REC_INFO;
      }

  return QIO_SUCCESS;
}

QIO_ChecksumInfo *QIO_create_checksum_info(uint32_t suma, uint32_t sumb){
  QIO_ChecksumInfo templ = QIO_CHECKSUM_INFO_TEMPLATE;
  QIO_ChecksumInfo *checksum_info;
  
  checksum_info = (QIO_ChecksumInfo *)malloc(sizeof(QIO_ChecksumInfo));
  if(!checksum_info)return NULL;

  memcpy(checksum_info, &templ, sizeof(QIO_ChecksumInfo));
  QIO_insert_checksum_version(checksum_info,QIO_CHECKSUMFORMATVERSION);
  QIO_insert_suma_sumb(checksum_info,suma,sumb);
  return checksum_info;
}

void QIO_destroy_checksum_info(QIO_ChecksumInfo *checksum_info){
  free(checksum_info);
}

int QIO_compare_checksum_info(QIO_ChecksumInfo *found, 
			      QIO_ChecksumInfo *expect, 
			      char *myname, int this_node){

  if(!QIO_defined_suma(found) || !QIO_defined_sumb(found)){
    printf("%s(%d): checksum info missing\n",myname,this_node);
    return QIO_ERR_CHECKSUM_INFO;
  }

  if(QIO_get_suma(expect) != QIO_get_suma(found) || 
     QIO_get_sumb(expect) != QIO_get_sumb(found)){
    printf("%s(%d): Checksum mismatch.  Found %x %x.  Expected %x %x\n",
	   myname,this_node,QIO_get_suma(found),QIO_get_sumb(found),
	   QIO_get_suma(expect),QIO_get_sumb(expect) );
    return QIO_CHECKSUM_MISMATCH;
  }
  else
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): Checksums %x %x OK\n",myname,this_node,
	     QIO_get_suma(found),QIO_get_sumb(found));
  
  return QIO_SUCCESS;
}


