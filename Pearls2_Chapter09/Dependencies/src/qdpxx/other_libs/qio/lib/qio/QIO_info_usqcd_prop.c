/* Read and write file and record info for USQCD propagator file formats */

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


void QIO_encode_usqcd_propfile_info(QIO_String *file_string, 
				    QIO_USQCDPropFileInfo *file_info)
{
  char myname[] = "QIO_encode_usqcd_propfile_info";
  char *buf;
  int remainder,n;
  char fileinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDPropFileInfoWrapper wrapper = QIO_USQCD_PROPFILE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = fileinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&file_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->type, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_usqcd_propfile_tag_string(&wrapper, fileinfo_tags);

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
    printf("%s: file_string overflow\n",myname);
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdpropfileinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info,
				   QIO_String *file_string)
{
  char *parse_pt = QIO_string_ptr(file_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDPropFileInfoWrapper wrapper = QIO_USQCD_PROPFILE_INFO_WRAPPER;
  QIO_USQCDPropFileInfo templ = QIO_USQCD_PROPFILE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize file info structure from a template */
  memcpy(file_info, &templ, sizeof(QIO_USQCDPropFileInfo));

  /* Start parsing file_string */
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdpropfileinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdpropfileinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_propfile_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&file_info->version);
    QIO_decode_as_string(tag,value_string,&file_info->type);
    QIO_decode_as_string(tag,value_string,&file_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&file_info->version);
  errors += QIO_check_string_occur(&file_info->type);
  errors += QIO_check_string_occur(&file_info->info);

  return errors;
}

void QIO_encode_usqcd_propsource_info(QIO_String *record_string, 
				      QIO_USQCDPropSourceInfo *record_info)
{
  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDPropSourceInfoWrapper wrapper = QIO_USQCD_PROPSOURCE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->spin, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into record wrapper structure */
  QIO_insert_usqcd_propsource_tag_string(&wrapper, recordinfo_tags);

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
    printf("QIO_encode_usqcd_propsource_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdpropsourceinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info,
				     QIO_String *record_string)
{
  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDPropSourceInfoWrapper wrapper = QIO_USQCD_PROPSOURCE_INFO_WRAPPER;
  QIO_USQCDPropSourceInfo templ = QIO_USQCD_PROPSOURCE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_USQCDPropSourceInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdpropsourceinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdpropsourceinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_propsource_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int   (tag,value_string,&record_info->spin);
    QIO_decode_as_int   (tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;
}

void QIO_encode_usqcd_proprecord_info(QIO_String *record_string, 
				      QIO_USQCDPropRecordInfo *record_info)
{
  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDPropRecordInfoWrapper wrapper = QIO_USQCD_PROPRECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->spin, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into record wrapper structure */
  QIO_insert_usqcd_proprecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("QIO_encode_usqcd_proprecord_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdproprecordinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info,
				     QIO_String *record_string)
{
  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDPropRecordInfoWrapper wrapper = QIO_USQCD_PROPRECORD_INFO_WRAPPER;
  QIO_USQCDPropRecordInfoWrapper wrapper_legacy = 
    QIO_USQCD_PROPRECORD_INFO_WRAPPER_LEGACY;
  QIO_USQCDPropRecordInfo templ = QIO_USQCD_PROPRECORD_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_USQCDPropRecordInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdproprecordinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdproprecordinfo_tags) == 0) {

    /* Otherwise start parsing the string of tags */
    parse_pt = QIO_get_usqcd_proprecord_info_tag_string(&wrapper);

  } else { 

    /* Try to decode with the old tagname */
    QIO_decode_as_string (tag, tags_string, 
			  &wrapper_legacy.usqcdproprecordinfo_tags);
    /* If outer wrapper has bad tag, exit with error status */
    if(QIO_check_string_occur(&wrapper_legacy.usqcdproprecordinfo_tags) ) {
      return QIO_BAD_XML;
    }
    
    /* Otherwise start parsing the string of tags */
    parse_pt = QIO_get_usqcd_proprecord_info_tag_string(&wrapper_legacy);
  }

  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);
    
    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int   (tag,value_string,&record_info->spin);
    QIO_decode_as_int   (tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }
  
  /* Check for completeness */
  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;
}

/********************************************************************/
/* Support for USQCD propagator file info */

/* Accessors */

/* Return integer code or negative value for failure */
int QIO_get_usqcd_propfile_type(QIO_USQCDPropFileInfo *file_info)
{
  char *string = file_info->type.value;

  if(strcmp(string,QIO_USQCDPROPFILETYPESTRING_C1D12) == 0)
    return QIO_USQCDPROPFILETYPE_C1D12;
  else if(strcmp(string,QIO_USQCDPROPFILETYPESTRING_DD_PAIRS) == 0)
    return QIO_USQCDPROPFILETYPE_DD_PAIRS;
  else if(strcmp(string,QIO_USQCDPROPFILETYPESTRING_CD_PAIRS) == 0)
    return QIO_USQCDPROPFILETYPE_CD_PAIRS;
  else if(strcmp(string,QIO_USQCDPROPFILETYPESTRING_LHPC) == 0)
    return QIO_USQCDPROPFILETYPE_LHPC;
  else
    return QIO_ERR_FILE_INFO;
}

char *QIO_get_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info)
{
  return file_info->info.value;
}

int QIO_defined_usqcd_propfile_type(QIO_USQCDPropFileInfo *file_info)
{
  return file_info->type.occur;
}

int QIO_defined_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info)
{
  return file_info->info.occur;
}

char *QIO_get_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->info.value;
}

int QIO_defined_usqcd_propsource_spin(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->spin.occur;
}

int QIO_get_usqcd_propsource_spin(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->spin.value;
}

int QIO_defined_usqcd_propsource_color(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->color.occur;
}

int QIO_get_usqcd_propsource_color(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->color.value;
}

int QIO_defined_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info)
{
  return record_info->info.occur;
}

/* added EES*/
int QIO_defined_usqcd_proprecord_spin(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->spin.occur;
}

int QIO_get_usqcd_proprecord_spin(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->spin.value;
}

/* added EES*/
int QIO_defined_usqcd_proprecord_color(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->color.occur;
}

int QIO_get_usqcd_proprecord_color(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->color.value;
}

char *QIO_get_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->info.value;
}

int QIO_defined_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info)
{
  return record_info->info.occur;
}

int QIO_insert_usqcd_propfile_version(QIO_USQCDPropFileInfo *file_info, char *version)
{
  file_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(file_info->version.value, version, QIO_MAXVALUESTRING-1);
  file_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

/* Takes one of the integer type codes and translates it to the string */
int QIO_insert_usqcd_propfile_type( QIO_USQCDPropFileInfo *file_info, int type)
{
  file_info->type.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  switch (type) {
  case QIO_USQCDPROPFILETYPE_C1D12:
    strncpy(file_info->type.value, 
	    QIO_USQCDPROPFILETYPESTRING_C1D12, QIO_MAXVALUESTRING-1);
    break;
  case QIO_USQCDPROPFILETYPE_DD_PAIRS:
    strncpy(file_info->type.value, 
	    QIO_USQCDPROPFILETYPESTRING_DD_PAIRS, QIO_MAXVALUESTRING-1);
    break;
  case QIO_USQCDPROPFILETYPE_CD_PAIRS:
    strncpy(file_info->type.value, 
	    QIO_USQCDPROPFILETYPESTRING_CD_PAIRS, QIO_MAXVALUESTRING-1);
    break;
  case QIO_USQCDPROPFILETYPE_LHPC:
    strncpy(file_info->type.value, 
	    QIO_USQCDPROPFILETYPESTRING_LHPC, QIO_MAXVALUESTRING-1);
    break;
  default:
    printf("QIO_insert_usqcd_propfile_type: Unknown type %d\n",type);
    return QIO_BAD_ARG;
  }
  file_info->type.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->type.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_propfile_info( QIO_USQCDPropFileInfo *file_info, char *info)
{
  file_info->info.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  strncpy(file_info->info.value, info, QIO_MAXVALUESTRING-1);
  file_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_propfile_tag_string(QIO_USQCDPropFileInfoWrapper *wrapper,
                                 char *fileinfo_tags){
  wrapper->usqcdpropfileinfo_tags.occur = 0;
  if(!fileinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdpropfileinfo_tags.value, fileinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdpropfileinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdpropfileinfo_tags.occur = 1;
  if(strlen(fileinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_propfile_info_tag_string(QIO_USQCDPropFileInfoWrapper *wrapper){
  return wrapper->usqcdpropfileinfo_tags.value;
}

int QIO_insert_usqcd_propsource_version(QIO_USQCDPropSourceInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_propsource_spin( QIO_USQCDPropSourceInfo *record_info, int spin)
{
  record_info->spin.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->spin.value =  spin;
  record_info->spin.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_propsource_color( QIO_USQCDPropSourceInfo *record_info, int color)
{
  record_info->color.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->color.value =  color;
  record_info->color.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_propsource_info( QIO_USQCDPropSourceInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_propsource_tag_string(QIO_USQCDPropSourceInfoWrapper *wrapper,
                                 char *recordinfo_tags){
  wrapper->usqcdpropsourceinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdpropsourceinfo_tags.value, recordinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdpropsourceinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdpropsourceinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_propsource_info_tag_string(QIO_USQCDPropSourceInfoWrapper *wrapper){
  return wrapper->usqcdpropsourceinfo_tags.value;
}

int QIO_insert_usqcd_proprecord_version(QIO_USQCDPropRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_proprecord_spin( QIO_USQCDPropRecordInfo *record_info, int spin)
{
  record_info->spin.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->spin.value =  spin;
  record_info->spin.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_proprecord_color( QIO_USQCDPropRecordInfo *record_info, int color)
{
  record_info->color.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->color.value =  color;
  record_info->color.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_proprecord_info( QIO_USQCDPropRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_proprecord_tag_string(QIO_USQCDPropRecordInfoWrapper *wrapper,
                                 char *recordinfo_tags){
  wrapper->usqcdproprecordinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdproprecordinfo_tags.value, recordinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdproprecordinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdproprecordinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_proprecord_info_tag_string(QIO_USQCDPropRecordInfoWrapper *wrapper){
  return wrapper->usqcdproprecordinfo_tags.value;
}

QIO_USQCDPropFileInfo *QIO_create_usqcd_propfile_info(int type, char *info)
{
  
  QIO_USQCDPropFileInfo templ = QIO_USQCD_PROPFILE_INFO_TEMPLATE;
  QIO_USQCDPropFileInfo *file_info;

  file_info = (QIO_USQCDPropFileInfo *)malloc(sizeof(QIO_USQCDPropFileInfo));
  if(!file_info)return NULL;

  memcpy(file_info, &templ, sizeof(QIO_USQCDPropFileInfo));
  QIO_insert_usqcd_propfile_version(file_info,QIO_USQCDPROPFILEFORMATVERSION);  
  
  QIO_insert_usqcd_propfile_type( file_info, type);
  QIO_insert_usqcd_propfile_info( file_info, info);

  return file_info;
}

void QIO_destroy_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info){
  free(file_info);
}

QIO_USQCDPropSourceInfo *QIO_create_usqcd_propsource_info(char *info)
{
  
  QIO_USQCDPropSourceInfo templ = QIO_USQCD_PROPSOURCE_INFO_TEMPLATE;
  QIO_USQCDPropSourceInfo *record_info;

  record_info = (QIO_USQCDPropSourceInfo *)malloc(sizeof(QIO_USQCDPropSourceInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDPropSourceInfo));
  QIO_insert_usqcd_propsource_version(record_info,QIO_USQCDPROPSOURCEFORMATVERSION);  
  
  QIO_insert_usqcd_propsource_info( record_info, info);

  return record_info;
}

/* Encoding when spin and color are required */

QIO_USQCDPropSourceInfo *QIO_create_usqcd_propsource_sc_info(int spin, 
     int color, char *info)
{
  
  QIO_USQCDPropSourceInfo templ = QIO_USQCD_PROPSOURCE_INFO_TEMPLATE;
  QIO_USQCDPropSourceInfo *record_info;

  record_info = (QIO_USQCDPropSourceInfo *)malloc(sizeof(QIO_USQCDPropSourceInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDPropSourceInfo));
  QIO_insert_usqcd_propsource_version(record_info,QIO_USQCDPROPSOURCEFORMATVERSION);  
  
  QIO_insert_usqcd_propsource_spin( record_info, spin);
  QIO_insert_usqcd_propsource_color( record_info, color);
  QIO_insert_usqcd_propsource_info( record_info, info);

  return record_info;
}

void QIO_destroy_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info){
  free(record_info);
}


/* Encoding when spin and color are not appropriate */

QIO_USQCDPropRecordInfo *QIO_create_usqcd_proprecord_info(char *info)
{
  
  QIO_USQCDPropRecordInfo templ = QIO_USQCD_PROPRECORD_INFO_TEMPLATE;
  QIO_USQCDPropRecordInfo *record_info;

  record_info = (QIO_USQCDPropRecordInfo *)malloc(sizeof(QIO_USQCDPropRecordInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDPropRecordInfo));
  QIO_insert_usqcd_proprecord_version(record_info,QIO_USQCDPROPRECORDFORMATVERSION);  
  
  QIO_insert_usqcd_proprecord_info( record_info, info);

  return record_info;
}

/* Encoding when spin and color are required */

QIO_USQCDPropRecordInfo *QIO_create_usqcd_proprecord_sc_info(int spin, 
     int color, char *info)
{
  
  QIO_USQCDPropRecordInfo templ = QIO_USQCD_PROPRECORD_INFO_TEMPLATE;
  QIO_USQCDPropRecordInfo *record_info;

  record_info = (QIO_USQCDPropRecordInfo *)malloc(sizeof(QIO_USQCDPropRecordInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDPropRecordInfo));
  QIO_insert_usqcd_proprecord_version(record_info,QIO_USQCDPROPRECORDFORMATVERSION);  
  
  QIO_insert_usqcd_proprecord_spin( record_info, spin);
  QIO_insert_usqcd_proprecord_color( record_info, color);
  QIO_insert_usqcd_proprecord_info( record_info, info);

  return record_info;
}

void QIO_destroy_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info){
  free(record_info);
}




