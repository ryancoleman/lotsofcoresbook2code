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


void QIO_encode_usqcd_kspropfile_info(QIO_String *file_string, 
				      QIO_USQCDKSPropFileInfo *file_info)
{
  char *buf;
  int remainder,n;
  char fileinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDKSPropFileInfoWrapper wrapper = QIO_USQCD_KSPROPFILE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = fileinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&file_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->type, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_usqcdkspropfile_tag_string(&wrapper, fileinfo_tags);

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
    printf("QIO_encode_usqcd_kspropfile_info: file_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdkspropfileinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info,
				   QIO_String *file_string)
{
  char *parse_pt = QIO_string_ptr(file_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDKSPropFileInfoWrapper wrapper = QIO_USQCD_KSPROPFILE_INFO_WRAPPER;
  QIO_USQCDKSPropFileInfo templ = QIO_USQCD_KSPROPFILE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize file info structure from a template */
  memcpy(file_info, &templ, sizeof(QIO_USQCDKSPropFileInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdkspropfileinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdkspropfileinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_kspropfile_info_tag_string(&wrapper);
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

void QIO_encode_usqcd_kspropsource_info(QIO_String *file_string, 
				    QIO_USQCDKSPropSourceInfo *file_info)
{
  char *buf;
  int remainder,n;
  char fileinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDKSPropSourceInfoWrapper wrapper = QIO_USQCD_KSPROPSOURCE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = fileinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&file_info->version, &remainder);
  buf = QIO_encode_as_int   (buf,&file_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_usqcdkspropsource_tag_string(&wrapper, fileinfo_tags);

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
    printf("QIO_encode_usqcd_kspropsource_info: file_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdkspropsourceinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *record_info,
				   QIO_String *record_string)
{
  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDKSPropSourceInfoWrapper wrapper = QIO_USQCD_KSPROPSOURCE_INFO_WRAPPER;
  QIO_USQCDKSPropSourceInfo templ = QIO_USQCD_KSPROPSOURCE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_USQCDKSPropSourceInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdkspropsourceinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdkspropsourceinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_kspropsource_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int   (tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;
}

void QIO_encode_usqcd_ksproprecord_info(QIO_String *record_string, 
				    QIO_USQCDKSPropRecordInfo *record_info)
{
  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDKSPropRecordInfoWrapper wrapper = QIO_USQCD_KSPROPRECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_int   (buf,&record_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into record wrapper structure */
  QIO_insert_usqcd_ksproprecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("QIO_encode_usqcd_ksproprecord_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdksproprecordinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info,
				   QIO_String *record_string)
{
  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDKSPropRecordInfoWrapper wrapper = QIO_USQCD_KSPROPRECORD_INFO_WRAPPER;
  QIO_USQCDKSPropRecordInfo templ = QIO_USQCD_KSPROPRECORD_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_USQCDKSPropRecordInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdksproprecordinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdksproprecordinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_ksproprecord_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int   (tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;
}

/********************************************************************/
/* Support for USQCD KS propagator file info */

/* Accessors */

/* Return integer code or negative value for failure */
int QIO_get_usqcd_kspropfile_type(QIO_USQCDKSPropFileInfo *file_info)
{
  char *string = file_info->type.value;

  if(strcmp(string,QIO_USQCDKSPROPFILETYPESTRING_C1V3) == 0)
    return QIO_USQCDKSPROPFILETYPE_C1V3;
  else if(strcmp(string,QIO_USQCDKSPROPFILETYPESTRING_VV_PAIRS) == 0)
    return QIO_USQCDKSPROPFILETYPE_VV_PAIRS;
  else if(strcmp(string,QIO_USQCDKSPROPFILETYPESTRING_CV_PAIRS) == 0)
    return QIO_USQCDKSPROPFILETYPE_CV_PAIRS;
  else
    return QIO_ERR_FILE_INFO;
}

char *QIO_get_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info)
{
  return file_info->info.value;
}

int QIO_defined_usqcd_kspropfile_type(QIO_USQCDKSPropFileInfo *file_info)
{
  return file_info->type.occur;
}

int QIO_defined_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info)
{
  return file_info->info.occur;
}

char *QIO_get_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *record_info)
{
  return record_info->info.value;
}

int QIO_defined_usqcd_kspropsource_color(QIO_USQCDKSPropSourceInfo *source_info)
{
  return source_info->color.occur;
}

int QIO_get_usqcd_kspropsource_color(QIO_USQCDKSPropSourceInfo *source_info)
{
  return source_info->color.value;
}

int QIO_defined_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *source_info)
{
  return source_info->info.occur;
}

int QIO_defined_usqcd_ksproprecord_color(QIO_USQCDKSPropRecordInfo *record_info)
{
  return record_info->color.occur;
}

int QIO_get_usqcd_ksproprecord_color(QIO_USQCDKSPropRecordInfo *record_info)
{
  return record_info->color.value;
}

char *QIO_get_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info)
{
  return record_info->info.value;
}

int QIO_defined_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info)
{
  return record_info->info.occur;
}

int QIO_insert_usqcdkspropfile_version(QIO_USQCDKSPropFileInfo *file_info, char *version)
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
int QIO_insert_usqcdkspropfile_type( QIO_USQCDKSPropFileInfo *file_info, int type)
{
  file_info->type.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  switch (type) {
  case QIO_USQCDKSPROPFILETYPE_C1V3:
    strncpy(file_info->type.value, 
	    QIO_USQCDKSPROPFILETYPESTRING_C1V3, QIO_MAXVALUESTRING-1);
    break;
  case QIO_USQCDKSPROPFILETYPE_VV_PAIRS:
    strncpy(file_info->type.value, 
	    QIO_USQCDKSPROPFILETYPESTRING_VV_PAIRS, QIO_MAXVALUESTRING-1);
    break;
  case QIO_USQCDKSPROPFILETYPE_CV_PAIRS:
    strncpy(file_info->type.value, 
	    QIO_USQCDKSPROPFILETYPESTRING_CV_PAIRS, QIO_MAXVALUESTRING-1);
    break;
  default:
    printf("QIO_insert_usqcdkspropfile_type: Unknown type %d\n",type);
    return QIO_BAD_ARG;
  }
  file_info->type.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->type.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcdkspropfile_info( QIO_USQCDKSPropFileInfo *file_info, char *info)
{
  file_info->info.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  strncpy(file_info->info.value, info, QIO_MAXVALUESTRING-1);
  file_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcdkspropfile_tag_string(QIO_USQCDKSPropFileInfoWrapper *wrapper,
                                 char *fileinfo_tags){
  wrapper->usqcdkspropfileinfo_tags.occur = 0;
  if(!fileinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdkspropfileinfo_tags.value, fileinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdkspropfileinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdkspropfileinfo_tags.occur = 1;
  if(strlen(fileinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_kspropfile_info_tag_string(QIO_USQCDKSPropFileInfoWrapper *wrapper){
  return wrapper->usqcdkspropfileinfo_tags.value;
}

int QIO_insert_usqcdkspropsource_version(QIO_USQCDKSPropSourceInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcdkspropsource_info( QIO_USQCDKSPropSourceInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_kspropsource_version(QIO_USQCDKSPropSourceInfo *source_info, char *version)
{
  source_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(source_info->version.value, version, QIO_MAXVALUESTRING-1);
  source_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  source_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_kspropsource_color( QIO_USQCDKSPropSourceInfo *source_info, int color)
{
  source_info->color.occur = 0;
  if(!source_info)return QIO_BAD_ARG;
  source_info->color.value =  color;
  source_info->color.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcdkspropsource_tag_string(QIO_USQCDKSPropSourceInfoWrapper *wrapper,
                                 char *sourceinfo_tags){
  wrapper->usqcdkspropsourceinfo_tags.occur = 0;
  if(!sourceinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdkspropsourceinfo_tags.value, sourceinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdkspropsourceinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdkspropsourceinfo_tags.occur = 1;
  if(strlen(sourceinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_kspropsource_info( QIO_USQCDKSPropSourceInfo *source_info, char *info)
{
  source_info->info.occur = 0;
  if(!source_info)return QIO_BAD_ARG;
  strncpy(source_info->info.value, info, QIO_MAXVALUESTRING-1);
  source_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  source_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_kspropsource_info_tag_string(QIO_USQCDKSPropSourceInfoWrapper *wrapper){
  return wrapper->usqcdkspropsourceinfo_tags.value;
}

int QIO_insert_usqcd_ksproprecord_version(QIO_USQCDKSPropRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_ksproprecord_color( QIO_USQCDKSPropRecordInfo *record_info, int color)
{
  record_info->color.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  record_info->color.value =  color;
  record_info->color.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcd_ksproprecord_info( QIO_USQCDKSPropRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcd_ksproprecord_tag_string(QIO_USQCDKSPropRecordInfoWrapper *wrapper,
                                 char *recordinfo_tags){
  wrapper->usqcdksproprecordinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdksproprecordinfo_tags.value, recordinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdksproprecordinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdksproprecordinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_ksproprecord_info_tag_string(QIO_USQCDKSPropRecordInfoWrapper *wrapper){
  return wrapper->usqcdksproprecordinfo_tags.value;
}

QIO_USQCDKSPropFileInfo *QIO_create_usqcd_kspropfile_info(int type, char *info)
{
  
  QIO_USQCDKSPropFileInfo templ = QIO_USQCD_KSPROPFILE_INFO_TEMPLATE;
  QIO_USQCDKSPropFileInfo *file_info;

  file_info = (QIO_USQCDKSPropFileInfo *)malloc(sizeof(QIO_USQCDKSPropFileInfo));
  if(!file_info)return NULL;

  memcpy(file_info, &templ, sizeof(QIO_USQCDKSPropFileInfo));
  QIO_insert_usqcdkspropfile_version(file_info,QIO_USQCDKSPROPFILEFORMATVERSION);  
  
  QIO_insert_usqcdkspropfile_type( file_info, type);
  QIO_insert_usqcdkspropfile_info( file_info, info);

  return file_info;
}

void QIO_destroy_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info){
  free(file_info);
}

QIO_USQCDKSPropSourceInfo *QIO_create_usqcd_kspropsource_info(char *info)
{
  
  QIO_USQCDKSPropSourceInfo templ = QIO_USQCD_KSPROPSOURCE_INFO_TEMPLATE;
  QIO_USQCDKSPropSourceInfo *source_info;

  source_info = (QIO_USQCDKSPropSourceInfo *)malloc(sizeof(QIO_USQCDKSPropSourceInfo));
  if(!source_info)return NULL;

  memcpy(source_info, &templ, sizeof(QIO_USQCDKSPropSourceInfo));
  QIO_insert_usqcdkspropsource_version(source_info,QIO_USQCDKSPROPSOURCEFORMATVERSION);  
  
  QIO_insert_usqcdkspropsource_info( source_info, info);

  return source_info;
}

/* Encoding when color is required */

QIO_USQCDKSPropSourceInfo *QIO_create_usqcd_kspropsource_c_info(
     int color, char *info)
{
  
  QIO_USQCDKSPropSourceInfo templ = QIO_USQCD_KSPROPSOURCE_INFO_TEMPLATE;
  QIO_USQCDKSPropSourceInfo *source_info;

  source_info = (QIO_USQCDKSPropSourceInfo *)malloc(sizeof(QIO_USQCDKSPropSourceInfo));
  if(!source_info)return NULL;

  memcpy(source_info, &templ, sizeof(QIO_USQCDKSPropSourceInfo));
  QIO_insert_usqcd_kspropsource_version(source_info,QIO_USQCDKSPROPSOURCEFORMATVERSION);  
  
  QIO_insert_usqcd_kspropsource_color( source_info, color);
  QIO_insert_usqcd_kspropsource_info( source_info, info);

  return source_info;
}

void QIO_destroy_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *source_info){
  free(source_info);
}


/* Encoding when color is not appropriate */

QIO_USQCDKSPropRecordInfo *QIO_create_usqcd_ksproprecord_info(char *info)
{
  
  QIO_USQCDKSPropRecordInfo templ = QIO_USQCD_KSPROPRECORD_INFO_TEMPLATE;
  QIO_USQCDKSPropRecordInfo *record_info;

  record_info = (QIO_USQCDKSPropRecordInfo *)malloc(sizeof(QIO_USQCDKSPropRecordInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDKSPropRecordInfo));
  QIO_insert_usqcd_ksproprecord_version(record_info,QIO_USQCDKSPROPRECORDFORMATVERSION);  
  
  QIO_insert_usqcd_ksproprecord_info( record_info, info);

  return record_info;
}

/* Encoding when color is required */

QIO_USQCDKSPropRecordInfo *QIO_create_usqcd_ksproprecord_c_info(
     int color, char *info)
{
  
  QIO_USQCDKSPropRecordInfo templ = QIO_USQCD_KSPROPRECORD_INFO_TEMPLATE;
  QIO_USQCDKSPropRecordInfo *record_info;

  record_info = (QIO_USQCDKSPropRecordInfo *)malloc(sizeof(QIO_USQCDKSPropRecordInfo));
  if(!record_info)return NULL;

  memcpy(record_info, &templ, sizeof(QIO_USQCDKSPropRecordInfo));
  QIO_insert_usqcd_ksproprecord_version(record_info,QIO_USQCDKSPROPRECORDFORMATVERSION);  
  
  QIO_insert_usqcd_ksproprecord_color( record_info, color);
  QIO_insert_usqcd_ksproprecord_info( record_info, info);

  return record_info;
}

void QIO_destroy_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info){
  free(record_info);
}




