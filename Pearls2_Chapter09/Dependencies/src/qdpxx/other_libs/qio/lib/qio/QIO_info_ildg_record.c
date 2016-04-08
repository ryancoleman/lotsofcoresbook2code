/* Read and write ILDG format info for LIME/ILDG binary lattice format */

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


int QIO_decode_ILDG_format_info(QIO_ILDGFormatInfo *ildg_info, 
				QIO_String *ildg_string){
  char *parse_pt = QIO_string_ptr(ildg_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_ILDGFormatInfoWrapper wrapper = QIO_ILDG_FORMAT_INFO_WRAPPER;
  QIO_ILDGFormatInfo templ = QIO_ILDG_FORMAT_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize ILDG format info structure from a template */
  memcpy(ildg_info, &templ, sizeof(QIO_ILDGFormatInfo));

  /* Start parsing ILDG format_string */
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
  QIO_decode_as_string (tag, tags_string, &wrapper.ildgformatinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.ildgformatinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_ildgformat_info_tag_string(&wrapper);

  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);
    
    QIO_decode_as_string(tag,value_string,&ildg_info->version);
    QIO_decode_as_string(tag,value_string,&ildg_info->field);
    QIO_decode_as_int   (tag,value_string,&ildg_info->precision);
    QIO_decode_as_int   (tag,value_string,&ildg_info->lx);
    QIO_decode_as_int   (tag,value_string,&ildg_info->ly);
    QIO_decode_as_int   (tag,value_string,&ildg_info->lz);
    QIO_decode_as_int   (tag,value_string,&ildg_info->lt);
  }

  /* Check for completeness */
  
  errors += QIO_check_string_occur(&ildg_info->version);
  errors += QIO_check_string_occur(&ildg_info->field);
  errors += QIO_check_int_occur   (&ildg_info->precision);
  errors += QIO_check_int_occur   (&ildg_info->lx);
  errors += QIO_check_int_occur   (&ildg_info->ly);
  errors += QIO_check_int_occur   (&ildg_info->lz);
  errors += QIO_check_int_occur   (&ildg_info->lt);

  return errors;
}

void QIO_encode_ILDG_format_info(QIO_String *ildg_string, 
				 QIO_ILDGFormatInfo *ildg_info){
  char *buf;
  int remainder,n;
  char ildgformatinfo_tags[QIO_MAXVALUESTRING];
  QIO_ILDGFormatInfoWrapper wrapper = QIO_ILDG_FORMAT_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = ildgformatinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&ildg_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&ildg_info->field, &remainder);
  buf = QIO_encode_as_int   (buf,&ildg_info->precision, &remainder);
  buf = QIO_encode_as_int   (buf,&ildg_info->lx, &remainder);
  buf = QIO_encode_as_int   (buf,&ildg_info->ly, &remainder);
  buf = QIO_encode_as_int   (buf,&ildg_info->lz, &remainder);
  buf = QIO_encode_as_int   (buf,&ildg_info->lt, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_ildgformat_tag_string(&wrapper, ildgformatinfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(ildg_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(ildg_string);
  remainder = QIO_string_length(ildg_string);
  
  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("QIO_encode_ildg_format_info: ildg_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.ildgformatinfo_tags, &remainder);
  }
}


/* Utilities for ILDG format info values */

int QIO_insert_ildgformat_tag_string(QIO_ILDGFormatInfoWrapper *wrapper, 
				      char *ildgformatinfo_tags){
  wrapper->ildgformatinfo_tags.occur = 0;
  if(!ildgformatinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->ildgformatinfo_tags.value, ildgformatinfo_tags, 
	  QIO_MAXVALUESTRING-1);
  wrapper->ildgformatinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->ildgformatinfo_tags.occur = 1;
  if(strlen(ildgformatinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_ildgformat_version(QIO_ILDGFormatInfo *ildg_info, 
				   char *version){
  ildg_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(ildg_info->version.value, version, QIO_MAXVALUESTRING-1);
  ildg_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  ildg_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_ildgformat_field(QIO_ILDGFormatInfo *ildg_info, 
				char *field_string){
  ildg_info->field.occur = 0;
  if(field_string == NULL)return QIO_BAD_ARG;
  strncpy(ildg_info->field.value, field_string, QIO_MAXVALUESTRING-1);
  ildg_info->field.value[QIO_MAXVALUESTRING-1] = '\0';
  ildg_info->field.occur = 1;
  if(strlen(field_string) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_ildgformat_precision(QIO_ILDGFormatInfo *ildg_info, int precision){
  ildg_info->precision.occur = 0;
  if(!ildg_info)return QIO_BAD_ARG;
  ildg_info->precision.value = precision;
  ildg_info->precision.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_ildgformat_lx(QIO_ILDGFormatInfo *ildg_info, int lx){
  ildg_info->lx.occur = 0;
  if(!ildg_info)return QIO_BAD_ARG;
  ildg_info->lx.value = lx;
  ildg_info->lx.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_ildgformat_ly(QIO_ILDGFormatInfo *ildg_info, int ly){
  ildg_info->ly.occur = 0;
  if(!ildg_info)return QIO_BAD_ARG;
  ildg_info->ly.value = ly;
  ildg_info->ly.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_ildgformat_lz(QIO_ILDGFormatInfo *ildg_info, int lz){
  ildg_info->lz.occur = 0;
  if(!ildg_info)return QIO_BAD_ARG;
  ildg_info->lz.value = lz;
  ildg_info->lz.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_ildgformat_lt(QIO_ILDGFormatInfo *ildg_info, int lt){
  ildg_info->lt.occur = 0;
  if(!ildg_info)return QIO_BAD_ARG;
  ildg_info->lt.value = lt;
  ildg_info->lt.occur = 1;
  return QIO_SUCCESS;
}

/* Accessors for ILDG record format info */

char *QIO_get_ildgformat_info_tag_string(QIO_ILDGFormatInfoWrapper *wrapper){
  return wrapper->ildgformatinfo_tags.value;
}

char *QIO_get_ildgformat_field(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->field.value;
}

int QIO_get_ildgformat_precision(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->precision.value;
}

int QIO_get_ildgformat_lx(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->lx.value;
}

int QIO_get_ildgformat_ly(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->ly.value;
}

int QIO_get_ildgformat_lz(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->lz.value;
}

int QIO_get_ildgformat_lt(QIO_ILDGFormatInfo *ildg_info){
  return ildg_info->lt.value;
}


QIO_ILDGFormatInfo *QIO_create_ildg_format_info(int precision, int *dims){
  QIO_ILDGFormatInfo templ = QIO_ILDG_FORMAT_INFO_TEMPLATE;
  QIO_ILDGFormatInfo *ildg_info;
  
  ildg_info = (QIO_ILDGFormatInfo *)malloc(sizeof(QIO_ILDGFormatInfo));
  if(!ildg_info)return NULL;

  memcpy(ildg_info, &templ, sizeof(QIO_ILDGFormatInfo));
  QIO_insert_ildgformat_version(ildg_info,QIO_ILDGFORMATVERSION);
  QIO_insert_ildgformat_field(ildg_info,"su3gauge");
  if(precision != 0)
    QIO_insert_ildgformat_precision(ildg_info,precision);
  if(dims != NULL){
    QIO_insert_ildgformat_lx(ildg_info,dims[0]);
    QIO_insert_ildgformat_ly(ildg_info,dims[1]);
    QIO_insert_ildgformat_lz(ildg_info,dims[2]);
    QIO_insert_ildgformat_lt(ildg_info,dims[3]);
  }
  return ildg_info;
}

void QIO_destroy_ildg_format_info(QIO_ILDGFormatInfo *ildg_info){
  free(ildg_info);
}

