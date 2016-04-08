/* Read and write record info for SciDAC/USQCD lattice file format */

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


void QIO_encode_usqcd_lattice_info(QIO_String *record_string, 
				     QIO_USQCDLatticeInfo *record_info)
{
  /* taken from QIO_encode_record_info */

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDLatticeInfoWrapper wrapper = QIO_USQCD_LATTICE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->plaq, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->linktr, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_usqcdlattice_tag_string(&wrapper, recordinfo_tags);

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
    printf("QIO_encode_usqcd_lattice_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdlatticeinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_lattice_info(QIO_USQCDLatticeInfo *record_info,
				    QIO_String *record_string)
{

  /* taken from QIO_decode_record_info */

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDLatticeInfoWrapper wrapper = QIO_USQCD_LATTICE_INFO_WRAPPER;
  QIO_USQCDLatticeInfo templ = QIO_USQCD_LATTICE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  memcpy(record_info, &templ, sizeof(QIO_USQCDLatticeInfo));
  
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
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdlatticeinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdlatticeinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_lattice_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_string(tag,value_string,&record_info->plaq);
    QIO_decode_as_string(tag,value_string,&record_info->linktr);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->plaq);
  errors += QIO_check_string_occur(&record_info->linktr);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}

/* Accessors for user record info in a USQCD lattice file */


char *QIO_get_plaq(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->plaq.value;
}

char *QIO_get_linktr(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->linktr.value;
}

char *QIO_get_info(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->info.value;
}

int QIO_defined_plaq(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->plaq.occur;
}

int QIO_defined_linktr(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->linktr.occur;
}

int QIO_defined_info(QIO_USQCDLatticeInfo *record_info)
{
  return record_info->info.occur;
}

int QIO_insert_usqcdlattice_version(QIO_USQCDLatticeInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}


int QIO_insert_usqcdlattice_plaq( QIO_USQCDLatticeInfo *record_info, char *plaq)
{
  record_info->plaq.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->plaq.value, plaq, QIO_MAXVALUESTRING-1);
  record_info->plaq.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->plaq.occur = 1;
  if(strlen(plaq) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcdlattice_linktr( QIO_USQCDLatticeInfo *record_info, char *linktr)
{
  record_info->linktr.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->linktr.value, linktr, QIO_MAXVALUESTRING-1);
  record_info->linktr.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->linktr.occur = 1;
  if(strlen(linktr) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcdlattice_info( QIO_USQCDLatticeInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int QIO_insert_usqcdlattice_tag_string(QIO_USQCDLatticeInfoWrapper *wrapper,
                                 char *recordinfo_tags){
  wrapper->usqcdlatticeinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdlatticeinfo_tags.value, recordinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdlatticeinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdlatticeinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_lattice_info_tag_string(QIO_USQCDLatticeInfoWrapper *wrapper){
  return wrapper->usqcdlatticeinfo_tags.value;
}

QIO_USQCDLatticeInfo *QIO_create_usqcd_lattice_info(char *plaq, char *linktr, char *info)
{
  
  /* taken from QIO_create_record_info and modified */

  QIO_USQCDLatticeInfo templ = QIO_USQCD_LATTICE_INFO_TEMPLATE;
  QIO_USQCDLatticeInfo *record_info;
  time_t cu_time;

  record_info = (QIO_USQCDLatticeInfo *)malloc(sizeof(QIO_USQCDLatticeInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  memcpy(record_info, &templ, sizeof(QIO_USQCDLatticeInfo));
  QIO_insert_usqcdlattice_version(record_info,QIO_USQCDLATTICEFORMATVERSION);  
  
  QIO_insert_usqcdlattice_plaq( record_info, plaq);
  QIO_insert_usqcdlattice_linktr( record_info, linktr);
  QIO_insert_usqcdlattice_info( record_info, info);

  return record_info;

}

void QIO_destroy_usqcd_lattice_info(QIO_USQCDLatticeInfo *record_info){
  free(record_info);
}
