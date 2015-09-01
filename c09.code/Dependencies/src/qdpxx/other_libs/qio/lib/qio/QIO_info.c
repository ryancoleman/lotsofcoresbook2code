/* Read and write private record, file, and checksum XML for SciDAC
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

/* Same as strpbrk: Find next occurrence of any character in tokens */
char *QIO_next_token(char *parse_pt, char *tokens){
  char *t;
  int found;
  for(; *parse_pt != '\0'; parse_pt++){
    found = 0;
    for(t = tokens; *t != '\0'; t++)if(*parse_pt == *t){
      found++;
      break;
    }
    if(found)break;
  }
  return parse_pt;
}

/* Opposite of strpbrk: Find next occurrence of any character not in tokens */
char *QIO_next_nontoken(char *parse_pt, char *tokens){
  char *t;
  int found;
  for(; *parse_pt != '\0'; parse_pt++){
    found = 0;
    for(t = tokens; *t != '\0'; t++)if(*parse_pt == *t){
      found++;
      break;
    }
    if(!found)break;
  }
  return parse_pt;
}

/* Almost same as strncat, but exits if n < 0 and decrements n
   for a successful copy */
char *QIO_strncat(char *s1, char *s2, int *n){
  char *c;
  if(*n < 0)return s1;
  c = strncat(s1,s2,*n);
  *n -= strlen(s2);
  if(*n < 0)*n = 0;
  return c;
}

/* Find the next <tag> or </tag> or <tag/> */
/* On return left_angle points to open "<" and return value
   to end ">" or to end of string, whichever comes first */
char *QIO_next_tag(char *parse_pt, char *tag, char **left_angle){
  char *begin_tag;
  int n;

  /* Initially results are null strings */
  tag[0] = '\0';

  /* Look for next "<" */
  parse_pt = QIO_next_token(parse_pt, "<");

  /* Mark left angle or '\0' */
  *left_angle = parse_pt;

  /* Exit if end of string was reached before finding < */
  if(!*parse_pt)return parse_pt;

  /* Move past '<' */
  parse_pt++;

  /* Tag starts at first nonwhite */
  parse_pt = QIO_next_nontoken(parse_pt, " \t\n\r");

  begin_tag = parse_pt;

  /* Exit if end of string was reached before finding start of tag */
  if(!*parse_pt)return parse_pt;

  /* Tag ends at next white or > */
  parse_pt = QIO_next_token(parse_pt, " >\t\n\r");

  /* Count characters in tag */
  n = parse_pt - begin_tag;
  
  /* Truncate if oversized */
  if(n > QIO_MAXTAG-1){
    n = QIO_MAXTAG-1;
    printf("QIO_next_tag: tag too long - truncated\n");
  }

  /* Copy tag and terminate with null */
  strncpy(tag, begin_tag, n);
  tag[n] = '\0';

  /* Scan to closing '>' */
  parse_pt = QIO_next_token(parse_pt, ">");

  if(!*parse_pt)return parse_pt;
  
  /* Move past '>' */
  parse_pt++;

  return parse_pt;
}

/* Starting from parse_pt, scan for next <tag> and then find matching
   end </tag>.  Copy the enclosed string to value_string with null
   termination. */

char *QIO_get_tag_value(char *parse_pt, char *tag, char *value_string){
  char tmptag[QIO_MAXTAG];

  char *begin_tag_ptr;
  char *begin_value_string;
  char *peek_pt;
  int end_tag = 0;
  int i,n;

  /* Initially results are null strings */
  *tag = '\0';
  *value_string = '\0';

  /* Get tag */

  parse_pt = QIO_next_tag(parse_pt, tag, &begin_tag_ptr);
  if(*tag == '\0')return parse_pt;

  /* Parse the value */

  /* Handle <tag/> case, signifying empty data */

  /* Look for a slash that is not the first character in "tag" */
  peek_pt = QIO_next_token(tag, "/");
  if(*peek_pt != '\0' && peek_pt != tag ){
    /* Replace the slash in "tag" with a terminating null */
    *peek_pt = '\0';
    /* At this point "value_string" is null and "parse_pt" is past the
       closing ">", so we are finished. */
    return parse_pt; 
  }

  /* Handle <tag>value</tag> case */

  /* Value starts at first nonwhite */
  parse_pt = QIO_next_nontoken(parse_pt, " \t\n\r");

  /* Exit if end of string was reached before finding beginning of value */
  if(!*parse_pt)return parse_pt;

  /* Mark beginning of value */
  begin_value_string = parse_pt;

  /* Look for matching end tag */
  while(!end_tag){
    parse_pt = QIO_next_tag(parse_pt, tmptag, &begin_tag_ptr);
    /* Stop when matching tag or end of string found */
    end_tag = !*parse_pt || (tmptag[0] == '/' && strcmp(tmptag+1,tag)==0);
  }

  /* Copy value as string */

  n = begin_tag_ptr - begin_value_string;
  if(n > QIO_MAXVALUESTRING-1){
    n = QIO_MAXVALUESTRING - 1;
    printf("QIO_get_tag_value: string truncated");
  }
  strncpy(value_string, begin_value_string, n );
  /* Terminate the string */
  value_string[n] = '\0';

  /* Trim trailing white space from value string */
  for(i = n-1; i >= 0; i--){
    if(value_string[i] != ' ' 
       && value_string[i] != '\t'
       && value_string[i] != '\n'
       && value_string[i] != '\r' )break;
    value_string[i] = '\0';
  }
  return parse_pt;
}

/* If tag matches, set value to the string */

void QIO_decode_as_string(char *tag, char *value_string, 
			    QIO_TagCharValue *tag_value){
  int n = strlen(value_string);

  if(strcmp(tag,tag_value->tag)==0){
    if(n > QIO_MAXVALUESTRING-1){
      n = QIO_MAXVALUESTRING - 1;
      printf("QIO_decode_as_string: string truncated");
    }
    strncpy(tag_value->value,value_string,QIO_MAXVALUESTRING-1);
    tag_value->value[QIO_MAXVALUESTRING-1] = '\0';
    tag_value->attr[0] = '\0';   /* Ignore attributes for now */
    tag_value->occur++;
  }
}

/* If tag matches, set value to int conversion of string */

void QIO_decode_as_int(char *tag, char *value_string, 
			 QIO_TagIntValue *tag_value){
  if(strcmp(tag,tag_value->tag)==0){
    tag_value->value = atoi(value_string);
    tag_value->attr[0] = '\0';   /* Ignore attributes for now */
    tag_value->occur++;
  }
}

/* If tag matches, set value to hex32 conversion of string */

void QIO_decode_as_hex32(char *tag, char *value_string, 
			    QIO_TagHex32Value *tag_value){
  if(strcmp(tag,tag_value->tag)==0){
    if(sscanf(value_string,"%x",&tag_value->value)==1){
      tag_value->occur++;
      tag_value->attr[0] = '\0';   /* Ignore attributes for now */
    }
  }
}

/* If tag matches, set array value to the list of integers */

void QIO_decode_as_intlist(char *tag, char *value_string, 
			     QIO_TagIntArrayValue *tag_value){
  int i;
  char *s;
  
  if(strcmp(tag,tag_value->tag)==0){
    for(s = strtok(value_string," "), i = 0; 
	s && *s && i < QIO_MAXINTARRAY; 
	s = strtok('\0'," "),i++)
      tag_value->value[i] = atoi(s);
    
    tag_value->n = i;
    tag_value->attr[0] = '\0';   /* Ignore attributes for now */

    /* Trouble if we didn't hit the end of the string and the array is full */
    if(s && i == QIO_MAXINTARRAY){
      printf("QIO_decode_as_intlist: exceeded internal array dimensions %d\n", QIO_MAXINTARRAY);
    }
    else if(tag_value->n > 0)tag_value->occur++;
  }
}

/* Write a tag with attributes, if specified */
char *QIO_write_tag(char *buf, char *tag, char *attr, int *remainder){
  QIO_strncat(buf,"<",remainder);
  QIO_strncat(buf,tag,remainder);
  if(strlen(attr) != 0){
    QIO_strncat(buf," ",remainder);
    QIO_strncat(buf,attr,remainder);
  }
  QIO_strncat(buf,">",remainder);
  return strchr(buf,'\0');
}

char *QIO_write_endtag(char *buf,char *tag,int *remainder){
  QIO_strncat(buf,"</",remainder);
  QIO_strncat(buf,tag,remainder);
  QIO_strncat(buf,">",remainder);
  return strchr(buf,'\0');
}

char *QIO_encode_as_string(char *buf, QIO_TagCharValue *tag_value,
			     int *remainder){

  /* Don't write value unless occurs */
  if(!tag_value->occur)return buf;
  buf = QIO_write_tag(buf, tag_value->tag, tag_value->attr, remainder);
  snprintf(buf,*remainder,"%s",tag_value->value);
  *remainder -= strlen(tag_value->value);
  if(*remainder <= 0){
    printf("QIO_encode_as_string: Buffer overflow\n");
    return buf;
  }
  buf = QIO_write_endtag(buf, tag_value->tag, remainder);
  return buf;
}

#define QIO_MAXINTSTRING 16
char *QIO_encode_as_int(char *buf, QIO_TagIntValue *tag_value,
			  int *remainder){
  char int_string[QIO_MAXINTSTRING];

  /* Don't write value unless occurs */
  if(!tag_value->occur)return buf;
  buf = QIO_write_tag(buf, tag_value->tag, tag_value->attr, remainder);
  snprintf(int_string,QIO_MAXINTSTRING,"%d",tag_value->value);
  *remainder -= strlen(int_string);
  if(*remainder <= 0){
    printf("QIO_encode_as_int: Buffer overflow\n");
    return buf;
  }
  buf = strcat(buf, int_string);
  buf = QIO_write_endtag(buf, tag_value->tag, remainder);
  return buf;
}

char *QIO_encode_as_hex32(char *buf, QIO_TagHex32Value *tag_value,
			  int *remainder){
  char int_string[QIO_MAXINTSTRING];

  /* Don't write value unless occurs */
  if(!tag_value->occur)return buf;
  buf = QIO_write_tag(buf, tag_value->tag, tag_value->attr, remainder);
  snprintf(int_string,QIO_MAXINTSTRING,"%x",tag_value->value);
  *remainder -= strlen(int_string);
  if(*remainder <= 0){
    printf("QIO_encode_as_hex32: Buffer overflow\n");
    return buf;
  }
  buf = strcat(buf, int_string);
  buf = QIO_write_endtag(buf, tag_value->tag, remainder);
  return buf;
}

char *QIO_encode_as_intlist(char *buf, 
			      QIO_TagIntArrayValue *tag_value, 
			      int n, int *remainder){
  int i;
  char int_string[QIO_MAXINTSTRING];
  
  /* Don't write value unless occurs */
  if(!tag_value->occur)return buf;
  buf = QIO_write_tag(buf, tag_value->tag, tag_value->attr, remainder);
  if(*remainder <= 0){
    printf("QIO_encode_as_intlist: Buffer overflow\n");
    return buf;
  }
  
  for(i = 0; i < n ; i++){
    snprintf(int_string, QIO_MAXINTSTRING, "%d ", tag_value->value[i]);
    *remainder -= strlen(int_string);
    if(*remainder <= 0){
      printf("QIO_encode_as_intlist: Buffer overflow\n");
      return buf;
    }
    buf = strcat(buf, int_string);
  }
  buf = QIO_write_endtag(buf, tag_value->tag, remainder);
  return buf;
}

int QIO_check_string_occur(QIO_TagCharValue *tag_value){
  if(tag_value->occur != 1){
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("QIO_check_string_occur: error %s tag count %d\n",
	     tag_value->tag, tag_value->occur);
    return 1;
  }
  return 0;
}

int QIO_check_int_occur(QIO_TagIntValue *tag_value){
  if(tag_value->occur != 1){
    printf("QIO_check_int_occur: error %s tag count %d\n",
	   tag_value->tag, tag_value->occur);
    return 1;
  }
  return 0;
}

int QIO_check_intarray_occur(QIO_TagIntArrayValue *tag_value){
  if(tag_value->occur != 1){
    printf("QIO_check_intarray_occur: error %s tag count %d\n",
	   tag_value->tag, tag_value->occur);
    return 1;
  }
  return 0;
}

int QIO_check_hex32_occur(QIO_TagHex32Value *tag_value){
  if(tag_value->occur != 1){
    printf("QIO_check_hex32_occur: error %s tag count %d\n",
	   tag_value->tag, tag_value->occur);
    return 1;
  }
  return 0;
}

