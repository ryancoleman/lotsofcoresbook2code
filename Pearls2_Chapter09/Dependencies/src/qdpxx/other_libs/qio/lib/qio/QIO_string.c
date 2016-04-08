/* Safe strings for QIO */

#include <qio_config.h>
#include <qio_string.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif

/*
 *  QIO_String manipulation utilities
 */

/* String creation */
QIO_String *QIO_string_create(void)
{
  QIO_String *qs;

  qs = (QIO_String *)malloc(sizeof(QIO_String));
  if(qs == NULL) return NULL;
  qs->string = NULL;
  qs->length = 0;
  return qs;
}

void QIO_string_destroy(QIO_String *qs)
{
  if(qs != NULL) {
    if(qs->length>0) { free(qs->string); qs->string=NULL; }
    free(qs);
  }
}

/* set String */
void QIO_string_set(QIO_String *qs, const char *const string)
{
  char myname[] = "QIO_string_set";
  if(qs != NULL) {
 
    /* If a NULL string is passed in set qs to be its default
       state */

    if(string==NULL) {
      if(qs->length>0 && (qs->string != NULL)) free(qs->string);
      qs->string = NULL;
      qs->length = 0;
    } else {

      /* string is not NULL, copy it in */
      size_t len = strlen(string) + 1;
#if 1
      /* Do explicit realloc ing */
      if( qs->string!=NULL ) { free(qs->string); }
      qs->string = (char *)malloc(len);
      if( qs->string == NULL ) {
	fprintf(stderr, "%s: Oops Malloc failed\n",myname);
	exit(-1);
      }
#else
      qs->string = (char *)realloc(qs->string, len);
#endif
      memcpy(qs->string, string, len);
      qs->length = len;
    }
  }
  else { 
    fprintf(stderr,"%s: Attempt to set NULL QIO_String\n",myname);
    fflush(stderr);
    exit(-1);
  }
}

/* Size of string */
size_t QIO_string_length(const QIO_String *const qs)
{
  char myname[] = "QIO_string_length";
  if( qs != NULL ) {
    return qs->length;
  }
  else { 
    fprintf(stderr, "%s: Attempt to get length of NULL QIO_String*\n",myname);
    fflush(stderr);
    exit(-1);
  }
}

/* Return pointer to string data */
char *QIO_string_ptr(QIO_String *qs)
{
  if( qs != NULL ) { 
    return qs->string;
  }
  else { 
    fprintf(stderr, "QIO_string_ptr: Attempt to get length of NULL QIO_String*\n");
    fflush(stderr);
    exit(-1);
  }
}

void QIO_string_copy(QIO_String *dest, QIO_String *src)
{
  if(src == NULL)return;
  if(dest == NULL)return;
  if(src == dest ) return;

  if( src->length > 0 ) { 
    size_t len = src->length;
#if 1
    if( dest->string ) { free(dest->string); }
    dest->string = (char *)malloc(len);
    if(dest->string == NULL ) {
      fprintf(stderr, "QIO_string_copy: Bugger\n");
      fflush(stderr);
      exit(-1);
    }
#else
    dest->string = (char *)realloc(dest->string, len);
#endif
    memcpy(dest->string, src->string, len);
    dest->length = len;
  }
  else { 
    if( dest->string != NULL ) { 
      free(dest->string); 
      dest->string=NULL;
    }
    dest->length=0;
  }
    
}

/* Non-destructive string reallocation */
void QIO_string_realloc(QIO_String *qs, int length)
{
  int i,min;
  char *tmp;

  /* If NULL string is passed just return */
  if(qs == NULL) return;

  /* If length == 0, freeing the string is "non destructive" */
  if( length == 0 ) { 
    qs->length = 0;
    if( qs->string != NULL ) { free(qs->string); qs->string=NULL; return; }
  }


  tmp = (char *)malloc(length);
  if(tmp == NULL) {
    /* Returning here is a bad idea as the string is not really
       realloced */
    printf("QIO_string_realloc: Can't malloc size %d\n",length);
    fflush(stdout);
    exit(-1);
    /* return */
  }
  
  for(i = 0; i < length; i++) tmp[i] = '\0';

  if(qs->string != NULL)
    {
      /* Follow semantics of realloc - shrink or grow/copy */
      min = length > qs->length ? qs->length : length;
      strncpy(tmp,qs->string,min);
      tmp[min-1] = '\0';  /* Assure null termination */
      free(qs->string);
    }      

  qs->length = length;
  qs->string = tmp;
}

void QIO_string_append(QIO_String *qs, const char *const string){
  int len, totlen;

  if(qs == NULL || string == NULL)return;
  len = strlen(string);
  if(len == 0)return;

  /* If qs string is uninitialized, set it to blank */
  if(qs->length == 0)
    QIO_string_set(qs, "");

  /* Reallocate to the desired total length if necessary */
  totlen = qs->length + len;
  if(totlen > qs->length)
    QIO_string_realloc(qs, totlen);

  /* Finally, append */
  strncat(qs->string, string, len + 1);
}

