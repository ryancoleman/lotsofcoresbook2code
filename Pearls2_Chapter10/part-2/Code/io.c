
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"

static double kmax_cut;
static double tau0;

int load_txt_dbl(char* filename, int columns, double* values, int* size){
	FILE* fp;
	char line[MAXLEN];
	int n;
	long int i = -1;
	char** cptr = (char**) malloc(sizeof(char*));

	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
		return 0;
	}
	
	while (fgets(line, MAXLEN, fp)) {
		if (*line == '#') continue;
		
		values[++i] = strtod( line, cptr);

		for(n=1;n<columns;n++){
			if (cptr!=NULL ) values[++i] = strtod( *cptr, cptr);
		}
		
		if ( i>(columns*(MAXLINES-1)) ) break;
	}
	*size = (i/columns)+1;
// 	*size = i;
	return 1;
}

int load_two(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = d2 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
  
    if (*line == '#') continue;	
    d1 = d2 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    if ( i>(2*MAXLINES-2) ) break;
  }
  *size = i/2 + 1;
  
  return 1;
}
int load_three_int(char* filename, int* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    values[++i] = (int)d1;
    values[++i] = (int)d2;
    values[++i] = (int)d3;
    if ( i>(3*MAXLINES-3) ) break;
  }
  *size = i/3+1;

  return 1;
}

int load_three(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    values[++i] = d3;
    if ( i>(3*MAXLINES-3) ) break;
  }
  *size = i/3+1;

  return 1;
}

double get_kmax_cut(){
	return kmax_cut;
}

double get_tau0(){
	return tau0;
}
