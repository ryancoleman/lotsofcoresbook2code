#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/* Write out array */
void array_write(int *n, char filename[100], double grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "w" );
   if (fp == NULL) {
   		printf("**** Error opening file: %s. ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fwrite(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);

   return;
}

/* Read in array */
void array_read(int *n, char filename[100], double grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "r");
   if (fp == NULL) {
   		printf("**** Error opening file: %s ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fread(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);

   return;
}

/* Double vector */
_OFFLOADABLE
double *create_vector(int length){
  return (double*) _mm_malloc( (size_t) length*sizeof(double), 64);
}

/* Double integer vector */
int *create_ivector(int length){
   return (int*) malloc( (size_t) length*sizeof(int));
}

void ivector_read(int *n, char filename[100], int grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "r");
   if (fp == NULL) {
   		printf("**** Error opening file: %s ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fread(grid, (size_t) sizeof(int), (size_t) (*n), fp);
   fclose(fp);

   return;
}

/* Double 2-D array */

double **create_array(int dim_x, int dim_y){
   int i; double **array;
   
   long int size1;
   long int size2;

   size1 = dim_x;
   size2 = size1*dim_y;

   array = (double**) malloc(size1 * sizeof(double*));
   array[0] = (double*) malloc(size2 * sizeof(double));

   for(i=1;i< dim_x;i++)    /* loop over rows */
     array[i]=array[0] + i * dim_y;

   return array;
}

void destroy_array(double **array){
   free(array[0]);
   free(array);
}
