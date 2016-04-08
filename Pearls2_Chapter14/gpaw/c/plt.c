/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "extensions.h"

int write_plt_file(char *fname,
		   int nx, int ny, int nz,
		   double x0, double y0, double z0,
		   double dx, double dy, double dz,
		   double *grid);

/* write grid to binary plt (gOpenMol) plot file */
PyObject* WritePLT(PyObject *self, PyObject *args)
{
  char* fname;       /* file name     */
  PyArrayObject* ho; /* grid spacings */
  PyArrayObject* go; /* grid to write */
  if (!PyArg_ParseTuple(args, "sOO", &fname, &ho, &go)) 
    return NULL; 

  /* must be 3D */
  if(PyArray_NDIM(go) != 3) return NULL; 

  double* g = DOUBLEP(go);
  double* h = DOUBLEP(ho);

  write_plt_file(fname,
		 PyArray_DIM(go, 0),
		 PyArray_DIM(go, 1),
		 PyArray_DIM(go, 2),
		 0.,0.,0.,
		 h[0],h[1],h[2],
		 g);
  Py_RETURN_NONE;
}

/* -----------------------------------------------------------------
 * write grid to binary plt (gOpenMol) plot file
 *
 * x0, dx etc are assumed to be atomic units
 * the grid is assumed to be in the format:
 * grid(ix,iy,iz) = grid[ ix + ( iy + iz*ny )*nx ];
 * where ix=0..nx-1 etc
 */

/* stolen from pltfile.c */
#define FWRITE(value , size)  { \
Items = fwrite(&value, size , 1L , Output_p);\
if(Items < 1) {\
  printf("?ERROR - in writing contour file (*)\n");\
  return(1);}}

int write_plt_file(char *fname,
		   int nx, int ny, int nz,
		   double x0, double y0, double z0,
		   double dx, double dy, double dz,
		   double *grid) {
  FILE *Output_p;
  static int Items;
  float scale,zmin,zmax,ymin,ymax,xmin,xmax,val;
  int rank,TypeOfSurface;
  int ix,iy,iz,indx;
  double norm,sum,dV;

  Output_p = fopen(fname,"wb");
  
  /* see http://www.csc.fi/gopenmol/developers/plt_format.phtml */

#define au_A 0.52917725
  scale = au_A; /* atomic length in Angstroem */

  rank=3; /* always 3 */
  FWRITE(rank , sizeof(int));
  TypeOfSurface=4; /* arbitrary */
  FWRITE(TypeOfSurface , sizeof(int));
  FWRITE(nz , sizeof(int));
  FWRITE(ny , sizeof(int));
  FWRITE(nx , sizeof(int));
  zmin= scale * ((float) z0);
  zmax= scale * ((float) z0+(nz-1)*dz);
  /* float zmax=(float) z0+nz*dz; */
  FWRITE(zmin , sizeof(float));
  FWRITE(zmax , sizeof(float));
  ymin= scale * ((float) y0);
  ymax= scale * ((float) y0+(ny-1)*dy);
  /* float ymax=(float) y0+ny*dy; */
  FWRITE(ymin , sizeof(float));
  FWRITE(ymax , sizeof(float));
  xmin= scale * ((float) x0);
  xmax= scale * ((float) x0+(nx-1)*dx);
  /* float xmax=(float) x0+nx*dx; */
  FWRITE(xmin , sizeof(float));
  FWRITE(xmax , sizeof(float));
      
  indx=0;
  norm = 0;
  sum=0;
  dV=dx*dy*dz;
  for(iz=0;iz<nz;iz++)
    for(iy=0;iy<ny;iy++)
      for(ix=0;ix<nx;ix++) {
	val = (float) grid[indx];
	sum += val;
	norm += val*val;
	FWRITE(val , sizeof(float));
	indx++;
      }

  fclose(Output_p);

  printf("#<write_plt_file> %s written (sum=%g,norm=%g)\n",
	 fname,sum*dV,norm*dV);
  
  return 0;
}

