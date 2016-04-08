/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#ifndef LFC_H
#define LFC_H

#include <Python.h>

typedef struct
{
  const double* A_gm;  // function values
  int nm;              // number of functions (2*l+1)
  int M;               // global number of first function
  int W;               // volume number
} LFVolume;


typedef struct 
{
  PyObject_HEAD
  double dv;                 // volume per grid point
  int nW;                    // number of volumes
  int nB;                    // number of boundary points
  double* work_gm;           // work space
  LFVolume* volume_W;        // pointers to volumes
  LFVolume* volume_i;        // pointers to volumes at current grid point
  int* G_B;                  // boundary grid points
  int* W_B;                  // volume numbers
  int* i_W;                  // mapping from all volumes to current volumes
  int* ngm_W;                // number of grid points per volume
  bool bloch_boundary_conditions;  // Gamma-point calculation?
  complex double* phase_kW;  // phase factors: exp(ik.R)
  complex double* phase_i;   // phase factors for current volumes
} LFCObject;


#define GRID_LOOP_START(lfc, k)                                    \
{                                                                  \
  int* G_B = lfc->G_B;                                             \
  int* W_B = lfc->W_B;                                             \
  int* i_W = lfc->i_W;                                             \
  complex double* phase_i = lfc->phase_i;                          \
  LFVolume* volume_i = lfc->volume_i;                              \
  LFVolume* volume_W = lfc->volume_W;                              \
  double complex* phase_W = lfc->phase_kW + k * lfc->nW;           \
  int Ga = 0;                                                      \
  int ni = 0;                                                      \
  for (int B = 0; B < lfc->nB; B++)                                \
    {                                                              \
      int Gb = G_B[B];                                             \
      int nG = Gb - Ga;                                            \
      if (nG > 0)                                                  \
        {

#define GRID_LOOP_STOP(lfc, k)                                     \
          for (int i = 0; i < ni; i++)                             \
            volume_i[i].A_gm += nG * volume_i[i].nm;               \
        }                                                          \
      int Wnew = W_B[B];                                           \
      if (Wnew >= 0)                                               \
        {                                                          \
          /* Entering new sphere: */                               \
          volume_i[ni] = volume_W[Wnew];                           \
          if (k >= 0)                                              \
            phase_i[ni] = phase_W[Wnew];                           \
          i_W[Wnew] = ni;                                          \
          ni++;                                                    \
        }                                                          \
      else                                                         \
        {                                                          \
          /* Leaving sphere: */                                    \
          int Wold = -1 - Wnew;                                    \
          int iold = i_W[Wold];                                    \
          volume_W[Wold].A_gm = volume_i[iold].A_gm;               \
          ni--;                                                    \
          volume_i[iold] = volume_i[ni];                           \
          if (k >= 0)                                              \
            phase_i[iold] = phase_i[ni];                           \
          int Wlast = volume_i[iold].W;                            \
          i_W[Wlast] = iold;                                       \
        }                                                          \
      Ga = Gb;                                                     \
    }                                                              \
  for (int W = 0; W < lfc->nW; W++)                                \
    volume_W[W].A_gm -= lfc->ngm_W[W];                             \
}

#endif
