/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#ifndef _XC_GPAW_H
#define _XC_GPAW_H

/*
BETA = 0.066725
MU = BETA * pi * pi / 3
C2 = (1 / (18 * pi)**(1 / 3))
C0I = 3 / (4 * pi)
C1 = -9 / (8 * pi) * (2 * pi / 3)**(1 / 3)
CC1 = 1 / (2**(4 / 3) - 2)
CC2 = 4 * CC1 / 3
IF2 = 3 / (2 * CC2);
C3 = pi * (4 / (9 * pi))**(1 / 3) / 16
C0 = 4 * pi / 3
*/

#define BETA   0.066725
#define GAMMA  0.031091
#define MU     0.2195164512208958
#define C2     0.26053088059892404
#define C0I    0.238732414637843
#define C1    -0.45816529328314287
#define CC1    1.9236610509315362
#define CC2    2.5648814012420482
#define IF2    0.58482236226346462
#define C3     0.10231023756535741
#define C0     4.1887902047863905
#define THIRD  0.33333333333333333
#define NMIN   1.0E-10

typedef int bool;

typedef struct
{
  bool gga;
  double kappa;
  int nparameters;
  double parameters[110];
} xc_parameters;

#endif /* _XC_GPAW_H */
