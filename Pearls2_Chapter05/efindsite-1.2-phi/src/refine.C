/*
===============================================================================
         ___________.__            .____________.__  __          
     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
        \/    \/            \/      \/       \/               \/ 

                                                  
   eFindSite - ligand-binding site prediction from meta-threading

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org wfeinstein@lsu.edu

   Copyright 2013 Michal Brylinski, Wei Feinstein

   This file is part of eFindSite.

   eFindSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eFindSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#include "refine.h"
#include "size.h"

using namespace std;

int refine_pockets( multimap< int, Template *, greater<int> > &template_dat, int ini1, int ini2[], int ini3, double ini4 )
{
 int xnl1 = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator xpl1;
 
 for ( xpl1 = template_dat.begin(); xpl1 != template_dat.end(); xpl1++ )
  for ( int xil1 = 0; xil1 < ((*xpl1).second)->getLigandsTotal(); xil1++ )
   xnl1++;
 
 double *ref1 = new double [xnl1*3];
 
 int xnl2 = 0;
 
 for ( xpl1 = template_dat.begin(); xpl1 != template_dat.end(); xpl1++ )
  for ( int xil1 = 0; xil1 < ((*xpl1).second)->getLigandsTotal(); xil1++ )
  {
   double xcen1[3];
   
   ((*xpl1).second)->getLigandCenter(xil1, xcen1, true);
   
   for ( int xil2 = 0; xil2 < 3; xil2++ )
    ref1[xnl2*3+xil2] = xcen1[xil2];
   
   xnl2++;
  }
 
 double *ref2 = new double [ini1*3];
 
 for ( int xil3 = 0; xil3 < ini1; xil3++ )
 {
  double xcen2[4];
  
  for ( int xil4 = 0; xil4 < 4; xil4++ )
   xcen2[xil4] = 0.0;
  
  for ( int xil5 = 0; xil5 < xnl2; xil5++ )
   if ( xil3 == ini2[xil5] )
   {
    for ( int xil6 = 0; xil6 < 3; xil6++ )
     xcen2[xil6] += ref1[xil5*3+xil6];
    
    xcen2[3] += 1.0;
   }
  
  for ( int xil7 = 0; xil7 < 3; xil7++ )
   ref2[xil3*3+xil7] = xcen2[xil7] / xcen2[3];
 }
 
 delete [] ref1;
 
 double *ref3 = new double [ini1*ini1];
 
 for ( int xil8 = 0; xil8 < ini1; xil8++ )
  for ( int xil9 = 0; xil9 < ini1; xil9++ )
   ref3[xil8*ini1+xil9] = sqrt(pow( ref2[xil8*3+0] - ref2[xil9*3+0], 2.0) + pow( ref2[xil8*3+1] - ref2[xil9*3+1], 2.0) + pow( ref2[xil8*3+2] - ref2[xil9*3+2], 2.0));
 
 delete [] ref2;
 
 int * clx1 = new int [ini1];
 
 int clx2 = cluster_avelink( ref3, clx1, ini1, ini4, "min" );
 
 delete [] ref3;
 
 for ( int xil10 = 0; xil10 < ini3; xil10++ )
  ini2[xil10] = clx1[ini2[xil10]];
 
 delete [] clx1;
 
 return clx2;
}
