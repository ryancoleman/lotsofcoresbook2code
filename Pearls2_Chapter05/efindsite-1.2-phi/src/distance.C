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


#include "distance.h"

using namespace std;

double getDistance( int ms1, Template * tpl1, int lig1, Template * tpl2, int lig2 )
{
 double tcen1[3];
 double tcen2[3];
 
 double tt[3];
 double uu[3][3];
 
 double dst1 = 0.0;
 
 bool ret1 = tpl1->getMatrix( tpl2->getProteinID(), tt, uu );
 
 if ( !ret1 )
 {
  tpl1->getLigandCenter(lig1, tcen1, false);
  tpl2->getLigandCenter(lig2, tcen2, false);
  
  double tcen3[3];
  
  for ( int i1 = 0; i1 < 3; i1++ )
   tcen3[i1] = tt[i1] + uu[i1][0] * tcen2[0] + uu[i1][1] * tcen2[1] + uu[i1][2] * tcen2[2];
  
  for ( int i1 = 0; i1 < 3; i1++ )
   tcen2[i1] = tcen3[i1];
  
  dst1 = sqrt(pow( tcen1[0] - tcen2[0], 2.0) + pow( tcen1[1] - tcen2[1], 2.0) + pow( tcen1[2] - tcen2[2], 2.0));
 }
 
 tpl1->getLigandCenter(lig1, tcen1, true);
 tpl2->getLigandCenter(lig2, tcen2, true);
 
 double dst2 = sqrt(pow( tcen1[0] - tcen2[0], 2.0) + pow( tcen1[1] - tcen2[1], 2.0) + pow( tcen1[2] - tcen2[2], 2.0));
 
 if ( dst2 > dst1 )
  dst1 = dst2;
 
 return dst1;
}
