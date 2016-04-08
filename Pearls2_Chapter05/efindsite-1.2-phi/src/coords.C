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


#include "coords.h"

using namespace std;

// ==================================================================================   CoordsProtein

CoordsProtein::CoordsProtein( int an, int ar, double ax, double ay, double az, string at, string aa )
{
 n = an;
 r = ar;
 x = ax;
 y = ay;
 z = az;
 t = at;
 a = aa;
}

CoordsProtein::CoordsProtein( void )
{
 n = 0;
 r = 0;
 x = 0.0;
 y = 0.0;
 z = 0.0;
 t = "";
 a = "";
}

CoordsProtein::~CoordsProtein() {}

int CoordsProtein::getAtomNumber( void )
{
 return n;
}

int CoordsProtein::getResidueNumber( void )
{
 return r;
}

double CoordsProtein::getCoords( int an )
{
 switch (an)
 {
  case  1  :  return x;
  case  2  :  return y;
  case  3  :  return z;

  default  :  return 0;
 }
}

void CoordsProtein::setCoords( double ax, double ay, double az)
{
 x = ax;
 y = ay;
 z = az;
}

string CoordsProtein::getResidueName( void )
{
 return t;
}

string CoordsProtein::getAtomName( void )
{
 return a;
}


// ==================================================================================   CoordsLigand

CoordsLigand::CoordsLigand( int an, double ax, double ay, double az, string aa )
{
 n = an;
 x = ax;
 y = ay;
 z = az;
 a = aa;
}

CoordsLigand::CoordsLigand( void )
{
 n = 0;
 x = 0.0;
 y = 0.0;
 z = 0.0;
 a = "";
}

CoordsLigand::~CoordsLigand() {}

int CoordsLigand::getAtomNumber( void )
{
 return n;
}

double CoordsLigand::getCoords( int an )
{
 switch (an)
 {
  case  1  :  return x;
  case  2  :  return y;
  case  3  :  return z;

  default  :  return 0;
 }
}

void CoordsLigand::setCoords( double ax, double ay, double az)
{
 x = ax;
 y = ay;
 z = az;
}

string CoordsLigand::getAtomName( void )
{
 return a;
}
