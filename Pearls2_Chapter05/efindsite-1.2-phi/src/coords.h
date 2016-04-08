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


#ifndef __COORDS_H_
#define __COORDS_H_

#include<string>
#include<deque>

using namespace std;

// ==================================================================================   CoordsProtein

class CoordsProtein {
        
  private:
    
    int    n; // atom number
    int    r; // residue number
    double x; // x coord
    double y; // y coord
    double z; // z coord
    string t; // residue name
    string a; // atom name

  public:

    CoordsProtein( int, int, double, double, double, string, string );
    
    CoordsProtein( void );

    ~CoordsProtein();

    int getAtomNumber();
    
    int getResidueNumber();
    
    double getCoords( int );
    
    void setCoords( double, double, double );
    
    string getResidueName();
    
    string getAtomName();
};

// ==================================================================================   CoordsLigand

class CoordsLigand {
        
  private:
    
    int    n; // atom number
    double x; // x coordinate
    double y; // y coordinate
    double z; // z coordinate
    string a; // atom name

  public:

    CoordsLigand( int, double, double, double, string );
    
    CoordsLigand( void );

    ~CoordsLigand();

    int getAtomNumber();
    
    double getCoords( int );
    
    void setCoords( double, double, double );
    
    string getAtomName();
};

#endif
