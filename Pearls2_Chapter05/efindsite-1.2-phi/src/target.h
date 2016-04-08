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


#ifndef __TARGET_H_
#define __TARGET_H_

#include<vector>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cstdio>
#include<cctype>
#include<string>

#include "coords.h"
#include "data.h"
#include "runsvm.h"
#include "size.h"

using namespace std;

class Target {
        
  private:
    
    vector<CoordsProtein> _protein_xyz;                  // protein coords
    
    int    _na;                           // number of protein atoms
    int    _nr;                           // number of protein residues
    
    string _protein_seq1;                 // aa sequence
    char   _protein_seq2[MAXPRO];         // aa sequence
    int    _protein_seq3[MAXPRO];         // aa sequence numbering
    double _protein_prf1[MAXPRO][MAXSV1]; // sequence profiles
    double _protein_svm1[MAXPRO];         // composition SVM
    
  public:
    
    Target( int, int );
    
    Target( void );
    
    ~Target();
    
    bool loadTarget( std::string );
    
    bool loadPsipred( std::string );
    
    bool loadSequence( std::string );
    
    std::string getProteinSequence( void );
    
    int getProteinResiduesTotal( void );
    
    int getProteinCoordsCA( double [][3] );
    
    void getProteinCoords1D( double [] );
    
    bool compositionSVM( ModelSVM * );
    
    double getCompositionScoreSVM( int );
    
    int getProteinNumbering( int );
};

#endif
