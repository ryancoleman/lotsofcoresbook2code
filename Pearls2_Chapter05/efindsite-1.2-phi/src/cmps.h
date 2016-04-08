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


#ifndef __CMPS_H_
#define __CMPS_H_

#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<bitset>
#include<cstdlib>

#include "size.h"

using namespace std;

class Cmps {
        
  private:
    
    int _nc;                              // number of auxiliary compounds
    
    bitset<MAXSMI> _cmps_fpt_smi[MAXCMP]; // SMILES fingerprints
    bitset<MAXMAC> _cmps_fpt_mac[MAXCMP]; // MACCS fingerprints
    
    double         _cmps_mw[MAXCMP];      // molecular weight
    double         _cmps_logp[MAXCMP];    // logp
    double         _cmps_psa[MAXCMP];     // polar surface area
    
    int            _cmps_hbd[MAXCMP];     // HB donors
    int            _cmps_hba[MAXCMP];     // HB acceptors
    
  public:
    
    Cmps( int );
    
    Cmps( void );
    
    ~Cmps();
    
    bool loadCompounds( std::string );
    
    int getCmpsTotal( void );
    
    void getSMILES( int, bitset<MAXSMI> & );
    void getMACCS( int, bitset<MAXMAC> & );
    
    double getMW( int );
    double getLOGP( int );
    double getPSA( int );
    
    int getHBD( int );
    int getHBA( int );
};

#endif
