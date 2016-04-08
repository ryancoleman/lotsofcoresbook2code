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


#include "cmps.h"

using namespace std;

Cmps::Cmps( int ac )
{
 _nc = ac;
}

Cmps::Cmps( void )
{
 _nc = 0;
}

Cmps::~Cmps() {}


// ==================================================================================   loadCompounds

bool Cmps::loadCompounds( std::string c1_name )
{
 string line1;
 
 std::vector<string> tmp_sdf_t;
 
 ifstream c1_file( c1_name.c_str() );
 
 if ( !c1_file.is_open() ) { cout << "Cannot open " << c1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(c1_file,line1))
 {
  tmp_sdf_t.push_back(line1);
  
  if ( line1.compare("$$$$") == 0 )
  {
   if ( _nc < MAXCMP )
   {
    int n1a = atoi(tmp_sdf_t[3].substr(0,3).c_str());
    int n1b = atoi(tmp_sdf_t[3].substr(3,3).c_str());
    
    bool w1 = false;
    bool w2 = false;
    bool w3 = false;
    bool w4 = false;
    bool w5 = false;
    bool w6 = false;
    bool w7 = false;
    
    for ( int i1 = 4 + n1a + n1b; i1 < (int) tmp_sdf_t.size() - 1; i1++ )
    {
     if ( tmp_sdf_t[i1].find("FINGERPRINT") != string::npos )
     {
      string::iterator i2;
      int i3;
      
      for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
       switch(*i2)
       {
        case '0': _cmps_fpt_smi[_nc][MAXSMI-1-i3] = 0; break;
        case '1': _cmps_fpt_smi[_nc][MAXSMI-1-i3] = 1; break;
       }
      
      w1 = true;
     }
     
     else if ( tmp_sdf_t[i1].find("MACCS166") != string::npos )
     {
      string::iterator i2;
      int i3;
      
      for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
       switch(*i2)
       {
        case '0': _cmps_fpt_mac[_nc][MAXMAC-1-i3] = 0; break;
        case '1': _cmps_fpt_mac[_nc][MAXMAC-1-i3] = 1; break;
       }
      
      w2 = true;
     }
     
     else if ( tmp_sdf_t[i1].find("OB_MW") != string::npos )
     {
      _cmps_mw[_nc] = atof(tmp_sdf_t[i1+1].c_str());
      
      w3 = true;
     }
     
     else if ( tmp_sdf_t[i1].find("OB_logP") != string::npos )
     {
      _cmps_logp[_nc] = atof(tmp_sdf_t[i1+1].c_str());
      
      w4 = true;
     }
     
     else if ( tmp_sdf_t[i1].find("OB_PSA") != string::npos )
     {
      _cmps_psa[_nc] = atof(tmp_sdf_t[i1+1].c_str());
      
      w5 = true;
     }
     
     else if ( tmp_sdf_t[i1].find("MCT_HBD") != string::npos )
     {
      _cmps_hbd[_nc] = atoi(tmp_sdf_t[i1+1].c_str());
      
      w6 = true;
     }
      
     else if ( tmp_sdf_t[i1].find("MCT_HBA") != string::npos )
     {
      _cmps_hba[_nc] = atoi(tmp_sdf_t[i1+1].c_str());
      
      w7 = true;
     }
    }
    
    if ( w1 && w2 && w3 && w4 && w5 && w6 && w7 )
     _nc++;
   }
   
   tmp_sdf_t.clear();
  }
 }
 
 c1_file.close();
 
 return EXIT_SUCCESS;
}


// ==================================================================================   getCmpsTotal

int Cmps::getCmpsTotal( void )
{
 return _nc;
}


// ==================================================================================   getSMILES

void Cmps::getSMILES( int ac, bitset<MAXSMI> &af )
{
 af = _cmps_fpt_smi[ac];
}


// ==================================================================================   getMACCS

void Cmps::getMACCS( int ac, bitset<MAXMAC> &af )
{
 af = _cmps_fpt_mac[ac];
}


// ==================================================================================   getMW

double Cmps::getMW( int ac )
{
 return _cmps_mw[ac];
}


// ==================================================================================   getLOGP

double Cmps::getLOGP( int ac )
{
 return _cmps_logp[ac];
}


// ==================================================================================   getPSA

double Cmps::getPSA( int ac )
{
 return _cmps_psa[ac];
}


// ==================================================================================   getHBD

int Cmps::getHBD( int ac )
{
 return _cmps_hbd[ac];
}


// ==================================================================================   getHBA

int Cmps::getHBA( int ac )
{
 return _cmps_hba[ac];
}
