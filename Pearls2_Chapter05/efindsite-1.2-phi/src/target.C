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


#include "target.h"

using namespace std;

Target::Target( int aa, int ar )
{
 _na      = aa;
 _nr      = ar;
}

Target::Target( void )
{
 _na      = 0;
 _nr      = 0;
}

Target::~Target() {}


// ==================================================================================   loadTarget

bool Target::loadTarget( std::string p1_name )
{
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  if ( line1.size() > 53 )
   if ( line1.substr(0,6) == "ATOM  ")
   {
    string residue1 = line1.substr(17,3);
    string atom1    = line1.substr(12,4);
    
    int residue2 = atoi(line1.substr(22,4).c_str());
    int atom2    = atoi(line1.substr(6,5).c_str());
    
    double x1 = atof(line1.substr(30,8).c_str());
    double y1 = atof(line1.substr(38,8).c_str());
    double z1 = atof(line1.substr(46,8).c_str());
    
    _protein_xyz.push_back( CoordsProtein( atom2, residue2, x1, y1, z1, residue1, atom1 ) );
    
    _na++;
    
    if ( atom1 == " CA " )
    {
     _protein_seq1.append( three2oneS( line1.substr(17,3) ) );
     
     _protein_seq3[_nr] = residue2;
     
     _nr++;
    }
   }
 
 p1_file.close();
 
 strcpy(_protein_seq2, _protein_seq1.c_str());
 
 return EXIT_SUCCESS;
}


// ==================================================================================   loadPsipred

bool Target::loadPsipred( std::string p1_name )
{
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  if ( line1.size() > 29 )
  {
   int res1 = atoi(line1.substr(0,4).c_str());
   
   _protein_prf1[res1-1][0] = atof(line1.substr(9,7).c_str());
   _protein_prf1[res1-1][1] = atof(line1.substr(16,7).c_str());
   _protein_prf1[res1-1][2] = atof(line1.substr(23,7).c_str());
  }
 
 p1_file.close();
 
 return EXIT_SUCCESS;
}


// ==================================================================================   loadSequence

bool Target::loadSequence( std::string p1_name )
{
 for ( int ii1 = 0; ii1 < _nr; ii1++ )
  for ( int ii2 = 0; ii2 < 20; ii2++ )
   _protein_prf1[ii1][ii2] = 0.0;
 
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() ) { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 const char * line4 = 0;
 
 int res1 = 0;
 int res2 = -2;
 
 while (getline(p1_file,line1))
 {
  char * line2, * line3;
  
  line2 = new char [line1.size()+1];
  
  strcpy (line2, line1.c_str());
  
  line3 = strtok (line2," ");
  
  while ( line3 != NULL )
  {
   if ( isalpha(line3[0]) )
   {
    res1 = atoi(line4);
    
    res2 = -2;
   }
   
   if ( res1 > 0 && res2 > -1 && res2 < 20  )
    _protein_prf1[res1-1][res2] = atof(line3);
   
   line4 = strdup(line3);
   
   res2++;
   
   line3 = strtok(NULL," ");
  }

  delete[] line2;
  delete[] line3;
 }
 
 p1_file.close();
 
 return EXIT_SUCCESS;
}


// ==================================================================================   getProteinSequence

std::string Target::getProteinSequence( void )
{
 return _protein_seq1;
}


// ==================================================================================   getProteinResiduesTotal

int Target::getProteinResiduesTotal( void )
{
 return _nr;
}


// ==================================================================================   getProteinCoordsCA

int Target::getProteinCoordsCA( double tabCA1[][3] )
{
 vector<CoordsProtein>::iterator ipca1;
 
 int ipca2 = 0;
 
 for ( ipca1 = _protein_xyz.begin(); ipca1 < _protein_xyz.end(); ipca1++ )
  if ( (*ipca1).getAtomName() == " CA " )
  {
   for ( int ipca3 = 0; ipca3 < 3; ipca3++ )
    tabCA1[ipca2][ipca3] = (*ipca1).getCoords(ipca3+1);
   
   ipca2++;
  }
 
 return ipca2;
}


// ==================================================================================   getProteinCoords1D

void Target::getProteinCoords1D( double tabCA1[] )
{
 vector<CoordsProtein>::iterator ipca1;
 
 int ipca2 = 0;
 
 for ( ipca1 = _protein_xyz.begin(); ipca1 < _protein_xyz.end(); ipca1++ )
  if ( (*ipca1).getAtomName() == " CA " )
   for ( int ipca3 = 0; ipca3 < 3; ipca3++ )
    tabCA1[ipca2++] = (*ipca1).getCoords(ipca3+1);
}


// ==================================================================================   compositionSVM

bool Target::compositionSVM( ModelSVM * fsvm1 )
{
 for ( int st1 = 0; st1 < _nr; st1++ )
 {
  double fsvm2[MAXSV1];
  
  for ( int st2 = 0; st2 < MAXSV1; st2++ )
   fsvm2[st2] = _protein_prf1[st1][st2];
  
  _protein_svm1[st1] = fsvm1->SVMpredict( 1, fsvm2 );
 }
 
 return EXIT_SUCCESS;
}


// ==================================================================================   getCompositionScore

double Target::getCompositionScoreSVM( int res1 )
{
 return _protein_svm1[res1];
}


// ==================================================================================   getProteinNumbering

int Target::getProteinNumbering ( int mapp1 )
{
 return _protein_seq3[mapp1];
}
